#ifndef ERF_FIRE_ATM_STATE_H
#define ERF_FIRE_ATM_STATE_H

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_REAL.H>
#include <ERF_MicrophysicsUtils.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_Constants.H>
#include <ERF_SurfacePrecipitation.H>

/**
 * @file ERF_FireAtmState.H
 * @brief Optional atmospheric inputs of a fire model, computed only when the
 *        model asks for them (FireModelBase::RequestedInputs): the air
 *        temperature and relative humidity of the lowest cell, and the surface
 *        precipitation of every column with the rain-rate helpers a fuel
 *        moisture model uses.
 */

/**
 * @brief Compute relative humidity from conserved atmospheric variables at k=0.
 *
 * Computes RH from pressure, temperature, and water vapor mixing ratio at
 * the first atmospheric level (k=0). Uses the ERF microphysics saturation
 * vapor pressure functions.
 *
 * @param[out] RH_out        Output relative humidity MultiFab [fraction 0-1]
 *                           Must be allocated with ncomp=1 on same grid as inputs
 * @param[in]  S_atm         Conservative state at k=0
 *                           Contains RhoTheta (comp 1) and RhoQ1 (comp 4)
 * @param[in]  geom          Geometry at level 0
 *
 * @note
 * - Output is clamped to [0, 1] range
 * - Assumes S_atm contains standard conserved variables (Rho, RhoTheta, RhoQ1)
 * - Temperature extracted from RhoTheta using standard EOS
 * - Pressure computed from RhoTheta and qv using getPgivenRTh
 */
inline
void compute_rh_from_conservative(amrex::MultiFab& RH_out,
                                  const amrex::MultiFab& S_atm,
                                  const amrex::Geometry& /*geom*/)
{
    using namespace amrex;

    // Precompute rdOcp = R_d / Cp_d for use in Exner function
    // Cp_d = 1004.5 J/(kg-K), R_d = 287.0 J/(kg-K)
    constexpr Real rdOcp = R_d / Cp_d;

    // Iterate over all boxes
    for (MFIter mfi(RH_out); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real>       rh_arr = RH_out.array(mfi);
        Array4<const Real> S_arr  = S_atm.const_array(mfi);

        // Kernel to compute RH at k=0
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;  // First atmospheric level

            // Extract conserved variables at k=0
            Real rho     = S_arr(i, j, k, Rho_comp);           // Density
            Real rho_th  = S_arr(i, j, k, RhoTheta_comp);     // RhoTheta

            // Compute mixing ratios (guard against moisture not being present)
            Real qv = 0.0_rt;  // Water vapor mixing ratio [kg/kg]
            if (RhoQ1_comp < S_arr.nComp()) {
                Real rho_qv  = S_arr(i, j, k, RhoQ1_comp);    // RhoQ1 (water vapor)
                qv = rho_qv / rho;
            }

            // Compute potential temperature
            Real theta = rho_th / rho;  // Potential temperature [K]

            // Compute pressure from RhoTheta and qv using standard EOS
            Real p = getPgivenRTh(rho_th, qv);  // Pressure [Pa]

            // Convert pressure to mbar for saturation vapor pressure functions
            Real p_mbar = p * 0.01_rt;  // Pa to mbar (1 mbar = 100 Pa)

            // Compute temperature from potential temperature
            // T = theta * (p / p_0)^(R_d / c_p)
            // Use Exner function: pi = (p / p_0)^(R_d / c_p)
            Real exner = getExnergivenP(p, rdOcp);
            Real T_K = theta * exner;

            // Compute saturation vapor pressure [mbar]
            Real e_sat = erf_esatw(T_K, false);  // Use Flatau polynomial

            // Compute saturation mixing ratio [kg/kg]
            Real q_sat = erf_qsat_from_esat(e_sat, p_mbar);

            // Compute relative humidity as RH = qv / qsat
            Real RH = 0.0_rt;
            if (q_sat > 1.0e-10_rt) {
                RH = qv / q_sat;
            }

            // Clamp RH to valid range [0, 1]
            RH = amrex::min(1.0_rt, amrex::max(0.0_rt, RH));

            rh_arr(i, j, k) = RH;
        });
    }
}

/**
 * @brief Air temperature [K] at k = 0 from the conserved state.
 *
 * T = theta * (p / p_0)^(R_d / c_p), with p from RhoTheta and qv as in
 * compute_rh_from_conservative. Fuel-moisture curves and ignition criteria take
 * an air temperature, not the potential temperature of the first cell, which at
 * 840 hPa is 14-16 K above it.
 *
 * @param[out] T_out  2D MultiFab (one cell in z) on the level's boxes
 * @param[in]  S_atm  Conserved state (Rho, RhoTheta, RhoQ1 when moist)
 */
inline
void compute_t_from_conservative(amrex::MultiFab& T_out,
                                 const amrex::MultiFab& S_atm)
{
    using namespace amrex;
    constexpr Real rdOcp = R_d / Cp_d;
    for (MFIter mfi(T_out); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real>       t_arr = T_out.array(mfi);
        Array4<const Real> S_arr = S_atm.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
            const int i = iv[0];
            const int j = iv[1];
            const int k = 0;
            const Real rho    = S_arr(i, j, k, Rho_comp);
            const Real rho_th = S_arr(i, j, k, RhoTheta_comp);
            Real qv = 0.0_rt;
            if (RhoQ1_comp < S_arr.nComp()) {
                qv = S_arr(i, j, k, RhoQ1_comp) / rho;
            }
            const Real p     = getPgivenRTh(rho_th, qv);
            const Real exner = getExnergivenP(p, rdOcp);
            t_arr(i, j, k)   = (rho_th / rho) * exner;
        });
    }
}

// Rain per column for a fuel moisture model, from the microphysics surface
// precipitation accumulators, for a fire model that requests
// FireInput_SurfacePrecip.
//
// Every Eulerian microphysics scheme with precipitation keeps a cumulative
// surface accumulator per atmospheric column (Kessler's rain_accum, SAM's and
// Morrison's rain, snow and graupel accumulations, ...), exposed through
// Microphysics::Get_Surface_Precip_Accumulation_Ptrs() with a factor that
// turns the scheme's native unit into kg/m2 of liquid water. One kg/m2 of
// liquid water is one millimetre of depth, so the change of the accumulation
// over a step, divided by the step, is the rain rate in mm/hr that the
// time-lag and stick moisture updates take. The fire model keeps the
// accumulation it saw at the end of its previous step (checkpointed) and
// forms that difference every step; a fresh start begins at zero, which is
// where every scheme starts its accumulators.
//
// Snow, graupel and hail count with their liquid-water equivalent: the fuel
// models have no frozen-precipitation branch. There is no canopy interception.

/**
 * @brief Rain rate [mm/hr] from the change of the surface accumulation over one step.
 *
 * @param[in] accum_now_kgm2   accumulation at the end of the step [kg/m2 == mm]
 * @param[in] accum_prev_kgm2  accumulation at the end of the previous step [kg/m2]
 * @param[in] dt_s             the step [s]
 * @return (now - prev) * 3600 / dt; zero for a non-positive step or a non-positive
 *         change (an accumulator that was reset, or round-off), and zero for NaN.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real fire_precip_rate_mm_hr (amrex::Real accum_now_kgm2,
                                    amrex::Real accum_prev_kgm2,
                                    amrex::Real dt_s) noexcept
{
    if (!(dt_s > amrex::Real(0.0))) { return amrex::Real(0.0); }
    const amrex::Real d = accum_now_kgm2 - accum_prev_kgm2;
    if (!(d > amrex::Real(0.0))) { return amrex::Real(0.0); }
    return d * amrex::Real(3600.0) / dt_s;
}

/**
 * @brief Total liquid-water-equivalent surface precipitation accumulation [kg/m2]
 *        of every column, from the scheme-native accumulators.
 *
 * Fills accum_2d, a one-component MultiFab on the level's boxes coarsened to one
 * cell in z (box n of accum_2d is box n of the accumulators, with the same
 * owner), from the k = klo plane of the sources. Schemes that provide a total
 * give it directly; schemes that provide the species give their sum, so
 * nothing is counted twice (the same rule as the 2D plotfile precipitation
 * diagnostics). Without any source the field is zero.
 */
inline void fire_surface_precip_accum_k0 (amrex::MultiFab& accum_2d,
                                          const SurfacePrecipAccumulationSources& sources,
                                          int klo)
{
    using namespace amrex;
    accum_2d.setVal(Real(0.0));
    if (!surface_precip_has_any_source(sources)) { return; }

    const bool has_total   = surface_precip_has_source(sources.total);
    const bool has_rain    = surface_precip_has_source(sources.rain);
    const bool has_snow    = surface_precip_has_source(sources.snow);
    const bool has_graupel = surface_precip_has_source(sources.graupel);
    const bool has_hail    = surface_precip_has_source(sources.hail);
    const Real f_total   = sources.total.native_to_kg_m2;
    const Real f_rain    = sources.rain.native_to_kg_m2;
    const Real f_snow    = sources.snow.native_to_kg_m2;
    const Real f_graupel = sources.graupel.native_to_kg_m2;
    const Real f_hail    = sources.hail.native_to_kg_m2;

    for (MFIter mfi(accum_2d, false); mfi.isValid(); ++mfi) {
        Array4<Real> a = accum_2d.array(mfi);
        Array4<const Real> total   = has_total   ? sources.total.accum->const_array(mfi)   : Array4<const Real>{};
        Array4<const Real> rain    = has_rain    ? sources.rain.accum->const_array(mfi)    : Array4<const Real>{};
        Array4<const Real> snow    = has_snow    ? sources.snow.accum->const_array(mfi)    : Array4<const Real>{};
        Array4<const Real> graupel = has_graupel ? sources.graupel.accum->const_array(mfi) : Array4<const Real>{};
        Array4<const Real> hail    = has_hail    ? sources.hail.accum->const_array(mfi)    : Array4<const Real>{};
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept
        {
            const int i = iv[0], j = iv[1];
            Real v;
            if (has_total) {
                v = total(i, j, klo) * f_total;
            } else {
                v = Real(0.0);
                if (has_rain)    { v += rain(i, j, klo)    * f_rain; }
                if (has_snow)    { v += snow(i, j, klo)    * f_snow; }
                if (has_graupel) { v += graupel(i, j, klo) * f_graupel; }
                if (has_hail)    { v += hail(i, j, klo)    * f_hail; }
            }
            a(i, j, 0) = v;
        });
    }
}

/**
 * @brief Rain rate [mm/hr] of every column over the step just taken, and roll
 *        the snapshot forward.
 *
 * rate_2d = fire_precip_rate_mm_hr(now, prev, dt) per column when prev_valid,
 * zero otherwise (a restart from a checkpoint without the snapshot); then
 * prev_2d = now_2d. All three MultiFabs share boxes and distribution.
 */
inline void fire_precip_rate_from_accum (amrex::MultiFab& rate_2d,
                                         amrex::MultiFab& prev_2d,
                                         const amrex::MultiFab& now_2d,
                                         amrex::Real dt_s,
                                         bool prev_valid)
{
    using namespace amrex;
    for (MFIter mfi(rate_2d, false); mfi.isValid(); ++mfi) {
        Array4<Real> r = rate_2d.array(mfi);
        Array4<Real> p = prev_2d.array(mfi);
        Array4<const Real> n = now_2d.const_array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept
        {
            const int i = iv[0], j = iv[1];
            r(i, j, 0) = prev_valid ? fire_precip_rate_mm_hr(n(i, j, 0), p(i, j, 0), dt_s) : Real(0.0);
            p(i, j, 0) = n(i, j, 0);
        });
    }
}

#endif  // ERF_FIRE_ATM_STATE_H
