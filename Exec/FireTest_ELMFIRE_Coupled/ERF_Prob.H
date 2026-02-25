#ifndef ERF_PROB_H
#define ERF_PROB_H

#include <string>

#include "AMReX_REAL.H"
#include "AMReX_GpuMemory.H"

#include "ERF_ProbCommon.H"

/**
 * ProbParm: Problem-specific parameters for ELMFIRE atmospheric coupling test
 *
 * This test case provides a uniform atmosphere with initial wind.
 * Fire will extract wind from ERF atmosphere and return heat flux.
 */
struct ProbParm : ProbParmDefaults {
    // Initial atmospheric state
    amrex::Real rho_0 = 1.2;    // kg/m³
    amrex::Real T_0 = 300.0;    // K
    amrex::Real U_0 = 5.0;      // m/s (initial eastward wind)
    amrex::Real V_0 = 0.0;      // m/s
    amrex::Real W_0 = 0.0;      // m/s
};

/**
 * Problem: ELMFIRE atmospheric coupling test problem
 *
 * Provides a uniform atmosphere with initial wind. Fire wind comes
 * from ERF atmosphere (atmospheric coupling mode).
 */
class Problem : public ProblemBase
{
public:
    Problem();  // Constructor - reads parameters from input file

    // Use constant density HSE initialization
    #include "Prob/ERF_InitConstantDensityHSE.H"

    void init_custom_pert(
        const amrex::Box&  bx,
        const amrex::Box& xbx,
        const amrex::Box& ybx,
        const amrex::Box& zbx,
        amrex::Array4<amrex::Real const> const& state,
        amrex::Array4<amrex::Real      > const& state_pert,
        amrex::Array4<amrex::Real      > const& x_vel_pert,
        amrex::Array4<amrex::Real      > const& y_vel_pert,
        amrex::Array4<amrex::Real      > const& z_vel_pert,
        amrex::Array4<amrex::Real      > const& r_hse,
        amrex::Array4<amrex::Real      > const& p_hse,
        amrex::Array4<amrex::Real const> const& z_nd,
        amrex::Array4<amrex::Real const> const& z_cc,
        amrex::GeometryData const& geomdata,
        amrex::Array4<amrex::Real const> const& mf_m,
        amrex::Array4<amrex::Real const> const& mf_u,
        amrex::Array4<amrex::Real const> const& mf_v,
        const SolverChoice& sc,
        const int lev) override;

protected:
    std::string name() {
        return "FireTest_ELMFIRE_Coupled";
    }

    // Flat terrain (no custom terrain needed)
    void init_custom_terrain(
        const amrex::Geometry& geom,
        amrex::FArrayBox& z_phys_nd,
        const amrex::Real& time
    ) override;

private:
    ProbParm parms;
};

#endif
