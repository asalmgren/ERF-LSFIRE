#ifndef ERF_PROB_H
#define ERF_PROB_H

#include <string>

#include "AMReX_REAL.H"
#include "AMReX_GpuMemory.H"

#include "ERF_ProbCommon.H"

/**
 * ProbParm: Problem-specific parameters for ELMFIRE constant wind test
 *
 * This test case provides a uniform quiescent atmosphere for testing
 * ELMFIRE fire spread with constant wind driven by the fire model itself.
 */
struct ProbParm : ProbParmDefaults {
    // Initial atmospheric state (quiescent atmosphere)
    amrex::Real rho_0 = 1.2;    // kg/m³
    amrex::Real T_0 = 300.0;    // K
    amrex::Real U_0 = 0.0;      // m/s (no atmospheric wind)
    amrex::Real V_0 = 0.0;      // m/s
    amrex::Real W_0 = 0.0;      // m/s
};

/**
 * Problem: ELMFIRE constant wind test problem
 *
 * Provides a uniform quiescent atmosphere. Fire wind is specified
 * in the input file via fire.constant_wind_u and fire.constant_wind_v.
 */
class Problem : public ProblemBase
{
public:
    Problem();  // Constructor - reads parameters from input file

    // Use constant density HSE initialization
    #include "Prob/ERF_InitConstantDensityHSE.H"

protected:
    std::string name() {
        return "FireTest_ELMFIRE_Constant";
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
