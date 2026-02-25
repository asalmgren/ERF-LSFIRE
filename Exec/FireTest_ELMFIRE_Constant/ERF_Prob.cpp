#include "ERF_Prob.H"
#include "ERF_EOS.H"
#include "AMReX_ParmParse.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem()
{
    ParmParse pp("prob");

    // Read atmospheric initial conditions
    pp.query("rho_0", parms.rho_0);
    pp.query("T_0", parms.T_0);
    pp.query("U_0", parms.U_0);
    pp.query("V_0", parms.V_0);
    pp.query("W_0", parms.W_0);

    // Initialize base state for HSE
    init_base_parms(parms.rho_0, parms.T_0);

    // Print configuration
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "\n";
        Print() << "========================================\n";
        Print() << "ELMFIRE Constant Wind Test Configuration\n";
        Print() << "========================================\n";
        Print() << "  Initial atmosphere:\n";
        Print() << "    Density: " << parms.rho_0 << " kg/m³\n";
        Print() << "    Temperature: " << parms.T_0 << " K\n";
        Print() << "    Wind: (" << parms.U_0 << ", " << parms.V_0 << ", " << parms.W_0 << ") m/s\n";
        Print() << "\n";
        Print() << "  Fire wind is specified via fire.constant_wind_u/v\n";
        Print() << "  (see input file for fire parameters)\n";
        Print() << "========================================\n";
        Print() << "\n";
    }
}

void Problem::init_custom_terrain(
    const amrex::Geometry& /*geom*/,
    amrex::FArrayBox& /*z_phys_nd*/,
    const amrex::Real& /*time*/)
{
    // Flat terrain - nothing to do
}
