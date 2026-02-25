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
        Print() << "ELMFIRE Atmospheric Coupling Test Configuration\n";
        Print() << "========================================\n";
        Print() << "  Initial atmosphere:\n";
        Print() << "    Density: " << parms.rho_0 << " kg/m³\n";
        Print() << "    Temperature: " << parms.T_0 << " K\n";
        Print() << "    Wind: (" << parms.U_0 << ", " << parms.V_0 << ", " << parms.W_0 << ") m/s\n";
        Print() << "\n";
        Print() << "  Fire will extract wind from ERF atmosphere\n";
        Print() << "  Fire will return heat flux to ERF atmosphere\n";
        Print() << "  (Two-way coupling active)\n";
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

void Problem::init_custom_pert(
    const amrex::Box&  /*bx*/,
    const amrex::Box& xbx,
    const amrex::Box& ybx,
    const amrex::Box& zbx,
    amrex::Array4<amrex::Real const> const& /*state*/,
    amrex::Array4<amrex::Real      > const& /*state_pert*/,
    amrex::Array4<amrex::Real      > const& x_vel_pert,
    amrex::Array4<amrex::Real      > const& y_vel_pert,
    amrex::Array4<amrex::Real      > const& z_vel_pert,
    amrex::Array4<amrex::Real      > const& /*r_hse*/,
    amrex::Array4<amrex::Real      > const& /*p_hse*/,
    amrex::Array4<amrex::Real const> const& /*z_nd*/,
    amrex::Array4<amrex::Real const> const& /*z_cc*/,
    amrex::GeometryData const& /*geomdata*/,
    amrex::Array4<amrex::Real const> const& /*mf_m*/,
    amrex::Array4<amrex::Real const> const& /*mf_u*/,
    amrex::Array4<amrex::Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    // Initialize velocities to uniform values
    Real U_0 = parms.U_0;
    Real V_0 = parms.V_0;
    Real W_0 = parms.W_0;

    // Set x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        x_vel_pert(i,j,k) = U_0;
    });

    // Set y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        y_vel_pert(i,j,k) = V_0;
    });

    // Set z-velocity
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        z_vel_pert(i,j,k) = W_0;
    });
}
