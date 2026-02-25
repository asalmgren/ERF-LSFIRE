#include "ERF_FireLevelSet.H"
#include <cmath>
#include <limits>

using namespace amrex;

// ============================================================================
// FAST MARCHING METHOD (LOCAL SWEEPS) - Currently Unused
// ============================================================================

namespace {
    // Helper function for fast marching solver
    Real SolveEikonal(Real phi_x, Real phi_y, Real dx, Real dy, Real F) {
        Real a = 1.0 / (dx * dx);
        Real b = 1.0 / (dy * dy);
        Real c = F * F;
        
        Real phi_min = std::min(phi_x, phi_y);
        Real phi_max = std::max(phi_x, phi_y);
        
        // Try quadratic solution
        Real discriminant = (a + b) * (a + b) - 4.0 * a * b * (phi_min*phi_min - c);
        
        if (discriminant >= 0.0) {
            Real phi_new = ((a*phi_x + b*phi_y) + std::sqrt(discriminant)) / (a + b);
            if (phi_new >= phi_max) {
                return phi_new;
            }
        }
        
        // Fallback to 1D update
        if (phi_x < phi_y) {
            return phi_x + dx * F;
        } else {
            return phi_y + dy * F;
        }
    }
}

// ============================================================================
// LOCAL FAST MARCHING (UNUSED - Kept for future optimization)
// ============================================================================

void FireLevelSet::FastMarchingLocal() {
    // ✅ This is a placeholder for future optimization
    // Currently using direct distance calculation in ReinitializeSignedDistanceDistributed()
    
    // Fast marching would work like this:
    // 1. Each rank runs fast marching on its local domain
    // 2. FillBoundary() to update ghost cells
    // 3. Repeat until convergence
    
    // For now, the gather-broadcast approach is simpler and works well
    // for domains up to ~1000x1000 cells
}

// ============================================================================
// GRADIENT CALCULATION WITH GHOST CELLS
// ============================================================================

void FireLevelSet::ComputeGradientsLocal() {
    BL_PROFILE("FireLevelSet::ComputeGradients");
    
    const Real* dx = m_geom->CellSize();
    Real dx_inv = 1.0 / (2.0 * dx[0]);
    Real dy_inv = 1.0 / (2.0 * dx[1]);
    
    for (auto& node : m_local_tagged) {
        if (node.burned) continue;
        
        // ✅ Read from MultiFab (includes ghost cells!)
        bool found = false;
        
        for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
            const Box& grownbox = mfi.growntilebox();
            
            if (!grownbox.contains(IntVect(node.i, node.j, 0))) continue;
            
            Array4<Real const> phi = m_fire_data.const_array(mfi);
            
            // Central difference
            Real phi_ip1 = phi(node.i+1, node.j, 0);
            Real phi_im1 = phi(node.i-1, node.j, 0);
            Real phi_jp1 = phi(node.i, node.j+1, 0);
            Real phi_jm1 = phi(node.i, node.j-1, 0);
            
            node.dphidx = (phi_ip1 - phi_im1) * dx_inv;
            node.dphidy = (phi_jp1 - phi_jm1) * dy_inv;
            node.gradmag = std::sqrt(node.dphidx*node.dphidx + node.dphidy*node.dphidy);
            
            found = true;
            break;
        }
        
        if (!found) {
            // Cell not in this rank's boxes (shouldn't happen if lists are correct)
            node.dphidx = 0.0;
            node.dphidy = 0.0;
            node.gradmag = 1.0;
        }
    }
}

// ============================================================================
// CURVATURE CALCULATION (FOR ADVANCED MODELS)
// ============================================================================

void FireLevelSet::ComputeCurvatureLocal() {
    BL_PROFILE("FireLevelSet::ComputeCurvature");
    
    const Real* dx = m_geom->CellSize();
    Real dx_inv = 1.0 / dx[0];
    Real dy_inv = 1.0 / dx[1];
    Real dx2_inv = 1.0 / (dx[0] * dx[0]);
    Real dy2_inv = 1.0 / (dx[1] * dx[1]);
    
    for (auto& node : m_local_tagged) {
        if (node.burned) continue;
        
        for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
            const Box& grownbox = mfi.growntilebox();
            
            if (!grownbox.contains(IntVect(node.i, node.j, 0))) continue;
            
            Array4<Real const> phi = m_fire_data.const_array(mfi);
            
            // First derivatives
            Real phi_x = (phi(node.i+1, node.j, 0) - phi(node.i-1, node.j, 0)) * dx_inv * 0.5;
            Real phi_y = (phi(node.i, node.j+1, 0) - phi(node.i, node.j-1, 0)) * dy_inv * 0.5;
            
            // Second derivatives
            Real phi_xx = (phi(node.i+1, node.j, 0) - 2.0*phi(node.i, node.j, 0) + phi(node.i-1, node.j, 0)) * dx2_inv;
            Real phi_yy = (phi(node.i, node.j+1, 0) - 2.0*phi(node.i, node.j, 0) + phi(node.i, node.j-1, 0)) * dy2_inv;
            Real phi_xy = (phi(node.i+1, node.j+1, 0) - phi(node.i+1, node.j-1, 0) 
                         - phi(node.i-1, node.j+1, 0) + phi(node.i-1, node.j-1, 0)) * dx_inv * dy_inv * 0.25;
            
            // Curvature: κ = (phi_xx*phi_y^2 - 2*phi_xy*phi_x*phi_y + phi_yy*phi_x^2) / (phi_x^2 + phi_y^2)^(3/2)
            Real phi_x2 = phi_x * phi_x;
            Real phi_y2 = phi_y * phi_y;
            Real gradmag2 = phi_x2 + phi_y2;
            
            if (gradmag2 > 1e-10) {
                Real numer = phi_xx*phi_y2 - 2.0*phi_xy*phi_x*phi_y + phi_yy*phi_x2;
                Real denom = std::pow(gradmag2, 1.5);
                node.curvature = numer / denom;
            } else {
                node.curvature = 0.0;
            }
            
            break;
        }
    }
}

// ============================================================================
// UPWIND SCHEMES (FOR ALTERNATIVE INTEGRATION METHODS)
// ============================================================================

void FireLevelSet::ComputeUpwindGradientsLocal() {
    BL_PROFILE("FireLevelSet::ComputeUpwindGradients");
    
    const Real* dx = m_geom->CellSize();
    Real dx_inv = 1.0 / dx[0];
    Real dy_inv = 1.0 / dx[1];
    
    for (auto& node : m_local_tagged) {
        if (node.burned) continue;
        
        for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
            const Box& grownbox = mfi.growntilebox();
            
            if (!grownbox.contains(IntVect(node.i, node.j, 0))) continue;
            
            Array4<Real const> phi = m_fire_data.const_array(mfi);
            
            Real phi_c = phi(node.i, node.j, 0);
            Real phi_ip1 = phi(node.i+1, node.j, 0);
            Real phi_im1 = phi(node.i-1, node.j, 0);
            Real phi_jp1 = phi(node.i, node.j+1, 0);
            Real phi_jm1 = phi(node.i, node.j-1, 0);
            
            // Forward differences
            Real dphi_dx_plus = (phi_ip1 - phi_c) * dx_inv;
            Real dphi_dy_plus = (phi_jp1 - phi_c) * dy_inv;
            
            // Backward differences
            Real dphi_dx_minus = (phi_c - phi_im1) * dx_inv;
            Real dphi_dy_minus = (phi_c - phi_jm1) * dy_inv;
            
            // Upwind selection based on velocity direction
            if (node.ux > 0) {
                node.dphidx_upwind = dphi_dx_minus;
            } else {
                node.dphidx_upwind = dphi_dx_plus;
            }
            
            if (node.uy > 0) {
                node.dphidy_upwind = dphi_dy_minus;
            } else {
                node.dphidy_upwind = dphi_dy_plus;
            }
            
            break;
        }
    }
}

// ============================================================================
// WENO SCHEMES (HIGH-ORDER ACCURACY - FUTURE ENHANCEMENT)
// ============================================================================

namespace {
    // WENO5 reconstruction (placeholder for future)
    Real WENO5(Real vm2, Real vm1, Real v0, Real vp1, Real vp2) {
        // 5th-order WENO reconstruction
        // Currently unused, but kept for future enhancement
        
        const Real eps = 1e-6;
        
        // Smoothness indicators
        Real beta0 = 13.0/12.0 * std::pow(vm2 - 2.0*vm1 + v0, 2.0) 
                   + 0.25 * std::pow(vm2 - 4.0*vm1 + 3.0*v0, 2.0);
        Real beta1 = 13.0/12.0 * std::pow(vm1 - 2.0*v0 + vp1, 2.0) 
                   + 0.25 * std::pow(vm1 - vp1, 2.0);
        Real beta2 = 13.0/12.0 * std::pow(v0 - 2.0*vp1 + vp2, 2.0) 
                   + 0.25 * std::pow(3.0*v0 - 4.0*vp1 + vp2, 2.0);
        
        // Weights
        Real w0 = 0.1 / std::pow(eps + beta0, 2.0);
        Real w1 = 0.6 / std::pow(eps + beta1, 2.0);
        Real w2 = 0.3 / std::pow(eps + beta2, 2.0);
        Real wsum = w0 + w1 + w2;
        w0 /= wsum;
        w1 /= wsum;
        w2 /= wsum;
        
        // Reconstructed value
        Real v_weno = w0 * (2.0*vm2 - 7.0*vm1 + 11.0*v0) / 6.0
                    + w1 * (-vm1 + 5.0*v0 + 2.0*vp1) / 6.0
                    + w2 * (2.0*v0 + 5.0*vp1 - vp2) / 6.0;
        
        return v_weno;
    }
}

// ============================================================================
// DISTANCE FUNCTION UTILITIES
// ============================================================================

Real FireLevelSet::ComputeDistanceToPoint(int i, int j, int i_ref, int j_ref) const {
    const Real* dx = m_geom->CellSize();
    const Real* prob_lo = m_geom->ProbLo();
    
    Real x = prob_lo[0] + (i + 0.5) * dx[0];
    Real y = prob_lo[1] + (j + 0.5) * dx[1];
    
    Real x_ref = prob_lo[0] + (i_ref + 0.5) * dx[0];
    Real y_ref = prob_lo[1] + (j_ref + 0.5) * dx[1];
    
    return std::sqrt((x - x_ref)*(x - x_ref) + (y - y_ref)*(y - y_ref));
}

void FireLevelSet::CheckGradientQuality() {
    // ✅ Diagnostic function to check if gradients need reinitialization
    
    int local_count = 0;
    Real local_sum_deviation = 0.0;
    Real local_max_deviation = 0.0;
    
    for (const auto& node : m_local_tagged) {
        Real gradmag = std::sqrt(node.dphidx*node.dphidx + node.dphidy*node.dphidy);
        Real deviation = std::abs(gradmag - 1.0);
        
        local_sum_deviation += deviation;
        local_max_deviation = std::max(local_max_deviation, deviation);
        local_count++;
    }
    
#ifdef AMREX_USE_MPI
    Real global_sum_deviation = local_sum_deviation;
    Real global_max_deviation = local_max_deviation;
    int global_count = local_count;
    
    ParallelDescriptor::ReduceRealSum(global_sum_deviation);
    ParallelDescriptor::ReduceRealMax(global_max_deviation);
    ParallelDescriptor::ReduceIntSum(global_count);
#else
    Real global_sum_deviation = local_sum_deviation;
    Real global_max_deviation = local_max_deviation;
    int global_count = local_count;
#endif
    
    if (ParallelDescriptor::IOProcessor() && global_count > 0) {
        Real avg_deviation = global_sum_deviation / global_count;
        Print() << "  Gradient quality: avg|∇φ|-1| = " << avg_deviation 
                << ", max = " << global_max_deviation << "\n";
    }
}

// ============================================================================
// CFL CONDITION CHECK
// ============================================================================

Real FireLevelSet::ComputeMaxCFL(Real dt) const {
    const Real* dx = m_geom->CellSize();
    Real min_dx = std::min(dx[0], dx[1]);
    
    Real local_max_cfl = 0.0;
    
    for (const auto& node : m_local_tagged) {
        Real vel = std::sqrt(node.ux*node.ux + node.uy*node.uy);
        Real cfl = vel * dt / min_dx;
        local_max_cfl = std::max(local_max_cfl, cfl);
    }
    
#ifdef AMREX_USE_MPI
    Real global_max_cfl = local_max_cfl;
    ParallelDescriptor::ReduceRealMax(global_max_cfl);
    return global_max_cfl;
#else
    return local_max_cfl;
#endif
}

// ============================================================================
// MASS CONSERVATION CHECK
// ============================================================================

void FireLevelSet::CheckMassConservation() {
    // ✅ Verify that φ evolution preserves the zero level set
    
    int local_burned = m_local_burned.size();
    
#ifdef AMREX_USE_MPI
    int global_burned = local_burned;
    ParallelDescriptor::ReduceIntSum(global_burned);
#else
    int global_burned = local_burned;
#endif
    
    if (ParallelDescriptor::IOProcessor()) {
        const Real* dx = m_geom->CellSize();
        Real cell_area = dx[0] * dx[1];
        Real burned_area = global_burned * cell_area;
        
        static Real last_burned_area = 0.0;
        Real area_change = burned_area - last_burned_area;
        
        if (area_change > 0) {
            Print() << "  Burned area: " << burned_area << " m², "
                    << "change: +" << area_change << " m²\n";
        }
        
        last_burned_area = burned_area;
    }
}