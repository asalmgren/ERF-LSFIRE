#include "ERF_FireLevelSet.H"
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <cmath>

using namespace amrex;

void FireLevelSet::Define(SolverChoice& sc)
{
    m_met_source = sc.fire_met_source;
    
    amrex::ParmParse pp("prob");
    
    pp.query("fuel_density", m_fuel_density);
    pp.query("fuel_moisture", m_fuel_moisture);
    pp.query("fuel_heat", m_fuel_heat);
    pp.query("fuel_depth", m_fuel_depth);
    pp.query("fuel_savr", m_fuel_savr);
    
    pp.query("fire_ros_base", m_ros_base);
    pp.query("fire_wind_factor", m_wind_factor);
    pp.query("fire_slope_factor", m_slope_factor);
    pp.query("fire_refine_distance", m_refine_distance);
    
    pp.query("fire_ignition_time", m_ignition_time);
    pp.query("fire_ignition_x", m_ignition_x);
    pp.query("fire_ignition_y", m_ignition_y);
    pp.query("fire_ignition_radius", m_ignition_radius);
    
    if (m_met_source == "Constant") {
        pp.query("fire_constant_u", m_constant_u);
        pp.query("fire_constant_v", m_constant_v);
    }
    
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "\n=== FireLevelSet Configuration ===\n";
        Print() << "  Met source: " << m_met_source << "\n";
        if (m_met_source == "Constant") {
            Print() << "    Constant U = " << m_constant_u << " m/s\n";
            Print() << "    Constant V = " << m_constant_v << " m/s\n";
        }
        Print() << "  Base ROS: " << m_ros_base << " m/s\n";
        Print() << "  Wind factor: " << m_wind_factor << "\n";
        Print() << "  Slope factor: " << m_slope_factor << "\n";
        Print() << "  Fuel density: " << m_fuel_density << " kg/m³\n";
        Print() << "  Fuel moisture: " << m_fuel_moisture*100 << " %\n";
        Print() << "  Ignition at: (" << m_ignition_x << ", " << m_ignition_y << ")\n";
        Print() << "  Ignition time: " << m_ignition_time << " s\n";
        Print() << "  Ignition radius: " << m_ignition_radius << " m\n";
        Print() << "===================================\n\n";
    }
}

void FireLevelSet::Init(int lev,
                       const MultiFab& cons_in,
                       const Geometry& geom,
                       const Real& dt_advance)
{
    m_lev = lev;
    m_dt = dt_advance;
    m_geom = &geom;
    
    // Create 2D box array from 3D
    const BoxArray& ba3d = cons_in.boxArray();
    BoxArray ba2d(ba3d.size());
    for (int i = 0; i < ba3d.size(); ++i) {
        Box bx = ba3d[i];
        bx.setSmall(2, 0);
        bx.setBig(2, 0);
        ba2d.set(i, bx);
    }
    
    const DistributionMapping& dm = cons_in.DistributionMap();
    IntVect nghost(1, 1, 0);
    
    m_fire_data.define(ba2d, dm, FireVar::NumVars, nghost);
    m_fire_data.setVal(0.0);
    
    // Initialize slope storage
    m_slope_x = std::make_unique<MultiFab>(ba2d, dm, 1, nghost);
    m_slope_y = std::make_unique<MultiFab>(ba2d, dm, 1, nghost);
    m_slope_x->setVal(0.0);
    m_slope_y->setVal(0.0);
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> fire_arr = m_fire_data.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            fire_arr(i,j,0,FireVar::fuel_frac) = 1.0;
            fire_arr(i,j,0,FireVar::fuel_time) = 0.0;
            fire_arr(i,j,0,FireVar::heat_flux) = 0.0;
            fire_arr(i,j,0,FireVar::phi) = 1.0e10;
        });
    }
    
    m_fire_data.FillBoundary(m_geom->periodicity());
    m_active = false;
}

void FireLevelSet::Advance(const Real& dt,
                          const MultiFab* U,
                          const MultiFab* V,
                          const MultiFab* terrain)
{
    m_wind_u = U;
    m_wind_v = V;
    m_terrain = terrain;

    static Real accumulated_time = 0.0;
    accumulated_time += dt;
 
    BL_PROFILE("FireLevelSet::Advance");

    if (ParallelDescriptor::IOProcessor()) {
        Print() << "\n=== FireLevelSet::Advance ===\n";
        Print() << "  Time: " << accumulated_time << " s\n";
        Print() << "  dt: " << dt << " s\n";
        Print() << "  Active: " << (m_active ? "YES" : "NO") << "\n";
    }
    
    // Check for ignition
    if (!m_active && accumulated_time >= m_ignition_time) {
        Vector<Real> ignition_loc(2);
        ignition_loc[0] = m_ignition_x;
        ignition_loc[1] = m_ignition_y;
        
        SetIgnition(ignition_loc, m_ignition_radius);
        m_active = true;
        
        if (ParallelDescriptor::IOProcessor()) {
            Print() << "🔥 Fire ignited at (" << m_ignition_x << ", " << m_ignition_y << ")\n";
        }
    }
    
    if (!m_active) {
        return;
    }
    
    // Compute terrain slopes once if available
    if (m_terrain) {
        const Real* dx = m_geom->CellSize();
        ComputeTerrainSlopes(m_terrain, dx);
    }
    
    m_fire_data.FillBoundary(m_geom->periodicity());
    
    // Main fire advancement
    AdvectLevelSet(dt);
    
    // Reinitialize every 5 steps
    // static int reinit_counter = 0;
    // reinit_counter++;
    // if (reinit_counter % 5 == 0) {
    //     ReinitializeLevelSet();
    // }
    
    BurnFuel(dt);
    ComputeHeatRelease();
    
    m_fire_data.FillBoundary(m_geom->periodicity());
    
    if (ParallelDescriptor::IOProcessor()) {
        Real fire_area = get_fire_area();
        Print() << "  Fire area: " << fire_area << " m²\n";
        Print() << "===============================\n\n";
    }
}

void FireLevelSet::ComputeTerrainSlopes(const MultiFab* terrain, const Real* dx)
{
    BL_PROFILE("FireLevelSet::ComputeTerrainSlopes");
    
    if (!terrain) return;
    
    for (MFIter mfi(*m_slope_x); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        
        Array4<Real const> terrain_arr = terrain->const_array(mfi);
        Array4<Real> slope_x_arr = m_slope_x->array(mfi);
        Array4<Real> slope_y_arr = m_slope_y->array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            slope_x_arr(i,j,0) = (terrain_arr(i+1,j,0) - terrain_arr(i-1,j,0)) / (2.0 * dx[0]);
            slope_y_arr(i,j,0) = (terrain_arr(i,j+1,0) - terrain_arr(i,j-1,0)) / (2.0 * dx[1]);
        });
    }
    
    m_slope_x->FillBoundary(m_geom->periodicity());
    m_slope_y->FillBoundary(m_geom->periodicity());
}

void FireLevelSet::ComputeSpreadRate(MultiFab& ros_out,
                                     MultiFab& ros_x_out,
                                     MultiFab& ros_y_out)
{
    BL_PROFILE("FireLevelSet::ComputeSpreadRate");
    
    const Real ros_base = m_ros_base;
    const Real wind_factor = m_wind_factor;
    const Real slope_factor = m_slope_factor;
    const bool use_constant_wind = (m_met_source == "Constant");
    const Real const_u = m_constant_u;
    const Real const_v = m_constant_v;
    const bool have_slopes = (m_slope_x && m_slope_y && m_terrain);
    
    const Box& domain = m_geom->Domain();
    const Real* dx = m_geom->CellSize();
    
    // ========== FIX: Much narrower fire front band ==========
    // Old: 3.0 × dx = 900m (way too wide, causes discontinuous fire)
    // New: 1.5 × dx = 450m (tight band around phi=0)
    const Real fire_front_band = 1.5 * amrex::max(dx[0], dx[1]);
    // ========== END FIX ==========
    
    ros_out.setVal(0.0);
    ros_x_out.setVal(0.0);
    ros_y_out.setVal(0.0);
    
    for (MFIter mfi(ros_out); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        
        Array4<Real const> phi_arr = m_fire_data.const_array(mfi);
        Array4<Real const> fuel_arr = m_fire_data.const_array(mfi, FireVar::fuel_frac);
        Array4<Real> ros_arr = ros_out.array(mfi);
        Array4<Real> ros_x_arr = ros_x_out.array(mfi);
        Array4<Real> ros_y_arr = ros_y_out.array(mfi);
        
        Array4<Real const> u_arr, v_arr;
        bool have_wind = false;
        
        if (!use_constant_wind && m_wind_u && m_wind_v) {
            u_arr = m_wind_u->const_array(mfi);
            v_arr = m_wind_v->const_array(mfi);
            have_wind = true;
        }
        
        Array4<Real const> slope_x_arr, slope_y_arr;
        if (have_slopes) {
            slope_x_arr = m_slope_x->const_array(mfi);
            slope_y_arr = m_slope_y->const_array(mfi);
        }
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            bool at_boundary = (i == domain.smallEnd(0) || i == domain.bigEnd(0) ||
                               j == domain.smallEnd(1) || j == domain.bigEnd(1));
            if (at_boundary) return;
            
            Real phi_c = phi_arr(i,j,0,FireVar::phi);
            Real fuel_c = fuel_arr(i,j,0);
            
            if (fuel_c < 0.01) return;
            
            // Only compute ROS in narrow band around phi=0
            bool at_fire_front = (phi_c >= -fire_front_band && phi_c <= fire_front_band);
            if (!at_fire_front) return;
            
            // ========== ADD: Extra check for gradient direction ==========
            // This prevents isolated "fire islands" from forming
            Real phi_e = phi_arr(i+1,j,0,FireVar::phi);
            Real phi_w = phi_arr(i-1,j,0,FireVar::phi);
            Real phi_n = phi_arr(i,j+1,0,FireVar::phi);
            Real phi_s = phi_arr(i,j-1,0,FireVar::phi);
            
            Real phi_x = (phi_e - phi_w) / (2.0*dx[0]);
            Real phi_y = (phi_n - phi_s) / (2.0*dx[1]);
            Real grad_mag = std::sqrt(phi_x*phi_x + phi_y*phi_y);
            
            if (grad_mag < 1.0e-10) return;
            
            // ========== ADDITIONAL CHECK: Must have unburned neighbor ==========
            // This ensures fire only spreads at the actual front, not interior
            bool has_unburned_neighbor = (phi_e > 0.0) || (phi_w > 0.0) || 
                                         (phi_n > 0.0) || (phi_s > 0.0);
            if (!has_unburned_neighbor) return;  // Interior point, don't spread
            // ========== END CHECK ==========
            
            Real nx = phi_x / grad_mag;
            Real ny = phi_y / grad_mag;
            
            Real u_wind, v_wind;
            if (use_constant_wind) {
                u_wind = const_u;
                v_wind = const_v;
            } else if (have_wind) {
                int k_wind = 0;
                u_wind = 0.5 * (u_arr(i,j,k_wind) + u_arr(i+1,j,k_wind));
                v_wind = 0.5 * (v_arr(i,j,k_wind) + v_arr(i,j+1,k_wind));
            } else {
                u_wind = 0.0;
                v_wind = 0.0;
            }
            
            Real wind_dot = (u_wind * nx + v_wind * ny);
            wind_dot = amrex::max(-50.0, amrex::min(wind_dot, 50.0));
            Real wind_effect = wind_factor * wind_dot;
            
            Real slope_effect = 0.0;
            if (have_slopes) {
                Real sx = slope_x_arr(i,j,0);
                Real sy = slope_y_arr(i,j,0);
                Real slope_dot = (sx * nx + sy * ny);
                slope_dot = amrex::max(-1.0, amrex::min(slope_dot, 1.0));
                slope_effect = slope_factor * slope_dot;
            }
            
            Real ros_total = ros_base * (1.0 + wind_effect + slope_effect);
            ros_total = amrex::max(0.0, amrex::min(ros_total, 30.0));
            
            if (!std::isfinite(ros_total)) {
                ros_total = 0.0;
            }
            
            ros_arr(i,j,0) = ros_total;
            ros_x_arr(i,j,0) = ros_total * nx;
            ros_y_arr(i,j,0) = ros_total * ny;
        });
    }
    
    if (ParallelDescriptor::IOProcessor()) {
        Real ros_min = ros_out.min(0);
        Real ros_max = ros_out.max(0);
        Print() << "  ROS: min=" << ros_min << ", max=" << ros_max << " m/s\n";
    }
}

void FireLevelSet::AdvectLevelSet(const Real& dt)
{
    BL_PROFILE("FireLevelSet::AdvectLevelSet");
    
    const Real* dx = m_geom->CellSize();
    
    MultiFab ros(m_fire_data.boxArray(), m_fire_data.DistributionMap(), 1, 1);
    MultiFab ros_x(m_fire_data.boxArray(), m_fire_data.DistributionMap(), 1, 1);
    MultiFab ros_y(m_fire_data.boxArray(), m_fire_data.DistributionMap(), 1, 1);
    
    ComputeSpreadRate(ros, ros_x, ros_y);
    
    MultiFab phi_new(m_fire_data.boxArray(), m_fire_data.DistributionMap(), 1, 0);
    MultiFab::Copy(phi_new, m_fire_data, FireVar::phi, 0, 1, 0);
    
    const Box& domain = m_geom->Domain();
    
    for (MFIter mfi(phi_new); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        
        Array4<Real const> phi_old = m_fire_data.const_array(mfi);
        Array4<Real const> ros_arr = ros.const_array(mfi);
        Array4<Real> phi = phi_new.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            bool at_boundary = (i == domain.smallEnd(0) || i == domain.bigEnd(0) ||
                               j == domain.smallEnd(1) || j == domain.bigEnd(1));
            
            if (at_boundary) {
                phi(i,j,k,0) = phi_old(i,j,0,FireVar::phi);
                return;
            }
            
            Real phi_c = phi_old(i,j,0,FireVar::phi);
            Real ros_local = ros_arr(i,j,0);
            
            // Skip if no ROS here
            if (ros_local < 1.0e-10) {
                phi(i,j,k,0) = phi_c;
                return;
            }
            
            // Get neighbors
            Real phi_e = phi_old(i+1,j,0,FireVar::phi);
            Real phi_w = phi_old(i-1,j,0,FireVar::phi);
            Real phi_n = phi_old(i,j+1,0,FireVar::phi);
            Real phi_s = phi_old(i,j-1,0,FireVar::phi);
            
            // ========== CORRECT: Simple central difference ==========
            Real phi_x = (phi_e - phi_w) / (2.0 * dx[0]);
            Real phi_y = (phi_n - phi_s) / (2.0 * dx[1]);
            Real grad_mag = std::sqrt(phi_x*phi_x + phi_y*phi_y);
            
            if (grad_mag < 1.0e-10) {
                phi(i,j,k,0) = phi_c;
                return;
            }
            
            // Standard level-set advection: ∂φ/∂t = -ROS × |∇φ|
            // Negative sign because fire spreads from phi < 0 to phi > 0
            Real phi_new_val = phi_c - dt * ros_local * grad_mag;
            
            phi(i,j,k,0) = phi_new_val;
        });
    }
    
    if (ParallelDescriptor::IOProcessor()) {
        Real phi_min = phi_new.min(0);
        Real phi_max = phi_new.max(0);
        Real phi_min_old = m_fire_data.min(FireVar::phi);
        Real phi_change = phi_min - phi_min_old;
        
        Print() << "  Phi: min=" << phi_min << " (Δ=" << phi_change << "), max=" << phi_max << "\n";
    }
    
    MultiFab::Copy(m_fire_data, phi_new, 0, FireVar::phi, 1, 0);
}

void FireLevelSet::SetIgnition(const Vector<Real>& location, const Real& radius)
{
    const Real* dx = m_geom->CellSize();
    const Real* prob_lo = m_geom->ProbLo();
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> fire_arr = m_fire_data.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real x = prob_lo[0] + (i + 0.5) * dx[0];
            Real y = prob_lo[1] + (j + 0.5) * dx[1];
            
            Real dist = std::sqrt((x - location[0])*(x - location[0]) + 
                                 (y - location[1])*(y - location[1]));
            
            fire_arr(i,j,k,FireVar::phi) = dist - radius;
            
            if (dist < radius) {
                fire_arr(i,j,k,FireVar::fuel_time) = 0.0;
            }
        });
    }
    
    m_fire_data.FillBoundary(m_geom->periodicity());
}

void FireLevelSet::ReinitializeLevelSet()
{
    BL_PROFILE("FireLevelSet::ReinitializeLevelSet");
    
    const Real* dx = m_geom->CellSize();
    const Real dx_min = std::min(dx[0], dx[1]);
    const Box& domain = m_geom->Domain();
    
    if (ParallelDescriptor::IOProcessor()) {
        Real phi_before = m_fire_data.min(FireVar::phi);
        Print() << "  Reinitializing level-set (phi_min_before=" << phi_before << ")\n";
    }
    
    // ========== SIMPLIFIED REINITIALIZATION ==========
    // Use fast sweeping to restore |∇φ| ≈ 1 near fire front
    // Only operate in narrow band around phi=0
    // ========== END ==========
    
    MultiFab phi_old(m_fire_data.boxArray(), m_fire_data.DistributionMap(), 1, 1);
    MultiFab::Copy(phi_old, m_fire_data, FireVar::phi, 0, 1, 1);
    phi_old.FillBoundary(m_geom->periodicity());
    
    // Identify fire front cells (narrow band)
    const Real band_width = 5.0 * dx_min;  // Only reinitialize ±5 cells from front
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        
        Array4<Real const> phi_arr = phi_old.const_array(mfi);
        Array4<Real> phi_new = m_fire_data.array(mfi, FireVar::phi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            bool at_boundary = (i == domain.smallEnd(0) || i == domain.bigEnd(0) ||
                               j == domain.smallEnd(1) || j == domain.bigEnd(1));
            if (at_boundary) {
                return;  // Don't touch boundaries
            }
            
            Real phi_c = phi_arr(i,j,k);
            
            // ========== ONLY REINITIALIZE NEAR FIRE FRONT ==========
            // Far from front: leave phi unchanged
            if (std::abs(phi_c) > band_width) {
                return;  // Too far from front, don't touch
            }
            
            // Get neighbors
            Real phi_e = phi_arr(i+1,j,k);
            Real phi_w = phi_arr(i-1,j,k);
            Real phi_n = phi_arr(i,j+1,k);
            Real phi_s = phi_arr(i,j-1,k);
            
            // Compute gradient magnitude
            Real phi_x = (phi_e - phi_w) / (2.0 * dx[0]);
            Real phi_y = (phi_n - phi_s) / (2.0 * dx[1]);
            Real grad_mag = std::sqrt(phi_x*phi_x + phi_y*phi_y);
            
            // If gradient is too small or too large, adjust phi slightly
            if (grad_mag > 1.0e-10) {
                Real target_grad = 1.0;  // Want |∇φ| = 1
                
                // Simple correction: adjust phi to match target gradient
                // This is a single Gauss-Seidel-like iteration
                Real correction = (grad_mag - target_grad) / grad_mag;
                Real dtau = 0.3 * dx_min;  // Small pseudo-timestep
                
                // Move phi toward restoring |∇φ| = 1, preserving sign
                Real sign_phi = (phi_c >= 0.0) ? 1.0 : -1.0;
                phi_new(i,j,0) = phi_c - dtau * sign_phi * correction;
                
                // ========== CRITICAL: Preserve sign ==========
                // Ensure we don't flip burned->unburned or vice versa
                if ((phi_c < 0.0 && phi_new(i,j,0) > 0.0) ||
                    (phi_c > 0.0 && phi_new(i,j,0) < 0.0)) {
                    phi_new(i,j,0) = phi_c * 0.9;  // Reduce magnitude but keep sign
                }
            }
        });
    }
    
    m_fire_data.FillBoundary(m_geom->periodicity());
    
    if (ParallelDescriptor::IOProcessor()) {
        Real phi_min = m_fire_data.min(FireVar::phi);
        Real phi_max = m_fire_data.max(FireVar::phi);
        Print() << "    ✓ Reinitialized: φ ∈ [" << phi_min << ", " << phi_max << "]\n";
    }
}

void FireLevelSet::BurnFuel(const Real& dt)
{
    BL_PROFILE("FireLevelSet::BurnFuel");
    
    const Real fuel_burn_rate = m_fuel_density * m_fuel_depth / 60.0;
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> fire_arr = m_fire_data.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real phi = fire_arr(i,j,k,FireVar::phi);
            
            if (phi < 0.0 && phi > -10.0) {
                Real fuel_frac = fire_arr(i,j,k,FireVar::fuel_frac);
                
                if (fuel_frac > 0.01) {
                    Real fuel_burned = fuel_burn_rate * dt;
                    fire_arr(i,j,k,FireVar::fuel_frac) = amrex::max(0.0, fuel_frac - fuel_burned/(m_fuel_density*m_fuel_depth));
                    fire_arr(i,j,k,FireVar::fuel_time) += dt;
                }
            }
        });
    }
}

void FireLevelSet::ComputeHeatRelease()
{
    BL_PROFILE("FireLevelSet::ComputeHeatRelease");
    
    const Real heat_per_area = m_fuel_density * m_fuel_depth * m_fuel_heat;
    const Real burn_duration = 60.0;
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> fire_arr = m_fire_data.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real phi = fire_arr(i,j,k,FireVar::phi);
            Real fuel_time = fire_arr(i,j,k,FireVar::fuel_time);
            
            if (phi < 0.0 && fuel_time < burn_duration && fuel_time > 0.0) {
                fire_arr(i,j,k,FireVar::heat_flux) = heat_per_area / burn_duration;
            } else {
                fire_arr(i,j,k,FireVar::heat_flux) = 0.0;
            }
        });
    }
}

void FireLevelSet::GetHeatFlux(MultiFab& heat_flux)
{
    MultiFab::Copy(heat_flux, m_fire_data, FireVar::heat_flux, 0, 1, 0);
}

void FireLevelSet::GetMoistureFlux(MultiFab& moisture_flux)
{
    moisture_flux.setVal(0.0);
}

Real FireLevelSet::get_fire_area() const
{
    if (!m_active) return 0.0;
    
    const Real* dx = m_geom->CellSize();
    const Real cell_area = dx[0] * dx[1];
    
    Real local_area = 0.0;
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> fire_arr = m_fire_data.const_array(mfi);
        
        Real box_area = 0.0;
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            Real phi = fire_arr(i,j,k,FireVar::phi);
            if (phi < 0.0) {
                box_area += cell_area;
            }
        }}}
        local_area += box_area;
    }
    
    ParallelDescriptor::ReduceRealSum(local_area);
    return local_area;
}

void FireLevelSet::GetPlotVarNames(Vector<std::string>& names) const
{
    names.push_back("fire_phi");
    names.push_back("fire_fuel_frac");
    names.push_back("fire_heat_flux");
    names.push_back("fire_fuel_time");
}

void FireLevelSet::GetPlotVar(const std::string& name,
                              MultiFab& mf,
                              int lev) const
{
    if (!m_active) {
        mf.setVal(0.0);
        return;
    }
    
    int fire_comp = -1;
    if (name == "fire_phi") fire_comp = FireVar::phi;
    else if (name == "fire_fuel_frac") fire_comp = FireVar::fuel_frac;
    else if (name == "fire_heat_flux") fire_comp = FireVar::heat_flux;
    else if (name == "fire_fuel_time") fire_comp = FireVar::fuel_time;
    else {
        mf.setVal(0.0);
        return;
    }
    
    const BoxArray& ba3d = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();
    
    BoxArray ba2d(ba3d.size());
    for (int i = 0; i < ba3d.size(); ++i) {
        Box bx = ba3d[i];
        bx.setSmall(2, 0);
        bx.setBig(2, 0);
        ba2d.set(i, bx);
    }

    MultiFab mf_2d(ba2d, dm, 1, 0);
    mf_2d.setVal(0.0);
    MultiFab::Copy(mf_2d, m_fire_data, fire_comp, 0, 1, 0);

    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx3d = mfi.tilebox();
        const auto& fire_2d = mf_2d.const_array(mfi);
        const auto& mf_arr = mf.array(mfi);
        
        ParallelFor(bx3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mf_arr(i,j,k) = fire_2d(i,j,0);
        });
    }
}