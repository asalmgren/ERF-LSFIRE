#include "ERF_FireLevelSet.H"
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <cmath>

using namespace amrex;

void FireLevelSet::Define(SolverChoice& sc)
{
    // Get met source from SolverChoice
    m_met_source = sc.fire_met_source;
    
    amrex::ParmParse pp("prob");
    
    // Fuel properties
    pp.query("fuel_density", m_fuel_density);
    pp.query("fuel_moisture", m_fuel_moisture);
    pp.query("fuel_heat", m_fuel_heat);
    pp.query("fuel_depth", m_fuel_depth);
    pp.query("fuel_savr", m_fuel_savr);
    
    // Fire spread
    pp.query("fire_ros_base", m_ros_base);
    pp.query("fire_wind_factor", m_wind_factor);
    pp.query("fire_slope_factor", m_slope_factor);
    pp.query("fire_refine_distance", m_refine_distance);
    
    // Ignition 
    pp.query("fire_ignition_time", m_ignition_time);
    pp.query("fire_ignition_x", m_ignition_x);
    pp.query("fire_ignition_y", m_ignition_y);
    pp.query("fire_ignition_radius", m_ignition_radius);
    
    // Wind source (for constant wind case)
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
    
    // ========== FIX: Keep 2D with proper indexing ==========
    const BoxArray& ba3d = cons_in.boxArray();
    BoxArray ba2d(ba3d.size());
    for (int i = 0; i < ba3d.size(); ++i) {
        Box bx = ba3d[i];
        // Flatten to single vertical layer at k=0
        bx.setSmall(2, 0);
        bx.setBig(2, 0);
        ba2d.set(i, bx);
    }
    
    const DistributionMapping& dm = cons_in.DistributionMap();
    IntVect nghost(1, 1, 0);  // Only horizontal ghost cells
    
    // Allocate fire state variables as 2D
    m_fire_data.define(ba2d, dm, FireVar::NumVars, nghost);
    // ========== END FIX ==========
    
    // Initialize all components
    m_fire_data.setVal(0.0);
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> fire_arr = m_fire_data.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // All at k=0
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
                          const MultiFab* U,  // ← Changed to pointer
                          const MultiFab* V,  // ← Changed to pointer
                          const MultiFab* terrain)  // ← Added terrain pointer
{
    // Store pointers (can be nullptr if using constant wind)
    m_wind_u = U;
    m_wind_v = V;
    m_terrain = terrain;

    // Check if we need to ignite
    static Real accumulated_time = 0.0;
    accumulated_time += dt;
    
    if (!m_active && accumulated_time >= m_ignition_time) {
        Vector<Real> ignition_loc(2);
        ignition_loc[0] = m_ignition_x;
        ignition_loc[1] = m_ignition_y;
        
        SetIgnition(ignition_loc, m_ignition_radius);
        m_active = true;
        
        if (ParallelDescriptor::IOProcessor()) {
            Print() << "\n🔥 Fire ignited at t=" << accumulated_time 
                   << " at (" << m_ignition_x << ", " << m_ignition_y << ")" 
                   << " with radius " << m_ignition_radius << " m\n\n";
        }
    }
    
    if (!m_active) return;
    
    BL_PROFILE("FireLevelSet::Advance");
    
    // ========== ADD: Print wind source info ==========
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "\n=== Fire Advance ===\n";
        Print() << "  dt = " << dt << " s\n";
        Print() << "  Met source: " << m_met_source << "\n";
        
        if (m_met_source == "Constant") {
            Print() << "  Using constant wind: U=" << m_constant_u 
                   << ", V=" << m_constant_v << " m/s\n";
        } else if (U && V) {
            Real u_max = U->max(0);
            Real u_min = U->min(0);
            Real v_max = V->max(0);
            Real v_min = V->min(0);
            Print() << "  Atmospheric wind U: [" << u_min << ", " << u_max << "] m/s\n";
            Print() << "  Atmospheric wind V: [" << v_min << ", " << v_max << "] m/s\n";
        }
    }
    // ========== END ==========
    
    // Fill ghost cells before computations
    m_fire_data.FillBoundary(m_geom->periodicity());
    
    // Step 1: Advect level set function
    AdvectLevelSet(dt);
    
    // Step 2: Reinitialize level set to maintain signed distance property
    ReinitializeLevelSet();
    
    // Step 3: Burn fuel where fire front passes
    BurnFuel(dt);
    
    // Step 4: Compute heat release
    ComputeHeatRelease();
    
    // Fill ghost cells after all updates
    m_fire_data.FillBoundary(m_geom->periodicity());
    
    // ========== ADD: Print fire area ==========
    if (ParallelDescriptor::IOProcessor()) {
        Real fire_area = get_fire_area();
        Print() << "  Fire area: " << fire_area << " m²\n";
        Print() << "====================\n\n";
    }
    // ========== END ==========
}

void FireLevelSet::SetIgnition(const Vector<Real>& location, const Real& radius)
{
    const Real* dx = m_geom->CellSize();
    const Real* prob_lo = m_geom->ProbLo();
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> fire_arr = m_fire_data.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // Compute cell center coordinates
            Real x = prob_lo[0] + (i + 0.5) * dx[0];
            Real y = prob_lo[1] + (j + 0.5) * dx[1];
            
            // Distance from ignition point
            Real dist = std::sqrt((x - location[0])*(x - location[0]) + 
                                 (y - location[1])*(y - location[1]));
            
            // Set level set: negative inside fire, positive outside
            fire_arr(i,j,k,FireVar::phi) = dist - radius;
            
            // Mark as burning if inside ignition zone
            if (dist < radius) {
                fire_arr(i,j,k,FireVar::fuel_time) = 0.0;
            }
        });
    }
    
    // Fill ghost cells after ignition
    m_fire_data.FillBoundary(m_geom->periodicity());
}

void FireLevelSet::AdvectLevelSet(const Real& dt)
{
    BL_PROFILE("FireLevelSet::AdvectLevelSet");
    
    const Real* dx = m_geom->CellSize();
    const Real ros_base = m_ros_base;
    const Real wind_factor = m_wind_factor;

    // Determine if using constant wind
    const bool use_constant_wind = (m_met_source == "Constant");
    const Real const_u = m_constant_u;
    const Real const_v = m_constant_v;
    
    m_fire_data.FillBoundary(m_geom->periodicity());
    
    // Create temporary for new phi values
    MultiFab phi_new(m_fire_data.boxArray(), m_fire_data.DistributionMap(), 1, 0);
    MultiFab::Copy(phi_new, m_fire_data, FireVar::phi, 0, 1, 0);
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        
        Array4<Real const> phi_old = m_fire_data.const_array(mfi);
        Array4<Real> phi = phi_new.array(mfi);
        
        // ========== DECLARE WIND ARRAYS (may be unused if constant) ==========
        Array4<Real const> u_arr, v_arr;
        int klo_u = 0, klo_v = 0;
        
        if (!use_constant_wind && m_wind_u && m_wind_v) {
            // Get atmospheric wind arrays
            u_arr = m_wind_u->const_array(mfi);
            v_arr = m_wind_v->const_array(mfi);
            
            // Get lowest vertical index for face-centered data
            Box ubx = mfi.tilebox(IntVect(1,0,0)); // x-face
            Box vbx = mfi.tilebox(IntVect(0,1,0)); // y-face
            klo_u = ubx.smallEnd(2);
            klo_v = vbx.smallEnd(2);
        }
        // ========== END ==========
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // Compute horizontal gradients at k=0
            Real phi_x = (phi_old(i+1,j,0,FireVar::phi) - phi_old(i-1,j,0,FireVar::phi)) / (2.0*dx[0]);
            Real phi_y = (phi_old(i,j+1,0,FireVar::phi) - phi_old(i,j-1,0,FireVar::phi)) / (2.0*dx[1]);
            
            Real grad_mag = std::sqrt(phi_x*phi_x + phi_y*phi_y);
            
            if (grad_mag > 1.0e-10) {
                // ========== GET WIND BASED ON SOURCE ==========
                Real u_wind, v_wind;
                
                if (use_constant_wind) {
                    // Use constant winds
                    u_wind = const_u;
                    v_wind = const_v;
                } else {
                    // Use atmospheric winds from lowest level
                    // Average from face-centered to cell-centered
                    u_wind = 0.5 * (u_arr(i,j,klo_u) + u_arr(i+1,j,klo_u));
                    v_wind = 0.5 * (v_arr(i,j,klo_v) + v_arr(i,j+1,klo_v));
                }
                // ========== END ==========
                
                // Normal vector (points outward from fire)
                Real nx = -phi_x / grad_mag;
                Real ny = -phi_y / grad_mag;
                
                // Wind component in normal direction (headwind positive)
                Real wind_dot = (u_wind * nx + v_wind * ny);
                
                // Rate of spread (wind increases ROS in direction of wind)
                Real ros = ros_base * (1.0 + wind_factor * wind_dot);
                
                // Ensure minimum ROS (fire always spreads a bit)
                ros = amrex::max(0.1 * ros_base, ros);
                
                // Update level set (phi decreases where fire spreads)
                phi(i,j,k,0) = phi_old(i,j,0,FireVar::phi) - dt * ros * grad_mag;
            }
        });
    }
    
    // Copy updated phi back
    MultiFab::Copy(m_fire_data, phi_new, 0, FireVar::phi, 1, 0);
    m_fire_data.FillBoundary(m_geom->periodicity());
}

void FireLevelSet::ReinitializeLevelSet()
{
    BL_PROFILE("FireLevelSet::ReinitializeLevelSet");
    
    // Simplified reinitialization: just ensure smoothness
    // Full implementation would use PDE-based reinitialization
    // For now, we skip this (will be refined in later phases)
}

void FireLevelSet::BurnFuel(const Real& dt)
{
    BL_PROFILE("FireLevelSet::BurnFuel");
    
    const Real fuel_burn_rate = m_fuel_density * m_fuel_depth / 60.0;  // kg/m²/s (1 min burn time)
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> fire_arr = m_fire_data.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real phi = fire_arr(i,j,k,FireVar::phi);
            
            // If fire front is at this location (phi near zero and negative)
            if (phi < 0.0 && phi > -10.0) {  // Within 10m of front
                Real fuel_frac = fire_arr(i,j,k,FireVar::fuel_frac);
                
                if (fuel_frac > 0.01) {  // Still has fuel
                    // Burn fuel
                    Real fuel_burned = fuel_burn_rate * dt;
                    fire_arr(i,j,k,FireVar::fuel_frac) = amrex::max(0.0, fuel_frac - fuel_burned/(m_fuel_density*m_fuel_depth));
                    
                    // Update burn time
                    fire_arr(i,j,k,FireVar::fuel_time) += dt;
                }
            }
        });
    }
}

void FireLevelSet::ComputeHeatRelease()
{
    BL_PROFILE("FireLevelSet::ComputeHeatRelease");
    
    const Real heat_per_area = m_fuel_density * m_fuel_depth * m_fuel_heat;  // J/m²
    const Real burn_duration = 60.0;  // seconds
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> fire_arr = m_fire_data.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real phi = fire_arr(i,j,k,FireVar::phi);
            Real fuel_time = fire_arr(i,j,k,FireVar::fuel_time);
            
            // Heat flux only where actively burning
            if (phi < 0.0 && fuel_time < burn_duration && fuel_time > 0.0) {
                // Peak heat release rate (W/m²)
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
    // Not implemented yet - will add in Phase 2
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
        
        // Simple loop instead of Reduce (compatible with all AMReX versions)
        Real box_area = 0.0;
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            Real phi = fire_arr(i,j,k,FireVar::phi);
            if (phi < 0.0) {  // Inside fire region
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
    
    if (name == "fire_phi") {
        MultiFab::Copy(mf, m_fire_data, FireVar::phi, 0, 1, 0);
    }
    else if (name == "fire_fuel_frac") {
        MultiFab::Copy(mf, m_fire_data, FireVar::fuel_frac, 0, 1, 0);
    }
    else if (name == "fire_heat_flux") {
        MultiFab::Copy(mf, m_fire_data, FireVar::heat_flux, 0, 1, 0);
    }
    else if (name == "fire_fuel_time") {
        MultiFab::Copy(mf, m_fire_data, FireVar::fuel_time, 0, 1, 0);
    }
}

void FireLevelSet::ComputeTerrainSlopes(const amrex::MultiFab* terrain, const amrex::Real* dx)
{
    BL_PROFILE("FireLevelSet::ComputeTerrainSlopes");
    
    // Check if terrain data is available
    if (!terrain) {
        // No terrain, set slopes to zero
        return;
    }
    
    // Create MultiFabs to store slopes
    MultiFab slope_x(m_fire_data.boxArray(), m_fire_data.DistributionMap(), 1, 0);
    MultiFab slope_y(m_fire_data.boxArray(), m_fire_data.DistributionMap(), 1, 0);
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        
        Array4<Real const> terrain_arr = terrain->const_array(mfi);
        Array4<Real> slope_x_arr = slope_x.array(mfi);
        Array4<Real> slope_y_arr = slope_y.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // Compute slopes using central differences
            // dz/dx
            slope_x_arr(i,j,k) = (terrain_arr(i+1,j,0) - terrain_arr(i-1,j,0)) / (2.0 * dx[0]);
            
            // dz/dy
            slope_y_arr(i,j,k) = (terrain_arr(i,j+1,0) - terrain_arr(i,j-1,0)) / (2.0 * dx[1]);
        });
    }
    
    // Store slopes for use in ROS calculation
    // You can add these to m_fire_data or store separately
}