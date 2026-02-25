#include "ERF_FireLevelSet.H"
#include "ERF_FireDataGatherer.H"
#include <AMReX_ParmParse.H>

using namespace amrex;

std::map<int, FuelModelData> FireLevelSet::m_fuel_database;
bool FireLevelSet::m_fuel_database_loaded = false;

FireLevelSet::~FireLevelSet() = default;

// ============================================================================
// INITIALIZATION AND SETUP
// ============================================================================

void FireLevelSet::Define(SolverChoice& sc) {
    ParmParse pp("prob");
    
    LoadFuelDatabase();
    
    pp.query("bandthickness", m_bandthickness);
    pp.query("target_cfl", m_target_cfl);
    pp.query("simulation_dt_max", m_simulation_dt_max);
    pp.query("untag_interval", m_untag_interval);
    pp.query("fire_verbosity", m_verbosity);
    
    pp.query("fire_ignition_time", m_ignition_time);
    pp.query("fire_ignition_x", m_ignition_x);
    pp.query("fire_ignition_y", m_ignition_y);
    pp.query("fire_ignition_radius", m_ignition_radius);
    
    std::string met_source_str = "Constant";
    pp.query("fire_met_source", met_source_str);
    if (met_source_str == "Atmosphere") {
        m_met_source = "Atmosphere";
    }
    pp.query("fire_constant_u", m_constant_u);
    pp.query("fire_constant_v", m_constant_v);
    pp.query("fuel_moisture", m_fuel_moisture);
    
    pp.query("fuel_model", m_default_fuel_model);
    
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "Fire model: Level-set method (ELMFIRE-based)\n";
        Print() << "  Band thickness: " << m_bandthickness << "\n";
        Print() << "  Target CFL: " << m_target_cfl << "\n";
        Print() << "  Max dt: " << m_simulation_dt_max << " s\n";
        Print() << "  Ignition: t=" << m_ignition_time << "s, ("
                << m_ignition_x << ", " << m_ignition_y << ") m, r=" 
                << m_ignition_radius << " m\n";
        Print() << "  Wind: U=" << m_constant_u << " m/s, V=" << m_constant_v << " m/s\n";
        Print() << "  Fuel model: " << m_default_fuel_model;
        if (m_fuel_database.count(m_default_fuel_model)) {
            Print() << " (" << m_fuel_database[m_default_fuel_model].name << ")";
        }
        Print() << "\n";
        Print() << "  🔥 Fire module: DISTRIBUTED with ghost cells\n";
    }
}

void FireLevelSet::Init(int lev, const MultiFab& cons_in,
                        const Geometry& geom, const Real& dt_advance) {

    m_lev = lev;
    m_geom = &geom;

    const BoxArray& ba = cons_in.boxArray();
    const DistributionMapping& dm = cons_in.DistributionMap();
    
    // ✅ Calculate ghost cells based on band thickness
    const Real* dx = geom.CellSize();
    Real band_meters = (m_bandthickness + 1) * std::max(dx[0], dx[1]);
    int nghost = m_bandthickness + 5;  // +3 for gradient stencil + safety
    
    // ✅ Create distributed MultiFab WITH GHOST CELLS
    m_fire_data.define(ba, dm, FireVar::NumVars, nghost);
    m_fire_data.setVal(1000.0, FireVar::phi, 1, nghost);
    m_fire_data.setVal(1.0, FireVar::fuel_frac, 1, nghost);
    m_fire_data.setVal(0.0, FireVar::fuel_time, 1, nghost);
    m_fire_data.setVal(0.0, FireVar::heat_flux, 1, nghost);
    m_fire_data.setVal(0.0, FireVar::burned_area, 1, nghost);
    m_fire_data.FillBoundary(m_geom->periodicity());
    
    m_slope_x = std::make_unique<MultiFab>(ba, dm, 1, nghost);
    m_slope_y = std::make_unique<MultiFab>(ba, dm, 1, nghost);
    m_slope_x->setVal(0.0);
    m_slope_y->setVal(0.0);
    
    // ✅ Each rank has LOCAL lists
    m_local_tagged.clear();
    m_local_burned.clear();
    
    // ✅ Create FireDataGatherer
    m_data_gatherer = std::make_unique<FireDataGatherer>(m_geom);
    m_data_gatherer->SetConstantWind(m_constant_u, m_constant_v);
    
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "🔥 FireLevelSet initialized:\n";
        Print() << "   nghost = " << nghost << " cells\n";
        Print() << "   band_thickness = " << band_meters << " m\n";
        Print() << "   Each rank manages its local fire cells\n";
    }
    
    ParallelDescriptor::Barrier();
}

// ============================================================================
// IGNITION SETUP (ALL RANKS PARTICIPATE)
// ============================================================================

// Lines 93-147: Replace SetIgnition function

void FireLevelSet::SetIgnition(const Vector<Real>& location, const Real& radius) {
    BL_PROFILE("FireLevelSet::SetIgnition");
    
    const Real* dx = m_geom->CellSize();
    const Real* prob_lo = m_geom->ProbLo();

    int ic = static_cast<int>((location[0] - prob_lo[0]) / dx[0]);
    int jc = static_cast<int>((location[1] - prob_lo[1]) / dx[1]);
    
    // Define narrow band width (should be > 5*dx to avoid gradient issues!)
    Real narrow_band_width = 30.0 * std::max(dx[0], dx[1]);  // ← INCREASED from 5 to 10!
    
    int local_ignited_count = 0;
    
    // ========================================
    // STEP 1: ALL RANKS initialize phi in their boxes
    // ========================================
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Box& grown_box = mfi.growntilebox();  // ✅ Include ghost cells!
        const Real* prob_lo = m_geom->ProbLo();
        
        Array4<Real> phi = m_fire_data.array(mfi);
        Array4<Real> fuel_time = m_fire_data.array(mfi, FireVar::fuel_time);
        
        int k = 0;  // Fire is 2D
        
        // ✅ Loop over GROWN box (includes ghost cells for smooth gradients!)
        for (int j = grown_box.smallEnd(1); j <= grown_box.bigEnd(1); j++) {
            for (int i = grown_box.smallEnd(0); i <= grown_box.bigEnd(0); i++) {
                
                // Compute global coordinates
                Real x = prob_lo[0] + (i + 0.5) * dx[0];
                Real y = prob_lo[1] + (j + 0.5) * dx[1];
                
                // Distance from ignition center
                Real dist = std::sqrt((x - location[0]) * (x - location[0]) +
                                     (y - location[1]) * (y - location[1]));
                
                // ========================================
                // ✅ CRITICAL: Proper signed distance convention
                // ========================================
                // φ < 0: INSIDE burned region (already burned)
                // φ = 0: AT fire front
                // φ > 0: OUTSIDE (unburned fuel)
                
                Real phi_val = dist - radius;
                phi(i, j, k ) = phi_val;

                if (phi_val < 0.0) {
                    fuel_time(i, j, k) = m_current_time;
                    local_ignited_count++;
                }
                
                // if (dist < radius) {
                //     // Inside ignition circle: φ < 0 (burned)
                //     phi_val = dist - radius;  // Negative value
                    
                //     // Only set fuel_time in valid box (not ghost cells!)
                //     if (box.contains(IntVect(i, j, k))) {
                //         fuel_time(i, j, k) = m_current_time;
                //         local_ignited_count++;
                //     }
                    
                // } else if (dist < radius + narrow_band_width) {
                //     // Outside but within narrow band: φ > 0 (smooth transition)
                //     phi_val = dist - radius;  // Positive value
                    
                // } else {
                //     // Far outside: constant large value
                //     phi_val = narrow_band_width;  // ← NOT 100.0!
                // }
                
                // phi(i, j, k) = phi_val;
            }
        }
    }
    
    // ========================================
    // STEP 2: Exchange ghost cells between ranks
    // ========================================
    m_fire_data.FillBoundary(m_geom->periodicity());
    
    // ========================================
    // STEP 3: Global ignition count
    // ========================================
    int global_ignited_count = local_ignited_count;
#ifdef AMREX_USE_MPI
    ParallelDescriptor::ReduceIntSum(global_ignited_count);
#endif
    
    // ========================================
    // STEP 4: Diagnostics
    // ========================================
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "🔥 Ignition set at (" << location[0] << ", " << location[1] 
                << ") m, radius=" << radius << " m\n";
        // Print() << "   Cell size: " << dx[0] << " × " << dx[1] << " m\n";
        // Print() << "   Narrow band width: " << narrow_band_width << " m\n";
        Print() << "   Expected cells: ~" << M_PI * std::pow(radius/std::max(dx[0],dx[1]), 2.0) << "\n";
        Print() << "   Ignited " << global_ignited_count << " cells\n";
        Print() << "   phi convention: negative inside (burned), positive outside (unburned)\n";
    }
    
    m_needs_initial_tagging = true;

    TagBandLocal(ic, jc, m_current_time);

    int local_tagged = m_local_tagged.size();
    int global_tagged = local_tagged;
#ifdef AMREX_USE_MPI
    ParallelDescriptor::ReduceIntSum(global_tagged);
#endif
        
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "  Tagged " << global_tagged << " cells total (distributed)\n";
    }
}

// ============================================================================
// MAIN TIME ADVANCEMENT (DISTRIBUTED WITH GHOST CELLS)
// ============================================================================

void FireLevelSet::Advance(const Real& dt, const MultiFab* U, 
                           const MultiFab* V, const MultiFab* terrain,
                           const iMultiFab* luindex) {    
    BL_PROFILE("FireLevelSet::Advance");
    
    m_wind_u = U;
    m_wind_v = V;
    m_terrain = terrain;
    m_luindex = luindex;
    
    // Compute slopes (once)
    static bool slopes_computed = false;
    if (!slopes_computed && terrain) {
        const Real* dx = m_geom->CellSize();
        ComputeTerrainSlopes(terrain, dx);
        slopes_computed = true;
        
        if (ParallelDescriptor::IOProcessor()) {
            Print() << "  ⛰️  Terrain slopes computed\n";
        }
    }

    static bool first_advance = true;
    if (first_advance && ParallelDescriptor::IOProcessor()) {
        Print() << "  🔍 Advance() MultiFab pointers:\n";
        Print() << "     m_wind_u: " << (m_wind_u ? "✅ set" : "❌ null") << "\n";
        Print() << "     m_wind_v: " << (m_wind_v ? "✅ set" : "❌ null") << "\n";
        Print() << "     m_terrain: " << (m_terrain ? "✅ set" : "❌ null") << "\n";
        Print() << "     m_luindex: " << (m_luindex ? "✅ set" : "❌ null") << "\n";
        first_advance = false;
    }
    
    // ========================================
    // STEP 1: Check for ignition
    // ========================================
    static bool ignition_set = false;
    if (!ignition_set && m_current_time >= m_ignition_time) {
        Vector<Real> location = {m_ignition_x, m_ignition_y, 0.0};
        SetIgnition(location, m_ignition_radius);
        ignition_set = true;
    }
    
    // ========================================
    // STEP 2: Initial tagging (EACH RANK)
    // ========================================
    // Moved to SetIgnition()
    
    // ========================================
    // STEP 3: Check if fire exists
    // ========================================
    int local_has_fire = m_local_tagged.empty() ? 0 : 1;
    int global_has_fire = local_has_fire;
#ifdef AMREX_USE_MPI
    ParallelDescriptor::ReduceIntMax(global_has_fire);
#endif
    
    if (!global_has_fire) {
        m_current_time += dt;
        return;
    }
    
    // ========================================
    // STEP 4: Gather wind data (ALL RANKS)
    // ========================================
    // No longer needed
    
    // ========================================
    // STEP 5: Update weather (EACH RANK LOCAL)
    // ========================================
    // if (ParallelDescriptor::IOProcessor()) {
    //     Print() << "  📍 About to call UpdateNodeWeatherDataLocal()\n";
    //     Print() << "     m_wind_u: " << (m_wind_u ? "✅ set" : "❌ null") << "\n";
    //     Print() << "     m_local_tagged.size(): " << m_local_tagged.size() << "\n";
    // }
    
    UpdateNodeWeatherDataLocal();
    
    // ✅ DIAGNOSTIC: Check if nodes got data
    if (!m_local_tagged.empty() && ParallelDescriptor::IOProcessor()) {
        const auto& sample = m_local_tagged[0];
        Print() << "  📊 After UpdateNodeWeatherDataLocal():\n";
        Print() << "     Sample node (" << sample.i << "," << sample.j << "):\n";
        Print() << "       u_wind: " << sample.u_wind << " m/s\n";
        Print() << "       v_wind: " << sample.v_wind << " m/s\n";
        Print() << "       ws20: " << sample.ws20 << " m/s\n";
        Print() << "       wd20: " << sample.wd20 << "°\n";
        Print() << "       fuel: " << sample.ifbfm << "\n";
    }
    
    // ========================================
    // STEP 6: Fire physics (EACH RANK LOCAL)
    // ========================================
    SurfaceSpreadRateLocal();

    if (!m_local_tagged.empty()) {
        const auto& sample = m_local_tagged[0];
        AllPrint() << "  ✅ VERIFY immediately after SurfaceSpreadRateLocal():\n";
        AllPrint() << "     Cell (" << sample.i << "," << sample.j << "):\n";
        AllPrint() << "     phiw_surface = " << sample.phiw_surface << " (expect ~0.002)\n";
        AllPrint() << "     phis_surface = " << sample.phis_surface << " (expect 0.0)\n";
        AllPrint() << "     velocity_dms_surface = " << sample.velocity_dms_surface << " ft/min\n";
}
    
    Real local_max_velocity = GetMaxFireVelocityLocal();
    
    // ✅ Global max velocity
    Real max_velocity = local_max_velocity;
#ifdef AMREX_USE_MPI
    ParallelDescriptor::ReduceRealMax(max_velocity);
#endif
    
    Real fire_dt = ComputeAdaptiveTimestep(max_velocity);
    fire_dt = std::min(fire_dt, dt);
    
    // ========================================
    // STEP 7: Logging (IO PROCESSOR)
    // ========================================
    static int step_counter = 0;
    if (++step_counter % 1 == 0) {
        int local_tagged = m_local_tagged.size();
        int local_burned = m_local_burned.size();
        
        // Each rank prints IF it has cells
        if (local_tagged > 0) {
            Real local_max_vel = 0.0;
            for (const auto& node : m_local_tagged) {
                Real vel = std::sqrt(node.ux*node.ux + node.uy*node.uy);
                local_max_vel = std::max(local_max_vel, vel);
            }
            
            Print() << "  [Rank " << ParallelDescriptor::MyProc() << "] "
                    << "t=" << m_current_time << "s, "
                    << "tagged=" << local_tagged << ", "
                    << "max_vel=" << local_max_vel << " m/s, "
                    << "wind: u=" << m_local_tagged.front().u_wind 
                    << ", v=" << m_local_tagged.front().v_wind << "\n";
        }
        
        // ParallelDescriptor::Barrier();  // Sync output
        
        // Global summary
        int global_tagged = local_tagged;
        int global_burned = local_burned;
#ifdef AMREX_USE_MPI
        ParallelDescriptor::ReduceIntSum(global_tagged);
        ParallelDescriptor::ReduceIntSum(global_burned);
#endif
        
        if (ParallelDescriptor::IOProcessor()) {
            Print() << "  📊 GLOBAL SUMMARY: t=" << m_current_time 
                    << "s, tagged=" << global_tagged 
                    << ", burned=" << global_burned << "\n";
        }
    }
    
    // ========================================
    // STEP 8: RK2 integration (EACH RANK)
    // ========================================
    for (int step = 1; step <= 2; step++) {
        CalcNormalVectorsLocal(step);
        UxAndUyEllipticalLocal(step, m_current_time);     
        ApplyFluxLimiterLocal();
        Rk2IntegrateLocal(fire_dt, step);
    }

    // ========================================
    // STEP 9: Update burned cells (EACH RANK)
    // ========================================
    UpdateBurnedCellsLocal(m_current_time);
    
    // ========================================
    // STEP 10: Reinitialization (if needed)
    // ========================================
    static Real distance_since_reinit = 0.0;
    static int steps_since_reinit = 0;
    
    distance_since_reinit += max_velocity * fire_dt;
    steps_since_reinit++;
    
    // ✅ Each rank checks gradient degradation
    Real local_sum_gradmag = 0.0;
    int local_count = 0;
    for (const auto& node : m_local_tagged) {
        Real gradmag = std::sqrt(node.dphidx_limited * node.dphidx_limited + 
                                node.dphidy_limited * node.dphidy_limited);
        local_sum_gradmag += gradmag;
        local_count++;
    }
    
    Real global_sum_gradmag = local_sum_gradmag;
    int global_count = local_count;
#ifdef AMREX_USE_MPI
    ParallelDescriptor::ReduceRealSum(global_sum_gradmag);
    ParallelDescriptor::ReduceIntSum(global_count);
#endif
    
    Real avg_gradmag = (global_count > 0) ? global_sum_gradmag / global_count : 1.0;
    
    const Real* dx = m_geom->CellSize();
    Real reinit_distance = 3.0 * std::max(dx[0], dx[1]);
    bool gradients_degraded = std::abs(avg_gradmag - 1.0) > 0.3;
    bool traveled_far = distance_since_reinit > reinit_distance;

    int local_size = m_local_burned.size();
    int global_size = local_size;
    #ifdef AMREX_USE_MPI
    ParallelDescriptor::ReduceIntSum(global_size);
    #endif    
    
    if ((gradients_degraded || traveled_far) && steps_since_reinit > 50) {
        if (ParallelDescriptor::IOProcessor()) {
            Print() << "  🔄 Reinit triggered (|∇φ|=" << avg_gradmag
                    << "     m_local_burned.size() = " << global_size << " (global)\n"
                    << ", dist=" << distance_since_reinit << "m)\n";
        }
        ReinitializeSignedDistanceDistributed();
        distance_since_reinit = 0.0;
        steps_since_reinit = 0;
    }

    // ========================================
    // STEP 10.5: Clean up old burned cells (keep last 5 minutes only)
    // ========================================
    // Real cleanup_threshold = m_current_time - 300.0;  // 5 minutes
    // m_local_burned.erase(
    //     std::remove_if(m_local_burned.begin(), m_local_burned.end(),
    //         [cleanup_threshold](const FireNode& node) {
    //             return node.time_added < cleanup_threshold;
    //         }),
    //     m_local_burned.end()
    // );    
    
    BurnFuelLocal(fire_dt);
    ComputeHeatReleaseLocal();

    // ========================================
    // STEP 11: Untag cells (EACH RANK)
    // ========================================
    static int untag_counter = 0;
    if (++untag_counter >= m_untag_interval) {
        UntagCellsLocal(m_current_time);
        untag_counter = 0;
    }
    
    // ========================================
    // STEP 12: Burn fuel & heat release
    // ========================================
    
    m_current_time += dt;
    
    // ParallelDescriptor::Barrier();
}

// ============================================================================
// LOCAL OPERATIONS (EACH RANK ON ITS BOXES)
// ============================================================================

void FireLevelSet::TagBandLocal(int i_center, int j_center, Real time) {
    const Box& domain = m_geom->Domain();
    const Real* dx = m_geom->CellSize();
    
    Real cell_size = std::max(dx[0], dx[1]);
    Real max_tag_distance = (m_bandthickness + 1) * cell_size;
    
    int local_tagged = 0;
    
    // ✅ Search in grown box (includes ghost cells!)
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& validbox = mfi.validbox();
        const Box& grownbox = mfi.growntilebox();  // ✅ Includes ghosts!
        
        Array4<Real const> phi = m_fire_data.const_array(mfi);
        
        for (int j = j_center - m_bandthickness; j <= j_center + m_bandthickness; j++) {
            for (int i = i_center - m_bandthickness; i <= i_center + m_bandthickness; i++) {
                // Check if in grown box (can see via ghost cells!)
                if (!grownbox.contains(IntVect(i, j, 0))) continue;
                if (!domain.contains(IntVect(i, j, 0))) continue;
                
                if (IsNonburnable(i, j)) continue;
                
                // Check if already tagged
                bool already_tagged = false;
                for (const auto& node : m_local_tagged) {
                    if (node.i == i && node.j == j) {
                        already_tagged = true;
                        break;
                    }
                }
                if (already_tagged) continue;  // Skip already tagged

                bool already_burned = false;
                for (const auto& node : m_local_burned) {
                    if (node.i == i && node.j == j) {
                        already_burned = true;
                        break;
                    }
                }
                if (already_burned) continue; // Skip already burned

                bool already_front = false;
                for (const auto& node : m_local_front) {
                    if (node.i == i && node.j == j) {
                        already_front = true;
                        break;
                    }
                }
                if (already_front) continue; // Skip already front                
                
                Real phi_val = phi(i, j, 0);
                
                // if (std::abs(phi_val) > max_tag_distance) continue;
                
                // Only tag if in THIS rank's valid box
                // if (!validbox.contains(IntVect(i, j, 0))) continue;
                if (!grownbox.contains(IntVect(i, j, 0))) continue;
                
                FireNode new_node;
                new_node.i = i;
                new_node.j = j;
                new_node.time_added = time;
                new_node.burned = false;
                new_node.phip_old = phi_val;
                new_node.ifbfm = m_default_fuel_model;
                
                m_local_tagged.push_back(new_node);
                local_tagged++;
            }
        }
    }
}

void FireLevelSet::UpdateNodeWeatherDataLocal() {
    BL_PROFILE("FireLevelSet::UpdateNodeWeatherDataLocal");
    
    // ✅ FIX: Use AllPrint() so ALL ranks print!
    AllPrint() << "  [Rank " << ParallelDescriptor::MyProc() 
               << "] UpdateNodeWeatherDataLocal() START: " 
               << m_local_tagged.size() << " cells\n";
    
    if (m_local_tagged.empty()) {
        AllPrint() << "  [Rank " << ParallelDescriptor::MyProc() 
                   << "] Early exit: no tagged cells\n";
        return;
    }
    
    const Real* dx = m_geom->CellSize();
    
    // ✅ All ranks print
    AllPrint() << "  [Rank " << ParallelDescriptor::MyProc() << "] MultiFab checks:\n";
    AllPrint() << "     m_wind_u: " << (m_wind_u ? "✅" : "❌") << "\n";
    AllPrint() << "     m_wind_v: " << (m_wind_v ? "✅" : "❌") << "\n";
    AllPrint() << "     m_slope_x: " << (m_slope_x ? "✅" : "❌") << "\n";
    AllPrint() << "     m_slope_y: " << (m_slope_y ? "✅" : "❌") << "\n";
    AllPrint() << "     m_luindex: " << (m_luindex ? "✅" : "❌") << "\n";
    
    int cells_with_wind = 0;
    int cells_with_slope = 0;
    int cells_with_fuel = 0;
    
    // ========================================
    // STEP 1: Read wind data
    // ========================================
    AllPrint() << "  [Rank " << ParallelDescriptor::MyProc() 
               << "] Reading wind data...\n";

    for (auto& node : m_local_tagged) {
        if (node.burned) continue;

        Real u_wind = m_constant_u;  // Default fallback
        Real v_wind = m_constant_v;

        if (m_met_source == "Atmosphere" && m_wind_u && m_wind_v) {
            for (MFIter mfi(*m_wind_u); mfi.isValid(); ++mfi) {
                const Box& grownbox = mfi.growntilebox();  // Includes ghosts!
                if (grownbox.contains(IntVect(node.i, node.j, 0))) {
                    Array4<Real const> u_arr = m_wind_u->const_array(mfi);
                    Array4<Real const> v_arr = m_wind_v->const_array(mfi);

                    // Read wind at k=0 (surface level)
                    u_wind = u_arr(node.i, node.j, 0);
                    v_wind = v_arr(node.i, node.j, 0);   
                    
                    break;
                }
            }
        }

        // Convert wind to fire spread parameters
        Real wind_speed = std::sqrt(u_wind*u_wind + v_wind*v_wind);
        
        Real wind_dir_to = std::atan2(u_wind, v_wind) * 180.0 / M_PI;
        Real wind_dir_from = wind_dir_to + 180.0;
        
        if (wind_dir_from >= 360.0) wind_dir_from -= 360.0;
        if (wind_dir_from < 0.0)    wind_dir_from += 360.0;
        
        node.ws20 = wind_speed;
        node.wd20 = wind_dir_from;
        
        // Adjust from 20m to 6.1m height
        node.ws20 *= std::pow(6.1 / 40.0, 0.11);
        
        // Apply fuel-specific wind reduction factor
        Real wind_adjustment = 0.3;  // Default for typical shrub/grass mix

        if (node.ifbfm >= 101 && node.ifbfm <= 103) {
            wind_adjustment = 0.42;  // ✅ GR1-GR3: Short grass (open, less drag)
        } else if (node.ifbfm >= 121 && node.ifbfm <= 124) {
            wind_adjustment = 0.38;  // ✅ GS1-GS4: Grass-shrub mix
        } else if (node.ifbfm >= 141 && node.ifbfm <= 149) {
            wind_adjustment = 0.27;  // ✅ SH1-SH9: Low-moderate shrub
        } else if (node.ifbfm >= 161 && node.ifbfm <= 165) {
            wind_adjustment = 0.32;  // ✅ SH1-SH5: High shrub (YOUR CASE!)
        } else if (node.ifbfm >= 181 && node.ifbfm <= 189) {
            wind_adjustment = 0.11;  // ✅ TU1-TU5: Timber understory (heavy sheltering)
        }
        
        node.wsmf = node.ws20 * wind_adjustment * 196.85;  // Convert to ft/min        

        // ✅ Get terrain slopes
        node.slope_x = 0.0;
        node.slope_y = 0.0;

        if (m_slope_x && m_slope_y) {
            for (MFIter mfi(*m_slope_x); mfi.isValid(); ++mfi) {
                const Box& grownbox = mfi.growntilebox();  // Includes ghosts!
                if (grownbox.contains(IntVect(node.i, node.j, 0))) {
                    Array4<Real const> sx_arr = m_slope_x->const_array(mfi);
                    Array4<Real const> sy_arr = m_slope_y->const_array(mfi);
                    
                    node.slope_x = sx_arr(node.i, node.j, 0);
                    node.slope_y = sy_arr(node.i, node.j, 0);
                    break;
                }
            }
        }

        if (node.i == 64 && node.j == 57) {
            AllPrint() << "🌬️  Cell (64,57) wind update:\n"
                        << "  U=" << u_wind << " m/s, V=" << v_wind << " m/s\n"
                        << "  wind_speed=" << wind_speed << " m/s\n"
                        << "  ws20=" << node.ws20 << " m/s (after height adj)\n"
                        << "  wind_adjustment=" << wind_adjustment << "\n"
                        << "  wsmf=" << node.wsmf << " ft/min\n";
        }        
        
        node.adj = 1.0;

    }

}

void FireLevelSet::UpdateBurnedCellsLocal(Real current_time) {
    const Real* dx = m_geom->CellSize();
    Real cell_size = std::min(dx[0], dx[1]);
    
    // Define thresholds
    Real front_threshold = 0.5 * cell_size;   // ±150m for 300m cells
    Real burn_threshold  = -0.5 * cell_size;  // -150m
    
    const Box& domain = m_geom->Domain();
    
    std::vector<FireNode> newly_burned;
    std::vector<FireNode> newly_front;
    
    // ========================================
    // STEP 1: Classify cells by φ value
    // ========================================
    auto it = m_local_tagged.begin();
    while (it != m_local_tagged.end()) {
        Real phi = it->phip_old;
        
        // Case 1: BURNED (φ < -150m)
        if (phi < burn_threshold) {
            if (!it->burned) {  // ← Only trigger if newly burned
                newly_burned.push_back(*it);
                it->burned = true;
                
                // Update MultiFab
                for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
                    const Box& validbox = mfi.validbox();
                    if (!validbox.contains(IntVect(it->i, it->j, 0))) continue;
                    
                    Array4<Real> data = m_fire_data.array(mfi);
                    if (data(it->i, it->j, 0, FireVar::fuel_time) == 0.0) {
                        data(it->i, it->j, 0, FireVar::fuel_time) = current_time;
                    }
                    data(it->i, it->j, 0, FireVar::burned_area) = 1.0;
                    break;
                }
                
                // Move to burned list and remove from tagged
                m_local_burned.push_back(*it);
                it = m_local_tagged.erase(it);
                continue;  // ← Skip ++it (erase returns next iterator)
            }
        }
        
        // Case 2: FRONT (-150m < φ < +150m)
        else if (phi >= burn_threshold && phi <= front_threshold) {
            if (!it->front) {  // ← Only trigger if newly at front
                newly_front.push_back(*it);
                it->front = true;
                
                // Update MultiFab (mark as active fire front)
                for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
                    const Box& validbox = mfi.validbox();
                    if (!validbox.contains(IntVect(it->i, it->j, 0))) continue;
                    
                    Array4<Real> data = m_fire_data.array(mfi);
                    data(it->i, it->j, 0, FireVar::burned_area) = 0.5;  // ← Mark as "burning"
                    break;
                }
                
                // Add to front list (but KEEP in tagged list!)
                m_local_front.push_back(*it);
            }
        }
        
        // Case 3: UNBURNED (φ > +150m) - do nothing
        
        ++it; 
    }
    
    // ========================================
    // STEP 2: Tag new cells around FRONT cells
    // ========================================
    if (!newly_front.empty()) {        
        for (const auto& node : newly_front) {
            TagBandLocal(node.i, node.j, current_time);
        }
    }
    
    // ========================================
    // STEP 3: Tag new cells around BURNED cells (backup)
    // ========================================
    if (!newly_burned.empty()) {
        for (const auto& node : newly_burned) {
            TagBandLocal(node.i, node.j, current_time);
        }
    }
}

void FireLevelSet::UntagCellsLocal(Real current_time) {
    BL_PROFILE("FireLevelSet::UntagCellsLocal");
    
    const Real* dx = m_geom->CellSize();
    const Box& domain = m_geom->Domain();
    Real max_tag_time = 604800.0;  // 1 week in seconds
    Real max_tag_distance = 1.0 * m_bandthickness * std::max(dx[0], dx[1]);

    
    int num_deleted = 0;
    
    // ✅ DON'T EARLY RETURN! All ranks must participate in MPI_Reduce!
    if (!m_local_tagged.empty()) {
        auto it = m_local_tagged.begin();
        while (it != m_local_tagged.end()) {
            int i = it->i;
            int j = it->j;
            bool should_untag = false;
            
            // Criterion 1: Cell tagged too long
            if (it->time_added > 0.0 && 
                (current_time - it->time_added) > max_tag_time) {
                should_untag = true;
            }
            
            // Criterion 2: Check if isolated
            if (!should_untag) {
                bool has_neighbor[4] = {false, false, false, false};
                
                for (const auto& other : m_local_tagged) {
                    if (other.i == i-1 && other.j == j) has_neighbor[0] = true;
                    if (other.i == i+1 && other.j == j) has_neighbor[1] = true;
                    if (other.i == i && other.j == j-1) has_neighbor[2] = true;
                    if (other.i == i && other.j == j+1) has_neighbor[3] = true;
                }
                
                if (!has_neighbor[0] && !has_neighbor[1] && 
                    !has_neighbor[2] && !has_neighbor[3]) {
                    should_untag = true;
                }
            }
            
            // Criterion 3: All neighbors in band are burned
            if (!should_untag) {
                int istart = std::max(domain.smallEnd(0), i - m_bandthickness);
                int istop  = std::min(domain.bigEnd(0),   i + m_bandthickness);
                int jstart = std::max(domain.smallEnd(1), j - m_bandthickness);
                int jstop  = std::min(domain.bigEnd(1),   j + m_bandthickness);
                
                bool all_neighbors_burned = true;
                
                for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
                    if (!mfi.validbox().contains(IntVect(i, j, 0))) continue;
                    
                    Array4<Real const> phi = m_fire_data.const_array(mfi, FireVar::phi);
                    
                    for (int ii = istart; ii <= istop; ii++) {
                        if (domain.contains(IntVect(ii, j, 0))) {
                            if (phi(ii, j, 0) >= 0.0) {
                                all_neighbors_burned = false;
                                break;
                            }
                        }
                    }
                    
                    if (all_neighbors_burned) {
                        for (int jj = jstart; jj <= jstop; jj++) {
                            if (domain.contains(IntVect(i, jj, 0))) {
                                if (phi(i, jj, 0) >= 0.0) {
                                    all_neighbors_burned = false;
                                    break;
                                }
                            }
                        }
                    }
                    
                    break;
                }
                
                if (all_neighbors_burned) {
                    should_untag = true;
                }
            }
            
            // Criterion 4: Too far away from fire front
            if (!should_untag) {
                for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox();
                    if (bx.contains(IntVect(i, j, 0))) {
                        Array4<Real const> phi = m_fire_data.const_array(mfi);
                        Real phi_val = phi(i, j, 0, FireVar::phi);
                        if (std::abs(phi_val) > max_tag_distance) {
                            should_untag = true;
                        }
                    }
                }
            }

            if (should_untag) {
                it = m_local_tagged.erase(it);
                num_deleted++;
            } else {
                ++it;
            }
        }
    }
    
    // ✅ ALL ranks MUST participate in this reduction!
    int global_deleted = num_deleted;
#ifdef AMREX_USE_MPI
    ParallelDescriptor::ReduceIntSum(global_deleted);
#endif
    
    // ✅ Safe to print (all ranks have completed MPI operation)
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "  🗑️  Untagged " << global_deleted << " cells\n";
    }
}

void FireLevelSet::BurnFuelLocal(const Real& dt) {
    const Real fuel_burn_rate = 0.01;
    
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> data = m_fire_data.array(mfi);
        
        for (const auto& node : m_local_burned) {
            if (!bx.contains(IntVect(node.i, node.j, 0))) continue;
            Real& fuel_frac = data(node.i, node.j, 0, FireVar::fuel_frac);
            fuel_frac = std::max(0.0, fuel_frac - fuel_burn_rate * dt);
        }
    }    
}

void FireLevelSet::ComputeHeatReleaseLocal() {
    // Heat release stored in node.ir
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> heat = m_fire_data.array(mfi, FireVar::heat_flux);
        
        for (const auto& node : m_local_burned) {
            if (bx.contains(IntVect(node.i, node.j, 0))) {
                heat(node.i, node.j, 0) = node.ir;
            }
        }
    }
}

Real FireLevelSet::GetMaxFireVelocityLocal() {
    Real max_vel = 0.0;
    for (const auto& node : m_local_tagged) {
        if (node.burned) continue;
        Real vel = std::sqrt(node.ux*node.ux + node.uy*node.uy);
        max_vel = std::max(max_vel, vel);
    }
    return max_vel;
}

// ============================================================================
// REINITIALIZATION (DISTRIBUTED WITH GATHER-BROADCAST)
// ============================================================================

// Replace ReinitializeSignedDistanceDistributed (lines 1130-1273)

void FireLevelSet::ReinitializeSignedDistanceDistributed() {
    BL_PROFILE("FireLevelSet::Reinit");
    
    const Real* dx = m_geom->CellSize();
    const Real* prob_lo = m_geom->ProbLo();
    
    // ========================================
    // STEP 1: Collect RECENT burned cells only (last 5 minutes = 300s)
    // ========================================
    std::vector<int> local_burned_i, local_burned_j;
    Real recent_time_threshold = m_current_time - 300.0;  // 5 minutes
    
    for (const auto& node : m_local_burned) {
        // ✅ Only include recently burned cells!
        if (node.time_added >= recent_time_threshold) {
            local_burned_i.push_back(node.i);
            local_burned_j.push_back(node.j);
        }
    }
    
    // If no recent burns, fall back to last 100 burned cells
    if (local_burned_i.empty() && !m_local_burned.empty()) {
        int start_idx = std::max(0, static_cast<int>(m_local_burned.size()) - 100);
        for (size_t n = start_idx; n < m_local_burned.size(); n++) {
            local_burned_i.push_back(m_local_burned[n].i);
            local_burned_j.push_back(m_local_burned[n].j);
        }
    }
    
    // ========================================
    // STEP 2: Gather to Rank 0
    // ========================================
    std::vector<int> all_burned_i, all_burned_j;
    int total_burned = 0;
    
#ifdef AMREX_USE_MPI
    int nranks = ParallelDescriptor::NProcs();
    std::vector<int> counts(nranks);
    
    int local_count = local_burned_i.size();
    
    // ✅ Safety check: limit to avoid huge MPI messages
    const int MAX_BURNED_PER_RANK = 10000;
    if (local_count > MAX_BURNED_PER_RANK) {
        if (ParallelDescriptor::IOProcessor()) {
            Print() << "  ⚠️  Rank has " << local_count << " burned cells, limiting to " 
                    << MAX_BURNED_PER_RANK << "\n";
        }
        local_burned_i.resize(MAX_BURNED_PER_RANK);
        local_burned_j.resize(MAX_BURNED_PER_RANK);
        local_count = MAX_BURNED_PER_RANK;
    }
    
    MPI_Gather(&local_count, 1, MPI_INT,
               counts.data(), 1, MPI_INT,
               0, ParallelDescriptor::Communicator());
    
    if (ParallelDescriptor::IOProcessor()) {
        std::vector<int> displs(nranks);
        total_burned = 0;
        for (int r = 0; r < nranks; r++) {
            displs[r] = total_burned;
            total_burned += counts[r];
        }
        
        // ✅ Check total size before allocating
        if (total_burned > 100000) {
            Print() << "  🚨 ERROR: Too many burned cells (" << total_burned 
                    << "), reinitialization may fail!\n";
            Print() << "     Consider cleaning m_local_burned or using simpler reinit\n";
        }
        
        all_burned_i.resize(total_burned);
        all_burned_j.resize(total_burned);
        
        MPI_Gatherv(local_burned_i.data(), local_count, MPI_INT,
                    all_burned_i.data(), counts.data(), displs.data(), MPI_INT,
                    0, ParallelDescriptor::Communicator());
        
        MPI_Gatherv(local_burned_j.data(), local_count, MPI_INT,
                    all_burned_j.data(), counts.data(), displs.data(), MPI_INT,
                    0, ParallelDescriptor::Communicator());
    } else {
        MPI_Gather(&local_count, 1, MPI_INT, nullptr, 0, MPI_INT,
                   0, ParallelDescriptor::Communicator());
        
        MPI_Gatherv(local_burned_i.data(), local_count, MPI_INT,
                    nullptr, nullptr, nullptr, MPI_INT,
                    0, ParallelDescriptor::Communicator());
        
        MPI_Gatherv(local_burned_j.data(), local_count, MPI_INT,
                    nullptr, nullptr, nullptr, MPI_INT,
                    0, ParallelDescriptor::Communicator());
    }
    
    // ========================================
    // STEP 3: Broadcast to ALL ranks
    // ========================================
    MPI_Bcast(&total_burned, 1, MPI_INT, 0, ParallelDescriptor::Communicator());
    
    if (!ParallelDescriptor::IOProcessor()) {
        all_burned_i.resize(total_burned);
        all_burned_j.resize(total_burned);
    }
    
    MPI_Bcast(all_burned_i.data(), total_burned, MPI_INT,
              0, ParallelDescriptor::Communicator());
    MPI_Bcast(all_burned_j.data(), total_burned, MPI_INT,
              0, ParallelDescriptor::Communicator());
#else
    total_burned = local_burned_i.size();
    all_burned_i = local_burned_i;
    all_burned_j = local_burned_j;
#endif
    
    if (total_burned == 0) {
        if (ParallelDescriptor::IOProcessor()) {
            Print() << "  ⚠️  No burned cells for reinitialization!\n";
        }
        return;
    }
    
    // ========================================
    // STEP 4: Each rank reinitializes ITS boxes
    // ========================================
    for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> phi_arr = m_fire_data.array(mfi);
        
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); k++) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++) {
                    
                    Real x = prob_lo[0] + (i + 0.5) * dx[0];
                    Real y = prob_lo[1] + (j + 0.5) * dx[1];
                    
                    Real phi_old = phi_arr(i, j, k);
                    Real sign_old = (phi_old >= 0.0) ? 1.0 : -1.0;
                    
                    Real min_dist = 1e10;
                    for (int n = 0; n < total_burned; n++) {
                        int bi = all_burned_i[n];
                        int bj = all_burned_j[n];
                        
                        Real bx = prob_lo[0] + (bi + 0.5) * dx[0];
                        Real by = prob_lo[1] + (bj + 0.5) * dx[1];
                        
                        Real dist = std::sqrt((x - bx)*(x - bx) + (y - by)*(y - by));
                        min_dist = std::min(min_dist, dist);
                    }
                    
                    Real phi_new_val = sign_old * min_dist;
                    
                    // Preserve progress
                    if (sign_old < 0 && phi_old < phi_new_val) {
                        phi_new_val = phi_old;
                    } else if (sign_old > 0 && phi_old < phi_new_val) {
                        phi_new_val = phi_old;
                    }
                    
                    phi_arr(i, j, k) = phi_new_val;
                }
            }
        }
    }
    
    m_fire_data.FillBoundary(m_geom->periodicity());
    
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "  ✅ Reinit complete (used " << total_burned << " recent burned cells)\n";
    }
}

// ============================================================================
// UTILITIES
// ============================================================================

void FireLevelSet::ComputeTerrainSlopes(const MultiFab* terrain, const Real* dx) {
    if (!terrain) return;
    
    Real dx_inv = 1.0 / (2.0 * dx[0]);
    Real dy_inv = 1.0 / (2.0 * dx[1]);
    
    for (MFIter mfi(*terrain); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> z = terrain->const_array(mfi);
        Array4<Real> sx = m_slope_x->array(mfi);
        Array4<Real> sy = m_slope_y->array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            if (i > bx.smallEnd(0) && i < bx.bigEnd(0)) {
                sx(i,j,k) = (z(i+1,j,k) - z(i-1,j,k)) * dx_inv;
            } else {
                sx(i,j,k) = 0.0;
            }
            
            if (j > bx.smallEnd(1) && j < bx.bigEnd(1)) {
                sy(i,j,k) = (z(i,j+1,k) - z(i,j-1,k)) * dy_inv;
            } else {
                sy(i,j,k) = 0.0;
            }
        });
    }
    
    m_slope_x->FillBoundary(m_geom->periodicity());
    m_slope_y->FillBoundary(m_geom->periodicity());
}

bool FireLevelSet::IsNonburnable(int i, int j) {
    return false;
}

Real FireLevelSet::ComputeAdaptiveTimestep(Real max_velocity) {
    if (max_velocity < 1e-10) {
        return m_simulation_dt_max;
    }
    
    const Real* dx = m_geom->CellSize();
    Real min_dx = std::min(dx[0], dx[1]);
    Real dt = m_target_cfl * min_dx / max_velocity;
    dt = std::min(dt, m_simulation_dt_max);
    
    return dt;
}

// ============================================================================
// OUTPUT
// ============================================================================

void FireLevelSet::GetPlotVar(const std::string& name, MultiFab& mf, int lev) const {
    mf.setVal(0.0);
    
    if (name == "fire_phi") {
        MultiFab::Copy(mf, m_fire_data, FireVar::phi, 0, 1, 0);
        return;
    }
    if (name == "fire_fuel_frac") {
        MultiFab::Copy(mf, m_fire_data, FireVar::fuel_frac, 0, 1, 0);
        return;
    }
    if (name == "fire_velocity") {
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            Array4<Real> vel = mf.array(mfi);
            
            for (const auto& node : m_local_tagged) {
                if (bx.contains(IntVect(node.i, node.j, 0))) {
                    Real vel_ms = node.velocity * 0.3048 / 60.0;
                    vel(node.i, node.j, 0) = vel_ms;
                }
            }
        }
        return;
    }

    if (name == "fire_tagged") {
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            Array4<Real> tagged_arr = mf.array(mfi);
            
            // Mark all tagged cells with value 1.0
            for (const auto& node : m_local_tagged) {
                if (bx.contains(IntVect(node.i, node.j, 0))) {
                    tagged_arr(node.i, node.j, 0) = 1.0;
                }
            }
        }
        return;
    }    
}

Real FireLevelSet::get_fire_area() const {
    const Real* dx = m_geom->CellSize();
    Real cell_area = dx[0] * dx[1];
    
    int local_burned = m_local_burned.size();
    int global_burned = local_burned;
#ifdef AMREX_USE_MPI
    ParallelDescriptor::ReduceIntSum(global_burned);
#endif
    
    return global_burned * cell_area;
}

void FireLevelSet::GetPlotVarNames(Vector<std::string>& names) const {
    names.push_back("fire_phi");
    names.push_back("fire_fuel_frac");
    names.push_back("fire_heat_flux");
    names.push_back("fire_burned_area");
    names.push_back("fire_fuel_time");
    names.push_back("fire_velocity");
    names.push_back("fire_flame_length");
    names.push_back("fire_tagged");
}

void FireLevelSet::GetHeatFlux(MultiFab& heat_flux) {
    MultiFab::Copy(heat_flux, m_fire_data, FireVar::heat_flux, 0, 1, 0);
}

void FireLevelSet::GetMoistureFlux(MultiFab& moisture_flux) {
    moisture_flux.setVal(0.0);
}