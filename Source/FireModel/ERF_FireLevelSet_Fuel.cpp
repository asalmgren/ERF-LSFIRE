#include "ERF_FireLevelSet.H"
#include <AMReX_ParmParse.H>
#include <fstream>
#include <sstream>
#include <map>

using namespace amrex;

// ============================================================================
// Map MODIS-IGBP LU_INDEX to Scott & Burgan 40 fuel model
// ============================================================================
int FireLevelSet::MapLUINDEXtoFuelModel(int lu_index) {
    switch(lu_index) {
        case 1:  return 165; // Evergreen Needleleaf Forest → TU5
        case 2:  return 145; // Evergreen Broadleaf Forest → TL5
        case 3:  return 161; // Deciduous Needleleaf Forest → TU1
        case 4:  return 161; // Deciduous Broadleaf Forest → TU1
        case 5:  return 162; // Mixed Forest → TU2
        case 6:  return 124; // Closed Shrublands → GS4
        case 7:  return 122; // Open Shrublands → GS2
        case 8:  return 123; // Woody Savannas → GS3
        case 9:  return 123; // Savannas → GS3
        case 10: return 102; // Grasslands → GR2
        case 11: return 0;   // Permanent Wetlands → Non-burnable
        case 12: return 101; // Croplands → GR1
        case 13: return 162; // Urban and Built-Up → TU2
        case 14: return 121; // Cropland/Natural Mosaic → GS1
        case 15: return 0;   // Snow and Ice → Non-burnable
        case 16: return 0;   // Barren or Sparsely Vegetated → Non-burnable
        case 17: return 0;   // Water Bodies → Non-burnable
        case 18: return 102; // Wooded Tundra → GR2
        case 19: return 101; // Mixed Tundra → GR1
        case 20: return 0;   // Barren Tundra → Non-burnable
        case 21: return 0;   // Lake → Non-burnable
        default:
            if (ParallelDescriptor::IOProcessor()) {
                Print() << "WARNING: Unknown MODIS LU_INDEX " << lu_index 
                        << ", using default fuel model " << m_default_fuel_model << "\n";
            }
            return m_default_fuel_model;
    }
}

// ============================================================================
// Remove non-burnable cells and update fuel models from LUINDEX (LOCAL VERSION)
// ============================================================================
void FireLevelSet::UpdateFuelModelsFromLUINDEX(const iMultiFab* luindex) {
    if (!luindex) return;
    
    // ✅ Each rank updates its LOCAL tagged cells
    int local_cells_updated = 0;
    int local_cells_removed = 0;
    
    for (auto it = m_local_tagged.begin(); it != m_local_tagged.end(); ) {
        if (it->burned) { 
            ++it; 
            continue; 
        }
        
        bool found = false;
        int fuel_model = m_default_fuel_model;
        
        // ✅ Search in this rank's boxes
        for (MFIter mfi(*luindex); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (!bx.contains(IntVect(it->i, it->j, 0))) continue;
            
            Array4<const int> lu_arr = luindex->const_array(mfi);
            int lu_index = lu_arr(it->i, it->j, 0);
            fuel_model = MapLUINDEXtoFuelModel(lu_index);
            
            found = true;
            break;
        }
        
        // Remove non-burnable cells (fuel_model == 0)
        if (fuel_model == 0) {
            it = m_local_tagged.erase(it);
            local_cells_removed++;
            continue;
        }
        
        // Update fuel model if changed
        if (found && fuel_model != it->ifbfm) {
            it->ifbfm = fuel_model;
            local_cells_updated++;
        }
        
        ++it;
    }
    
    // ✅ Reduce to get global counts
    int global_cells_updated = local_cells_updated;
    int global_cells_removed = local_cells_removed;
    
#ifdef AMREX_USE_MPI
    ParallelDescriptor::ReduceIntSum(global_cells_updated);
    ParallelDescriptor::ReduceIntSum(global_cells_removed);
#endif
    
    if (ParallelDescriptor::IOProcessor() && (global_cells_updated > 0 || global_cells_removed > 0)) {
        Print() << "  Updated " << global_cells_updated << " cells with LU_INDEX fuel models\n";
        Print() << "  Removed " << global_cells_removed << " non-burnable cells (water/ice/barren)\n";
    }
}

// ============================================================================
// Load fuel model database from CSV
// ============================================================================
void FireLevelSet::LoadFuelDatabase() {
    // ✅ Only IO processor loads the file, but ALL ranks get the data
    if (m_fuel_database_loaded) return;
    
    std::string csv_path = "fuel_models.csv";
    ParmParse pp("prob");
    pp.query("fire_fuel_database", csv_path);
    
    // ✅ ALL ranks need the fuel database (not just Rank 0!)
    std::ifstream infile(csv_path);
    if (!infile.is_open()) {
        if (ParallelDescriptor::IOProcessor()) {
            Print() << "WARNING: Could not open fuel database: " << csv_path << "\n";
            Print() << "Using default fuel model parameters\n";
        }
        m_fuel_database_loaded = true;
        return;
    }
    
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "Loading fuel database from: " << csv_path << "\n";
    }
    
    std::string line;
    int count = 0;
    
    // Skip header if present
    std::getline(infile, line);
    
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;
        
        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }
        
        if (tokens.size() < 14) continue;
        
        try {
            FuelModelData fm;
            fm.id = std::stoi(tokens[0]);
            fm.name = tokens[1];
            fm.has_live_fuels = (tokens[2] == ".TRUE." || tokens[2] == "TRUE" || tokens[2] == "true");
            fm.w_1hr = std::stod(tokens[3]);
            fm.w_10hr = std::stod(tokens[4]);
            fm.w_100hr = std::stod(tokens[5]);
            fm.w_lh = std::stod(tokens[6]);
            fm.w_lw = std::stod(tokens[7]);
            fm.sigma_1hr = std::stod(tokens[8]);
            fm.sigma_10hr = std::stod(tokens[9]);
            fm.sigma_100hr = std::stod(tokens[10]);
            fm.depth = std::stod(tokens[11]);
            fm.mex_dead = std::stod(tokens[12]) / 100.0;
            fm.hoc = std::stod(tokens[13]);
            
            m_fuel_database[fm.id] = fm;
            count++;
        } catch (const std::exception& e) {
            if (ParallelDescriptor::IOProcessor()) {
                Print() << "WARNING: Error parsing fuel model line: " << line << "\n";
                Print() << "  Exception: " << e.what() << "\n";
            }
        }
    }
    
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "Loaded " << count << " fuel models\n";
    }
    
    m_fuel_database_loaded = true;
    infile.close();
}

// ============================================================================
// Get fuel model table (raw fuel properties, ELMFIRE-style)
// ============================================================================
FuelModelTable FireLevelSet::GetFuelModelTable(int fuel_model_id, Real mlh) {
    FuelModelTable fmt;
    FuelModelData fm;
    
    // ✅ ALL ranks need access to fuel database
    if (m_fuel_database.find(fuel_model_id) != m_fuel_database.end()) {
        fm = m_fuel_database[fuel_model_id];
    } else {
        // Only print warning once
        static std::set<int> warned_models;
        if (warned_models.find(fuel_model_id) == warned_models.end()) {
            if (ParallelDescriptor::IOProcessor()) {
                Print() << "WARNING: Fuel model " << fuel_model_id << " not found in database\n";
                Print() << "Using hardcoded fallback\n";
            }
            warned_models.insert(fuel_model_id);
        }
        return GetHardcodedFuelModel(fuel_model_id, mlh);
    }

    // Convert tons/acre to lb/ft² (ELMFIRE expects lb/ft²)
    const Real TONPAC_TO_LBPFT2 = 1.0; // If needed, set to 20.0 for strict ELMFIRE
    Real w_1hr = fm.w_1hr * TONPAC_TO_LBPFT2;
    Real w_10hr = fm.w_10hr * TONPAC_TO_LBPFT2;
    Real w_100hr = fm.w_100hr * TONPAC_TO_LBPFT2;
    Real w_lh = fm.w_lh * TONPAC_TO_LBPFT2;
    Real w_lw = fm.w_lw * TONPAC_TO_LBPFT2;

    // Dynamic live fuel partitioning (ELMFIRE-style)
    Real livefrac = 0.0, deadfrac = 1.0;
    if (fm.has_live_fuels) {
        livefrac = std::min(std::max((mlh - 30.0) / (120.0 - 30.0), 0.0), 1.0);
        deadfrac = 1.0 - livefrac;
        w_1hr += deadfrac * w_lh;
        w_lh = livefrac * w_lh;
    }

    // Fuel loadings
    Real w_dead = w_1hr + w_10hr + w_100hr;
    Real w_live = w_lh + w_lw;
    Real w_total = w_dead + w_live;
    if (w_total < 1e-10) return GetHardcodedFuelModel(fuel_model_id, mlh);

    // Bulk density
    Real delta = fm.depth;
    fmt.rhob = w_total / delta;

    // Surface area to volume ratios
    Real sigma[6];
    sigma[0] = (w_1hr > 1e-10) ? fm.sigma_1hr : 0.0;
    sigma[1] = (w_10hr > 1e-10) ? fm.sigma_10hr : 0.0;
    sigma[2] = (w_100hr > 1e-10) ? fm.sigma_100hr : 0.0;
    sigma[3] = 0.0;
    sigma[4] = (w_lh > 1e-10) ? 1800.0 : 0.0;   // ELMFIRE default
    sigma[5] = (w_lw > 1e-10) ? 1600.0 : 0.0;   // ELMFIRE default

    // Fuel fractions (dead and live)
    Real a_dead = 0.0, a_live = 0.0;
    for (int i = 0; i < 4; i++) a_dead += sigma[i] * ((i==0)?w_1hr:(i==1)?w_10hr:(i==2)?w_100hr:0.0) / 32.0;
    a_live += sigma[4] * w_lh / 32.0;
    a_live += sigma[5] * w_lw / 32.0;
    Real a_total = std::max(a_dead + a_live, 1e-9);

    for (int i = 0; i < 4; i++) {
        Real w0 = (i==0)?w_1hr:(i==1)?w_10hr:(i==2)?w_100hr:0.0;
        fmt.f[i] = (a_dead > 1e-9 && sigma[i] > 0.0) ? (sigma[i]*w0/32.0)/a_dead : 0.0;
    }
    for (int i = 4; i < 6; i++) {
        Real w0 = (i==4)?w_lh:w_lw;
        fmt.f[i] = (a_live > 1e-9 && sigma[i] > 0.0) ? (sigma[i]*w0/32.0)/a_live : 0.0;
    }
    fmt.f_dead = a_dead / a_total;
    fmt.f_live = a_live / a_total;

    // feps (f × exp(-138/sigma))
    for (int i = 0; i < 6; i++) {
        fmt.feps[i] = (sigma[i] > 1e-10) ? fmt.f[i] * std::exp(-138.0 / sigma[i]) : 0.0;
    }

    // ========================================
    // Calculate weighted average SAV (ELMFIRE formula)
    // ========================================
    Real fsig[6];
    for (int i = 0; i < 6; i++) {
        fsig[i] = fmt.f[i] * sigma[i];
    }

    Real sig_dead = fsig[0] + fsig[1] + fsig[2] + fsig[3];
    Real sig_live = fsig[4] + fsig[5];

    Real sigma_bar = fmt.f_dead * sig_dead + fmt.f_live * sig_live;
    if (sigma_bar < 1.0) sigma_bar = 1.0;
    fmt.sigma_bar = sigma_bar;

    // Packing ratio
    Real rho_p = 32.0;
    Real beta = fmt.rhob / rho_p;
    Real beta_opt = 3.348 * std::pow(sigma_bar, -0.8189);
    fmt.beta = beta;
    fmt.beta_opt = beta_opt;

    // Propagating flux ratio (ELMFIRE/Rothermel)
    fmt.xi = std::exp((0.792 + 0.681 * std::sqrt(sigma_bar)) * (beta + 0.1)) / (192.0 + 0.2595 * sigma_bar);
    fmt.xi = std::max(0.01, std::min(fmt.xi, 0.99));

    // Residence time
    fmt.tr = 384.0 / sigma_bar;

    // Mineral damping coefficient (ELMFIRE)
    Real S_e = 0.01;
    fmt.eta_s = 0.174 / std::pow(S_e, 0.19);

    // Reaction intensity coefficients (ELMFIRE)
    Real hoc = fm.hoc;
    Real gamma_max = std::pow(sigma_bar, 1.5) / (495.0 + 0.0594 * std::pow(sigma_bar, 1.5));
    Real A = 133.0 * std::pow(sigma_bar, -0.7913);
    Real gamma = gamma_max * std::pow(beta / beta_opt, A) * std::exp(A * (1.0 - beta / beta_opt));
    fmt.gp_wnd_emd_es_hoc = gamma * w_dead * fmt.eta_s * hoc / delta;
    fmt.gp_wnl_eml_es_hoc = gamma * w_live * fmt.eta_s * hoc / delta;

    // Moisture of extinction
    fmt.mex_dead = fm.mex_dead;
    fmt.mex_live = 1.5; // Will be recalculated in physics

    // wprimenumer and fmex for MEX calculation
    for (int i = 0; i < 4; i++) {
        Real w0 = (i==0)?w_1hr:(i==1)?w_10hr:(i==2)?w_100hr:0.0;
        fmt.wprimenumer[i] = (sigma[i] > 1e-10) ? w0 * std::exp(-138.0 / sigma[i]) : 0.0;
        fmt.fmex[i] = fmt.f[i] * fmt.mex_dead;
    }
    Real mprimedenom_sum = 0.0;
    for (int i = 0; i < 4; i++) {
        Real w0 = (i==0)?w_1hr:(i==1)?w_10hr:(i==2)?w_100hr:0.0;
        if (sigma[i] > 1e-10) mprimedenom_sum += w0 * std::exp(-138.0 / sigma[i]);
    }
    fmt.r_mprimedenome14sum_mex_dead = (mprimedenom_sum > 1e-10) ? 1.0 / (mprimedenom_sum * fm.mex_dead) : 0.0;

    return fmt;
}

// ============================================================================
// Fallback: Hardcoded fuel model (GR1, short grass)
// ============================================================================
FuelModelTable FireLevelSet::GetHardcodedFuelModel(int fuel_model_id, Real mlh) {
    FuelModelTable fmt;
    
    // GR1 (short grass) fallback
    fmt.wprimenumer[0] = 0.00459; 
    fmt.wprimenumer[1] = 0.0; 
    fmt.wprimenumer[2] = 0.0; 
    fmt.wprimenumer[3] = 0.0;
    
    fmt.f[0] = 0.25; 
    fmt.f[1] = 0.0; 
    fmt.f[2] = 0.0; 
    fmt.f[3] = 0.0; 
    fmt.f[4] = 0.75; 
    fmt.f[5] = 0.0;
    
    fmt.feps[0] = 2200; 
    fmt.feps[1] = 0; 
    fmt.feps[2] = 0; 
    fmt.feps[3] = 0; 
    fmt.feps[4] = 2000; 
    fmt.feps[5] = 0;
    
    fmt.mex_dead = 0.15;
    fmt.mex_live = 1.5;
    fmt.rhob = 0.04;
    fmt.xi = 0.4;
    fmt.tr = 0.2;
    fmt.gp_wnd_emd_es_hoc = 50.0;
    fmt.gp_wnl_eml_es_hoc = 150.0;
    fmt.f_dead = 0.25;
    fmt.f_live = 0.75;
    
    for (int i = 0; i < 4; i++) {
        fmt.fmex[i] = fmt.f[i] * fmt.mex_dead;
    }
    
    Real sum_wprime = 0.0;
    for (int i = 0; i < 4; i++) {
        sum_wprime += fmt.wprimenumer[i];
    }
    fmt.r_mprimedenome14sum_mex_dead = (sum_wprime > 1e-10) ? 1.0 / (sum_wprime * fmt.mex_dead) : 0.0;
    
    fmt.sigma_bar = 2100.0;
    fmt.beta = 0.001;
    fmt.beta_opt = 0.02;
    fmt.eta_s = 0.41;
    
    return fmt;
}