#include "ERF_FireLevelSet.H"

using namespace amrex;

// ============================================================================
// ROTHERMEL SURFACE FIRE SPREAD MODEL (ELMFIRE-style) - LOCAL VERSION
// ============================================================================

void FireLevelSet::SurfaceSpreadRateLocal() {
    BL_PROFILE("FireLevelSet::SurfaceSpreadRateLocal");
    
    // // ✅ DIAGNOSTIC: Entry
    // AllPrint() << "  [Rank " << ParallelDescriptor::MyProc() 
    //            << "] SurfaceSpreadRateLocal() START: " 
    //            << m_local_tagged.size() << " cells\n";
    
    // if (m_local_tagged.empty()) {
    //     AllPrint() << "  [Rank " << ParallelDescriptor::MyProc() 
    //                << "] Early exit: no cells\n";
    //     return;
    // }
    
    const Real BTUPFT2MIN_TO_KWPM2 = 1.055 / (60.0 * 0.3048 * 0.3048);
    
    int cells_processed = 0;
    int sample_count = 0;
    
    Real max_velocity = 0.0;
    Real min_velocity = 1e10;
    Real sum_velocity = 0.0;

    for (auto& node : m_local_tagged) {
        if (node.burned) continue;
        
        // Skip non-burnable cells (0 = custom non-burnable, 91/98/99 = Scott & Burgan non-burnable)
        if (node.ifbfm == 0 || node.ifbfm == 91 || node.ifbfm == 98 || node.ifbfm == 99) {
            node.ir = 0.0;
            node.phiw_surface = 0.0;
            node.phis_surface = 0.0;
            node.vs0 = 0.0;
            node.velocity_dms_surface = 0.0;
            node.hpua_surface = 0.0;
            node.flin_dms_surface = 0.0;
            node.ux = 0.0;
            node.uy = 0.0;
            node.velocity = 0.0;
            continue;
        }
        
        // ========================================
        // Get fuel model parameters
        // ========================================
        FuelModelTable fmt = GetFuelModelTable(node.ifbfm, node.mlh);
        
        // ========================================
        // Calculate fuel-specific wind/slope coefficients
        // ========================================
        Real sigma_bar = fmt.sigma_bar;
        Real beta = fmt.beta;
        Real beta_opt = fmt.beta_opt;
        
        // Wind coefficient (ELMFIRE/Rothermel 1972)
        Real C = 7.47 * std::exp(-0.133 * std::pow(sigma_bar, 0.55));
        Real E = 0.715 * std::exp(-3.59e-4 * sigma_bar);
        Real b_coeff = 0.02526 * std::pow(sigma_bar, 0.54);
        Real phiwterm = C * std::pow(beta / beta_opt, -E);

        // Slope coefficient (ELMFIRE/Rothermel 1972)
        Real phisterm = 5.275 * std::pow(beta, -0.3);
        
        // Store in node for use in UxAndUyElliptical()
        node.b_coeff = b_coeff;
        node.phiwterm = phiwterm;
        node.phisterm = phisterm;
        
        // ========================================
        // Calculate live fuel moisture of extinction
        // ========================================
        Real m[6] = {node.m1, node.m10, node.m100, node.m1, node.mlh, node.mlw};
        
        Real mprimenumer[4];
        Real sum_mprimenumer = 0.0;
        for (int i = 0; i < 4; i++) {
            mprimenumer[i] = fmt.wprimenumer[i] * m[i];
            sum_mprimenumer += mprimenumer[i];
        }
        
        Real mex_live = fmt.mex_live * (1.0 - fmt.r_mprimedenome14sum_mex_dead * sum_mprimenumer) - 0.226;
        mex_live = std::max(mex_live, fmt.mex_dead);
        
        // ========================================
        // Calculate fuel moisture content (FMC)
        // ========================================
        Real fmc[6];
        for (int i = 0; i < 6; i++) {
            fmc[i] = fmt.f[i] * m[i];
        }
        
        // ========================================
        // Heat of ignition (weighted by fuel fractions)
        // ========================================
        Real qig_weighted = 0.0;
        for (int i = 0; i < 6; i++) {
            Real qig_i = 250.0 + 1116.0 * m[i];
            qig_weighted += fmt.f[i] * qig_i;
        }

        Real epsilon = std::exp(-138.0 / sigma_bar);
        Real rhobepsqig = fmt.rhob * epsilon * qig_weighted;
        
        // ========================================
        // Calculate dead fuel reaction intensity
        // ========================================
        Real m_dead = fmc[0] + fmc[1] + fmc[2] + fmc[3];
        Real momex_dead = m_dead / fmt.mex_dead;
        Real momex2 = momex_dead * momex_dead;
        Real momex3 = momex2 * momex_dead;
        
        Real etam_dead = 1.0 - 2.59*momex_dead + 5.11*momex2 - 3.52*momex3;
        etam_dead = std::max(0.0, std::min(etam_dead, 1.0));
        
        Real ir_dead = fmt.gp_wnd_emd_es_hoc * etam_dead;
        
        // ========================================
        // Calculate live fuel reaction intensity
        // ========================================
        Real m_live = fmc[4] + fmc[5];
        Real momex_live = m_live / mex_live;
        momex2 = momex_live * momex_live;
        momex3 = momex2 * momex_live;
        
        Real etam_live = 1.0 - 2.59*momex_live + 5.11*momex2 - 3.52*momex3;
        etam_live = std::max(0.0, std::min(etam_live, 1.0));
        
        Real ir_live = fmt.gp_wnl_eml_es_hoc * etam_live;

        // Total reaction intensity
        Real ir_total = ir_dead + ir_live;

        // ========================================
        // Wind factor
        // ========================================
        Real ws_limit = 0.9 * ir_total;
        Real wsmf_limited = std::min(node.wsmf, ws_limit);
        node.phiw_surface = phiwterm * std::pow(wsmf_limited, b_coeff);

        // ========================================
        // Slope factor
        // ========================================
        Real phis_max = phiwterm * std::pow(ws_limit, b_coeff);
        node.phis_surface = std::min(phisterm * node.tanslp2, phis_max);

        // TEMPORARILY DISABLE SLOPE EFFECTS
        node.phis_surface = 0.0;  // ✅ FORCE ZERO SLOPE EFFECT!
        node.phisx = 0.0;          // ✅ FORCE ZERO SLOPE COMPONENTS!
        node.phisy = 0.0;
        
        // ========================================
        // Base spread rate
        // ========================================
        if (rhobepsqig > 1e-10) {
            node.vs0 = node.adj * ir_total * fmt.xi / rhobepsqig;
        } else {
            node.vs0 = 0.0;
        }
        
        // ========================================
        // DMS velocity
        // ========================================
        node.velocity_dms_surface = node.vs0 * (1.0 + node.phis_surface + node.phiw_surface);
        
        // ========================================
        // Convert to SI units
        // ========================================
        node.ir = ir_total * BTUPFT2MIN_TO_KWPM2;
        node.hpua_surface = node.ir * fmt.tr * 60.0;
        node.flin_dms_surface = fmt.tr * node.ir * node.velocity_dms_surface * 0.3048;
        
        // ========================================
        // Slope components for ellipse
        // ========================================
        Real slope_mag = std::sqrt(node.slope_x*node.slope_x + node.slope_y*node.slope_y);
        if (slope_mag > 1e-10) {
            Real phis = 5.275 * std::pow(slope_mag, 2) * std::exp(-6.884 * slope_mag);
            node.phisx = phis * node.slope_x / slope_mag;
            node.phisy = phis * node.slope_y / slope_mag;
        } else {
            node.phisx = 0.0;
            node.phisy = 0.0;
        }

        // TEMPORARILY DISABLE SLOPE EFFECTS
        node.phisx = 0.0;  // ✅ FORCE ZERO!
        node.phisy = 0.0;  // ✅ FORCE ZERO!        
        
        // ========================================
        // Terrain projection factors
        // ========================================
        Real cosnslp = 1.0 / std::sqrt(1.0 + node.tanslp2);
        node.uxousx = cosnslp;
        node.uyousy = cosnslp;
    }
}

// ============================================================================
// ELLIPTICAL FIRE PROPAGATION - LOCAL VERSION
// ============================================================================

void FireLevelSet::UxAndUyEllipticalLocal(int step, Real time) {
    BL_PROFILE("FireLevelSet::UxAndUyEllipticalLocal");

    const Real FTPMIN_TO_MPS = 0.3048 / 60.0;
    const Real KWPM2_TO_BTUPFT2MIN = 60.0 * 0.3048 * 0.3048 / 1.055;
    const Real MAX_LOW = 8.0;

    for (auto& node : m_local_tagged) {
        if (node.burned) continue;

        // Skip non-burnable cells (already handled in SurfaceSpreadRateLocal, but double-check)
        if (node.ifbfm == 0 || node.ifbfm == 91 || node.ifbfm == 98 || node.ifbfm == 99) {
            node.ux = 0.0;
            node.uy = 0.0;
            node.velocity = 0.0;
            continue;
        }

        // ========================================
        // STEP 1: Wind/slope vector in DMS coords
        // ========================================
        Real aphiw = 1.0 * node.phiw_surface;
        Real phiwx = aphiw * std::sin((node.wd20 - 180.0) * M_PI / 180.0);
        Real phiwy = aphiw * std::cos((node.wd20 - 180.0) * M_PI / 180.0);
        Real phix = node.phisx + phiwx;
        Real phiy = node.phisy + phiwy;
        
        Real phimag = std::max(std::sqrt(phix*phix + phiy*phiy), 1e-10);

        if (phimag > 1e-10) {
            node.normx_dms = phix / phimag;
            node.normy_dms = phiy / phimag;
        } else {
            node.normx_dms = 1.0;
            node.normy_dms = 0.0;
        }
        
        // ========================================
        // STEP 2: DMS velocity
        // ========================================
        node.velocity_dms = node.vs0 * (1.0 + node.phis_surface + node.phiw_surface);

        // ========================================
        // STEP 3: Effective wind speed and L/W ratio
        // ========================================
        Real b_coeff_inv = 1.0 / node.b_coeff;
        Real wsmfeff_coeff = std::pow(1.0 / node.phiwterm, b_coeff_inv);
        Real wsmfeff = wsmfeff_coeff * std::pow(phimag, b_coeff_inv);
        Real ir_limit = 0.9 * KWPM2_TO_BTUPFT2MIN * node.ir;
        wsmfeff = std::min(wsmfeff, ir_limit);
        
        node.low = 0.936 * std::exp(0.2566 * wsmfeff) 
                 + 0.461 * std::exp(-0.1548 * wsmfeff) - 0.397;
        node.low = std::min(node.low, MAX_LOW);
        node.low = std::max(node.low, 1.0);
        
        // ========================================
        // STEP 4: Back fire velocity
        // ========================================
        Real boh;
        if (node.low > 0.999 && node.low < 1.001) {
            boh = 1.0;
        } else {
            Real sqrt_low2_m1 = std::sqrt(node.low * node.low - 1.0);
            boh = (node.low - sqrt_low2_m1) / (node.low + sqrt_low2_m1);
        }
        node.vback = boh * node.velocity_dms;
        
        // ========================================
        // STEP 5: Ellipse semi-axes
        // ========================================
        Real a = 0.5 * (node.velocity_dms + node.vback);
        Real b = 0.5 * (node.velocity_dms + node.vback) / node.low;
        
        // ========================================
        // STEP 6: Angle between gradient and DMS direction
        // ========================================
        Real cosang = node.normx * node.normx_dms + node.normy * node.normy_dms;  
        Real sinang = node.normx * node.normy_dms - node.normy * node.normx_dms;
        
        // ========================================
        // STEP 7: Velocity in rotated frame
        // ========================================
        Real denom = std::max(std::sqrt(a*a*cosang*cosang + b*b*sinang*sinang), 1e-10);
        Real dydt = (a*a*cosang / denom) + 0.5 * (node.velocity_dms - node.vback);
        Real dxdt = (b*b*sinang / denom);
        
        // ========================================
        // STEP 8: Rotate back to geographic coords
        // ========================================
        Real dxdt_rotated = dydt * node.normx_dms + dxdt * node.normy_dms;
        Real dydt_rotated = dydt * node.normy_dms - dxdt * node.normx_dms;
        
        // ========================================
        // STEP 9: Apply terrain projection and convert to m/s
        // ========================================
        node.ux = dxdt_rotated * node.uxousx * FTPMIN_TO_MPS;  
        node.uy = dydt_rotated * node.uyousy * FTPMIN_TO_MPS;

        node.velocity = std::sqrt(dxdt_rotated*dxdt_rotated + dydt_rotated*dydt_rotated);
    }
}

// ============================================================================
// NORMAL VECTOR CALCULATION - LOCAL VERSION
// ============================================================================

void FireLevelSet::CalcNormalVectorsLocal(int step) {
    const Real* dx = m_geom->CellSize();
    Real dx_inv = 1.0 / (2.0 * dx[0]);
    Real dy_inv = 1.0 / (2.0 * dx[1]);

    const Real BIG = 3.0e4;  // Maximum allowed gradient magnitude
    const Real EPSILON = 1.0e-10;  // Minimum for normalization    
    
    for (auto& node : m_local_tagged) {
        if (node.burned) continue;
        
        // ✅ Read from MultiFab (includes ghost cells!)
        Real phi_val = 0.0;
        Real phi_ip1 = 0.0, phi_im1 = 0.0;
        Real phi_jp1 = 0.0, phi_jm1 = 0.0;
        
        bool found = false;
        for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
            const Box& grownbox = mfi.growntilebox();
            
            if (grownbox.contains(IntVect(node.i, node.j, 0))) {
                Array4<Real const> phi = m_fire_data.const_array(mfi);
                
                phi_val = phi(node.i, node.j, 0);
                phi_ip1 = phi(node.i+1, node.j, 0);
                phi_im1 = phi(node.i-1, node.j, 0);
                phi_jp1 = phi(node.i, node.j+1, 0);
                phi_jm1 = phi(node.i, node.j-1, 0);
                
                found = true;
                break;
            }
        }
        
        if (!found) continue;  // Cell not in this rank's boxes
        
        // Central difference for gradient
        Real dphidx = (phi_ip1 - phi_im1) * dx_inv;
        Real dphidy = (phi_jp1 - phi_jm1) * dy_inv;
        
        // CLIP gradients to prevent instability (like ELMFIRE)
        dphidx = std::max(std::min(dphidx, BIG), -BIG);
        dphidy = std::max(std::min(dphidy, BIG), -BIG);
        
        Real gradmag = std::sqrt(dphidx*dphidx + dphidy*dphidy);
        
        if (gradmag > 1e-10) {
            node.normx = dphidx / gradmag;
            node.normy = dphidy / gradmag;
        } else {
            node.normx = 0.0;
            node.normy = 1.0;
        }
        
        node.dphidx = dphidx;
        node.dphidy = dphidy;
        node.gradmag = gradmag;
        
        if (step == 1) {
            node.phip_old = phi_val;
        }
    }
}

// ============================================================================
// FLUX LIMITER - LOCAL VERSION
// ============================================================================

void FireLevelSet::ApplyFluxLimiterLocal() {
    for (auto& node : m_local_tagged) {
        if (node.burned) continue;
        
        Real gradmag = node.gradmag;
        
        if (gradmag < 0.5) {
            node.dphidx_limited = node.dphidx * 0.5;
            node.dphidy_limited = node.dphidy * 0.5;
        } else if (gradmag > 1.5) {
            node.dphidx_limited = node.dphidx * (1.5 / gradmag);
            node.dphidy_limited = node.dphidy * (1.5 / gradmag);
        } else {
            node.dphidx_limited = node.dphidx;
            node.dphidy_limited = node.dphidy;
        }
    }
}

// ============================================================================
// RK2 INTEGRATION - LOCAL VERSION
// ============================================================================

void FireLevelSet::Rk2IntegrateLocal(const Real& dt, int step) {
    BL_PROFILE("FireLevelSet::Rk2IntegrateLocal");
    
    // ✅ DIAGNOSTIC
    static int call_count = 0;
    bool print_detail = (call_count < 3);
    
    if (print_detail && !m_local_tagged.empty()) {
        const auto& sample = m_local_tagged[0];
        AllPrint() << "  [Rank " << ParallelDescriptor::MyProc() 
                   << "] Rk2IntegrateLocal (step " << step << ", dt=" << dt << "):\n";
        AllPrint() << "     Sample cell (" << sample.i << "," << sample.j << "):\n";
        AllPrint() << "       phip_old (before): " << sample.phip_old << "\n";
        AllPrint() << "       ux=" << sample.ux << ", uy=" << sample.uy << " m/s\n";
        AllPrint() << "       dphidx_limited=" << sample.dphidx_limited 
                   << ", dphidy_limited=" << sample.dphidy_limited << "\n";
    }
    
    // ========================================
    // STEP 1: Update phip_old for all nodes
    // ========================================
    for (auto& node : m_local_tagged) {
        if (node.burned) continue;
        
        Real dphi = -(node.ux * node.dphidx_limited + node.uy * node.dphidy_limited) * dt;
        
        if (step == 1) {
            // First RK2 stage: store intermediate value
            node.phip_rk1 = node.phip_old + dphi;
        } else {
            // Second RK2 stage: compute final value
            node.phip_old = 0.5 * (node.phip_old + node.phip_rk1 + dphi);
        }
    }
    
    // ✅ DIAGNOSTIC: Show what changed
    if (print_detail && !m_local_tagged.empty()) {
        const auto& sample = m_local_tagged[0];
        Real dphi = -(sample.ux * sample.dphidx_limited + sample.uy * sample.dphidy_limited) * dt;
        AllPrint() << "       dphi: " << dphi << "\n";
        if (step == 1) {
            AllPrint() << "       phip_rk1 (after): " << sample.phip_rk1 << "\n";
        } else {
            AllPrint() << "       phip_old (after): " << sample.phip_old << "\n";
        }
    }
    
    // ========================================
    // STEP 2: Write back to MultiFab (only after step 2)
    // ========================================
    if (step == 2) {
        int total_written = 0;
        
        // ✅ OPTIMIZED: Loop over boxes first (not nodes!)
        for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
            const Box& validbox = mfi.validbox();
            Array4<Real> phi = m_fire_data.array(mfi);
            
            int nodes_written_this_box = 0;
            
            // For each box, check all nodes
            for (const auto& node : m_local_tagged) {
                if (validbox.contains(IntVect(node.i, node.j, 0))) {
                    phi(node.i, node.j, 0) = node.phip_old;
                    nodes_written_this_box++;
                }
            }
            
            total_written += nodes_written_this_box;
            
            if (print_detail && nodes_written_this_box > 0) {
                AllPrint() << "       Box " << mfi.index() << ": wrote " 
                           << nodes_written_this_box << " nodes\n";
            }
        }
        
        // ✅ DIAGNOSTIC: Verify write
        if (print_detail && !m_local_tagged.empty()) {
            const auto& sample = m_local_tagged[0];
            
            for (MFIter mfi(m_fire_data); mfi.isValid(); ++mfi) {
                const Box& validbox = mfi.validbox();
                if (validbox.contains(IntVect(sample.i, sample.j, 0))) {
                    Array4<Real const> phi = m_fire_data.const_array(mfi);
                    Real phi_after = phi(sample.i, sample.j, 0);
                    
                    AllPrint() << "     VERIFICATION:\n";
                    AllPrint() << "       Total nodes written: " << total_written << "\n";
                    AllPrint() << "       phi(" << sample.i << "," << sample.j 
                               << ") = " << phi_after << "\n";
                    AllPrint() << "       Expected: " << sample.phip_old << "\n";
                    AllPrint() << "       Match: " << (std::abs(phi_after - sample.phip_old) < 1e-6 ? "✅" : "❌") << "\n";
                    
                    call_count++;
                    break;
                }
            }
        }
        
        // ✅ Sync ghost cells after writing
        m_fire_data.FillBoundary(m_geom->periodicity());
    }
}