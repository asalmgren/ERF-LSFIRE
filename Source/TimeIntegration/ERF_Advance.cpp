#include <ERF.H>
#include <ERF_Utils.H>
#include <ERF_Constants.H>

// nvcc forbids __device__ lambdas inside private member functions (ERF::Advance is private).
// These static helpers host the ParallelFor kernels so the private method is not the
// enclosing function of any __device__ lambda.

static void LaunchExtractSurfaceWindKernel(
    const amrex::Box& bx,
    const amrex::Array4<const amrex::Real>& u_3d,
    const amrex::Array4<const amrex::Real>& v_3d,
    const amrex::Array4<amrex::Real>& u_2d,
    const amrex::Array4<amrex::Real>& v_2d)
{
    amrex::ParallelFor(bx, [u_3d, v_3d, u_2d, v_2d] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
        u_2d(i,j,0) = 0.5 * (u_3d(i,j,1) + u_3d(i+1,j,1));
        v_2d(i,j,0) = 0.5 * (v_3d(i,j,1) + v_3d(i,j+1,1));
    });
}

static void LaunchApplyHeatFluxKernel(
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& cons_arr,
    const amrex::Array4<const amrex::Real>& heat_arr,
    amrex::Real dt_lev, amrex::Real cp_d_val, amrex::Real dz, int theta_comp)
{
    amrex::ParallelFor(bx,
        [cons_arr, heat_arr, dt_lev, cp_d_val, dz, theta_comp] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (k == 0) {
                amrex::Real heat_source = heat_arr(i,j,0) * dt_lev / (cp_d_val * dz);
                cons_arr(i,j,k,theta_comp) += heat_source;
            }
        });
}

static void LaunchApplyMoistureFluxKernel(
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& cons_arr,
    const amrex::Array4<const amrex::Real>& moist_arr,
    amrex::Real dt_lev, amrex::Real dz, int q1_comp)
{
    amrex::ParallelFor(bx,
        [cons_arr, moist_arr, dt_lev, dz, q1_comp] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (k == 0) {
                amrex::Real moisture_source = moist_arr(i,j,0) * dt_lev / dz;
                cons_arr(i,j,k,q1_comp) += moisture_source;
            }
        });
}

#ifdef ERF_USE_WINDFARM
#include <ERF_WindFarm.H>
#endif

using namespace amrex;

/**
 * Function that advances the solution at one level for a single time step --
 * this does some preliminaries then calls erf_advance
 *
 * @param[in] lev level of refinement (coarsest level is 0)
 * @param[in] time start time for time advance
 * @param[in] dt_lev time step for this time advance
 */

void
ERF::Advance (int lev, Real time, Real dt_lev, int iteration, int /*ncycle*/)
{
    BL_PROFILE("ERF::Advance()");

    // We must swap the pointers so the previous step's "new" is now this step's "old"
    std::swap(vars_old[lev], vars_new[lev]);

    MultiFab& S_old = vars_old[lev][Vars::cons];
    MultiFab& S_new = vars_new[lev][Vars::cons];

    MultiFab& U_old = vars_old[lev][Vars::xvel];
    MultiFab& V_old = vars_old[lev][Vars::yvel];
    MultiFab& W_old = vars_old[lev][Vars::zvel];

    MultiFab& U_new = vars_new[lev][Vars::xvel];
    MultiFab& V_new = vars_new[lev][Vars::yvel];
    MultiFab& W_new = vars_new[lev][Vars::zvel];

    // We need to set these because otherwise in the first call to erf_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    U_new.setVal(1.e34,U_new.nGrowVect());
    V_new.setVal(1.e34,V_new.nGrowVect());
    W_new.setVal(1.e34,W_new.nGrowVect());

    //
    // NOTE: the momenta here are not fillpatched (they are only used as scratch space)
    // If lev == 0 we have already FillPatched this in ERF::TimeStep
    //
    if (lev > 0) {
        FillPatchFineLevel(lev, time, {&S_old, &U_old, &V_old, &W_old},
                           {&S_old, &rU_old[lev], &rV_old[lev], &rW_old[lev]},
                           base_state[lev], base_state[lev]);
    }

    //
    // So we must convert the fillpatched to momenta, including the ghost values
    //
    const MultiFab* c_vfrac = nullptr;
    if (solverChoice.terrain_type == TerrainType::EB) {
        c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
    }

    VelocityToMomentum(U_old, rU_old[lev].nGrowVect(),
                       V_old, rV_old[lev].nGrowVect(),
                       W_old, rW_old[lev].nGrowVect(),
                       S_old, rU_old[lev], rV_old[lev], rW_old[lev],
                       Geom(lev).Domain(),
                       domain_bcs_type, c_vfrac);

    // Update the inflow perturbation update time and amplitude
    if (solverChoice.pert_type == PerturbationType::Source ||
        solverChoice.pert_type == PerturbationType::Direct ||
        solverChoice.pert_type == PerturbationType::CPM)
    {
        turbPert.calc_tpi_update(lev, dt_lev, U_old, V_old, S_old);
    }

    // If PerturbationType::Direct or CPM is selected, directly add the computed perturbation
    // on the conserved field
    if (solverChoice.pert_type == PerturbationType::Direct ||
        solverChoice.pert_type == PerturbationType::CPM)
    {
        auto m_ixtype = S_old.boxArray().ixType(); // Conserved term
        for (MFIter mfi(S_old,TileNoZ()); mfi.isValid(); ++mfi) {
            Box bx  = mfi.tilebox();
            const Array4<Real> &cell_data  = S_old.array(mfi);
            const Array4<const Real> &pert_cell = turbPert.pb_cell[lev].array(mfi);
            turbPert.apply_tpi(lev, bx, RhoTheta_comp, m_ixtype, cell_data, pert_cell);
        }
    }

    // configure SurfaceLayer params if needed
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::surface_layer) {
        if (m_SurfaceLayer) {
            IntVect ng = Theta_prim[lev]->nGrowVect();
            MultiFab::Copy(  *Theta_prim[lev], S_old, RhoTheta_comp, 0, 1, ng);
            MultiFab::Divide(*Theta_prim[lev], S_old, Rho_comp     , 0, 1, ng);
            if (solverChoice.moisture_type != MoistureType::None) {
                ng = Qv_prim[lev]->nGrowVect();

                MultiFab::Copy(  *Qv_prim[lev], S_old, RhoQ1_comp, 0, 1, ng);
                MultiFab::Divide(*Qv_prim[lev], S_old, Rho_comp  , 0, 1, ng);

                if (solverChoice.moisture_indices.qr > -1) {
                    MultiFab::Copy(  *Qr_prim[lev], S_old, solverChoice.moisture_indices.qr, 0, 1, ng);
                    MultiFab::Divide(*Qr_prim[lev], S_old, Rho_comp  , 0, 1, ng);
                } else {
                    Qr_prim[lev]->setVal(0.0);
                }
            }
            // NOTE: std::swap above causes the field ptrs to be out of date.
            //       Reassign the field ptrs for MAC avg computation.
            m_SurfaceLayer->update_mac_ptrs(lev, vars_old, Theta_prim, Qv_prim, Qr_prim);
            m_SurfaceLayer->update_pblh(lev, vars_old, z_phys_cc[lev].get(),
                                        solverChoice.moisture_indices);
            m_SurfaceLayer->update_fluxes(lev, time, S_old, z_phys_nd[lev], walldist[lev]);
        }
    }

#if defined(ERF_USE_WINDFARM)
    // **************************************************************************************
    // Update the windfarm sources
    // **************************************************************************************
    if (solverChoice.windfarm_type != WindFarmType::None) {
        advance_windfarm(Geom(lev), dt_lev, S_old,
                         U_old, V_old, W_old, vars_windfarm[lev],
                         Nturb[lev], SMark[lev], time);
    }

#endif

    // **************************************************************************************
    // Update the radiation sources with the "old" state
    // **************************************************************************************
    advance_radiation(lev, S_old, dt_lev);

#ifdef ERF_USE_SHOC
    // **************************************************************************************
    // Update the "old" state using SHOC
    // **************************************************************************************
    if (solverChoice.use_shoc) {
        // Get SFC fluxes from SurfaceLayer
        if (m_SurfaceLayer) {
            Vector<const MultiFab*> mfs = {&S_old, &U_old, &V_old, &W_old};
            m_SurfaceLayer->impose_SurfaceLayer_bcs(lev, mfs, Tau[lev],
                                                    SFS_hfx1_lev[lev].get() , SFS_hfx2_lev[lev].get() , SFS_hfx3_lev[lev].get(),
                                                    SFS_q1fx1_lev[lev].get(), SFS_q1fx2_lev[lev].get(), SFS_q1fx3_lev[lev].get(),
                                                    z_phys_nd[lev].get());
        }

        // Get Shoc tendencies and update the state
        Real* w_sub = (solverChoice.custom_w_subsidence) ? d_w_subsid[lev].data() : nullptr;
        compute_shoc_tendencies(lev, &S_old, &U_old, &V_old, &W_old, w_sub,
                                Tau[lev][TauType::tau13].get(), Tau[lev][TauType::tau23].get(),
                                SFS_hfx3_lev[lev].get()       , SFS_q1fx3_lev[lev].get()      ,
                                eddyDiffs_lev[lev].get()      , z_phys_nd[lev].get()          ,
                                dt_lev);
    }
#endif

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    int nvars = S_old.nComp();

    // Source array for conserved cell-centered quantities -- this will be filled
    //     in the call to make_sources in ERF_TI_slow_rhs_pre.H
    MultiFab cc_source(ba,dm,nvars,1); cc_source.setVal(0.0);

    // Source arrays for momenta -- these will be filled
    //     in the call to make_mom_sources in ERF_TI_slow_rhs_pre.H
    BoxArray ba_x(ba); ba_x.surroundingNodes(0);
    MultiFab xmom_source(ba_x,dm,1,1); xmom_source.setVal(0.0);

    BoxArray ba_y(ba); ba_y.surroundingNodes(1);
    MultiFab ymom_source(ba_y,dm,1,1); ymom_source.setVal(0.0);

    BoxArray ba_z(ba); ba_z.surroundingNodes(2);
    MultiFab zmom_source(ba_z,dm,1,1); zmom_source.setVal(0.0);
    MultiFab    buoyancy(ba_z,dm,1,1); buoyancy.setVal(0.0);

    amrex::Vector<MultiFab> state_old;
    amrex::Vector<MultiFab> state_new;

    // **************************************************************************************
    // Here we define state_old and state_new which are to be advanced
    // **************************************************************************************
    // Initial solution
    // Note that "old" and "new" here are relative to each RK stage.
    state_old.push_back(MultiFab(S_old      , amrex::make_alias, 0, nvars)); // cons
    state_old.push_back(MultiFab(rU_old[lev], amrex::make_alias, 0,     1)); // xmom
    state_old.push_back(MultiFab(rV_old[lev], amrex::make_alias, 0,     1)); // ymom
    state_old.push_back(MultiFab(rW_old[lev], amrex::make_alias, 0,     1)); // zmom

    // Final solution
    // state_new at the end of the last RK stage holds the t^{n+1} data
    state_new.push_back(MultiFab(S_new      , amrex::make_alias, 0, nvars)); // cons
    state_new.push_back(MultiFab(rU_new[lev], amrex::make_alias, 0,     1)); // xmom
    state_new.push_back(MultiFab(rV_new[lev], amrex::make_alias, 0,     1)); // ymom
    state_new.push_back(MultiFab(rW_new[lev], amrex::make_alias, 0,     1)); // zmom

    // **************************************************************************************
    // Tests on the reasonableness of the solution before the dycore
    // **************************************************************************************
    // Test for NaNs after dycore
    if (check_for_nans > 1) {
        if (verbose > 1) {
            amrex::Print() << "Testing old state and vels for NaNs before dycore" << std::endl;
        }
        check_state_for_nans(S_old);
        check_vels_for_nans(rU_old[lev],rV_old[lev],rW_old[lev]);
    }

    // We only test on low temp if we have a moisture model because we are protecting against
    //    the test on low temp inside the moisture models
    if (solverChoice.moisture_type != MoistureType::None) {
        if (verbose > 1) {
            amrex::Print() << "Testing on low temperature before dycore" << std::endl;
        }
        check_for_low_temp(S_old);
    } else {
        if (verbose > 1) {
            amrex::Print() << "Testing on negative temperature before dycore" << std::endl;
        }
        check_for_negative_theta(S_old);
    }

    // **************************************************************************************
    // Update the dycore
    // **************************************************************************************
    advance_dycore(lev, state_old, state_new,
                   U_old, V_old, W_old,
                   U_new, V_new, W_new,
                   cc_source, xmom_source, ymom_source, zmom_source, buoyancy,
                   Geom(lev), dt_lev, time);

    // **************************************************************************************
    // Tests on the reasonableness of the solution after the dycore
    // **************************************************************************************
    // Test for NaNs after dycore
    if (check_for_nans > 0) {
        if (verbose > 1) {
            amrex::Print() << "Testing new state and vels for NaNs after dycore" << std::endl;
        }
        check_state_for_nans(S_new);
        check_vels_for_nans(rU_new[lev],rV_new[lev],rW_new[lev]);
    }

    // We only test on low temp if we have a moisture model because we are protecting against
    //    the test on low temp inside the moisture models
    if (solverChoice.moisture_type != MoistureType::None) {
        if (verbose > 1) {
            amrex::Print() << "Testing on low temperature after dycore" << std::endl;
        }
        check_for_low_temp(S_new);
    } else {
        // Otherwise we will test on negative (rhotheta) coming out of the dycore
        if (verbose > 1) {
            amrex::Print() << "Testing on negative temperature after dycore" << std::endl;
        }
        check_for_negative_theta(S_new);
    }

    // **************************************************************************************
    // Update the microphysics (moisture)
    // **************************************************************************************
    if (!solverChoice.moisture_tight_coupling)
    {
        advance_microphysics(lev, S_new, dt_lev, iteration, time);

        // Test for NaNs after microphysics
        if (check_for_nans > 0) {
            amrex::Print() << "Testing new state for NaNs after advance_microphysics" << std::endl;
            check_state_for_nans(S_new);
        }
    }

    // **************************************************************************************
    // Update the fire model
    // **************************************************************************************
#ifdef ERF_USE_FIRE
    if (solverChoice.fire_type != FireType::None) {     
        // ✅ Prepare wind data for fire model
        MultiFab* U_ptr = nullptr;
        MultiFab* V_ptr = nullptr;
        MultiFab* terrain_ptr = nullptr;
        const iMultiFab* luindex_ptr = nullptr;
        
        if (solverChoice.fire_met_source == "Atmosphere") {
            // Extract surface-level winds from 3D atmospheric data
            const BoxArray& ba = S_new.boxArray();
            const DistributionMapping& dm = S_new.DistributionMap();
            
            // Create 2D MultiFabs for surface winds (persistent static variables)
            static std::unique_ptr<MultiFab> U_surface;
            static std::unique_ptr<MultiFab> V_surface;
            
            // Reallocate if grid changed
            if (!U_surface || U_surface->boxArray() != ba) {
                BoxArray ba_2d(ba.size());
                for (int i = 0; i < ba.size(); ++i) {
                    Box bx = ba[i];
                    bx.setSmall(2, 0);
                    bx.setBig(2, 0);
                    ba_2d.set(i, bx);
                }
                U_surface = std::make_unique<MultiFab>(ba_2d, dm, 1, 1);
                V_surface = std::make_unique<MultiFab>(ba_2d, dm, 1, 1);
            }
            
            // MPI DEBUG: About to extract surface wind
            if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
                Print() << "Extracting surface wind from k=1 for fire model\n";
            }

            // Extract k=0 (surface level in terrain-following coords) and average to cell centers
            for (MFIter mfi(*U_surface); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                
                const auto& u_3d = U_new.const_array(mfi);
                const auto& v_3d = V_new.const_array(mfi);
                
                const auto& u_2d = U_surface->array(mfi);
                const auto& v_2d = V_surface->array(mfi);
                
                // Extract k=1 (not k=0 which has NoSlipWall BC) and average from faces to cell centers
                // k=1 is the first cell center above ground where velocity is non-zero
                LaunchExtractSurfaceWindKernel(bx, u_3d, v_3d, u_2d, v_2d);
            }
            
            // Fill boundary for MPI ghost cells
            U_surface->FillBoundary(geom[lev].periodicity());
            V_surface->FillBoundary(geom[lev].periodicity());

            U_ptr = U_surface.get();
            V_ptr = V_surface.get();

            // DEBUG: Print wind statistics being sent to fire model (first few timesteps, then every 30s)
            static int debug_counter = 0;
            if (verbose > 0) {
                if (debug_counter < 5 || int(time) % 30 == 0) {
                    // CRITICAL: min/max are collective MPI operations, ALL ranks must call them!
                    Real u_min = U_surface->min(0);
                    Real u_max = U_surface->max(0);
                    Real v_min = V_surface->min(0);
                    Real v_max = V_surface->max(0);

                    // Only IO processor prints the results
                    if (ParallelDescriptor::IOProcessor()) {
                        Print() << "DEBUG WIND at t=" << time << " s (extracted at k=1 for fire):\n";
                        Print() << "  U: min=" << u_min << ", max=" << u_max << " m/s\n";
                        Print() << "  V: min=" << v_min << ", max=" << v_max << " m/s\n";
                    }
                    debug_counter++;
                }
            }
        }
        // If fire_met_source == "Constant", U_ptr and V_ptr remain nullptr

        // Get terrain pointer
        terrain_ptr = z_phys_cc[lev].get();
        
        // Get land use index pointer if available
        if (lev < luindex_lev.size() && !luindex_lev[lev].empty() && luindex_lev[lev][0]) {
            luindex_ptr = luindex_lev[lev][0].get();
        }
        
        // Advance fire spread model
        fire_model.Advance(lev, dt_lev, U_ptr, V_ptr, terrain_ptr, luindex_ptr);

        // DEBUG: Print fire_heat_feedback status once
        if (verbose > 0 && time < 1.0 && ParallelDescriptor::IOProcessor()) {
            Print() << "DEBUG: fire_heat_feedback = " << (solverChoice.fire_heat_feedback ? "TRUE" : "FALSE") << "\n";
        }

        // Apply fire heat flux feedback to atmosphere if enabled
        if (solverChoice.fire_heat_feedback) {
            // DEBUG: Confirm we're entering this block
            if (verbose > 0 && time < 1.0 && ParallelDescriptor::IOProcessor()) {
                Print() << "DEBUG: ENTERING heat flux feedback block (this should NOT appear if fire_heat_feedback=false)\n";
            }
            // Get heat flux from fire and apply to atmospheric state
            if (fire_heat_flux[lev]) {
                fire_model.GetHeatFlux(lev, *fire_heat_flux[lev]);
                
                // Apply fire heat flux to lowest atmospheric level
                // This couples the fire model back to the atmosphere
                for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    const auto& cons_arr = S_new.array(mfi);
                    const auto& heat_arr = fire_heat_flux[lev]->const_array(mfi);
                    
                    // ✅ SIMPLEST FIX: Just get the value directly
                    Real dz = geom[lev].CellSize(2);
                    
                    // Convert heat flux (W/m²) to source term for RhoTheta
                    // Use Cp_d from ERF_Constants.H (1004.5 J/(kg·K))
                    LaunchApplyHeatFluxKernel(bx, cons_arr, heat_arr, dt_lev, Cp_d, dz, int(RhoTheta_comp));
                }
            }

            // Apply moisture flux (latent heat) to water vapor mixing ratio
            // Only if moisture_flux MultiFab is allocated
            if (fire_moisture_flux[lev]) {
                fire_model.GetMoistureFlux(lev, *fire_moisture_flux[lev]);

                for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    const auto& cons_arr = S_new.array(mfi);
                    const auto& moist_arr = fire_moisture_flux[lev]->const_array(mfi);
                    
                    Real dz = geom[lev].CellSize(2);

                    // Apply moisture flux as source term to RhoQ (water vapor density)
                    // Moisture flux [kg/m²/s] -> source [kg/m³/s]
                    // RhoQ += (moisture_flux / dz) * dt
                    LaunchApplyMoistureFluxKernel(bx, cons_arr, moist_arr, dt_lev, dz, int(RhoQ1_comp));
                }
            }
        }
        
        if (verbose > 1) {
            Real fire_area = fire_model.get_fire_area(lev);
            Real max_heat = fire_heat_flux[lev] ? fire_heat_flux[lev]->max(0) : 0.0;
            Real max_moisture = fire_moisture_flux[lev] ? fire_moisture_flux[lev]->max(0) : 0.0;
            amrex::Print() << "  Fire at t=" << time+dt_lev 
                          << ": area=" << fire_area << " m²"
                          << ", max_heat=" << max_heat << " W/m²"
                          << ", max_moisture=" << max_moisture << " kg/m²/s\n";
        }
    }
#endif
    // ========== END FIRE MODEL ADVANCE ==========

    // **************************************************************************************
    // Update the land surface model
    // **************************************************************************************
    advance_lsm(lev, S_new, U_new, V_new, dt_lev);

#ifdef ERF_USE_PARTICLES
    // **************************************************************************************
    // Update the particle positions
    // **************************************************************************************
   evolveTracers( lev, dt_lev, vars_new, z_phys_nd );
#endif

    // ***********************************************************************************************
    // Impose domain boundary conditions here so that in FillPatching the fine data we won't
    // need to re-fill these
    // ***********************************************************************************************
    if (lev < finest_level) {
         IntVect ngvect_vels = vars_new[lev][Vars::xvel].nGrowVect();
         (*physbcs_cons[lev])(vars_new[lev][Vars::cons], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              0,vars_new[lev][Vars::cons].nComp(),
                              vars_new[lev][Vars::cons].nGrowVect(),time,BCVars::cons_bc,true);
            (*physbcs_u[lev])(vars_new[lev][Vars::xvel], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              ngvect_vels,time,BCVars::xvel_bc,true);
            (*physbcs_v[lev])(vars_new[lev][Vars::yvel], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              ngvect_vels,time,BCVars::yvel_bc,true);
            (*physbcs_w[lev])(vars_new[lev][Vars::zvel], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              ngvect_vels,time,BCVars::zvel_bc,true);
    }

    // **************************************************************************************
    // Register old and new coarse data if we are at a level less than the finest level
    // **************************************************************************************
    if (lev < finest_level) {
        if (cf_width > 0) {
            // We must fill the ghost cells of these so that the parallel copy works correctly
            state_old[IntVars::cons].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::cons].FillBoundary(geom[lev].periodicity());
            FPr_c[lev].RegisterCoarseData({&state_old[IntVars::cons], &state_new[IntVars::cons]},
                                          {time, time+dt_lev});
        }

        if (cf_width >= 0) {
            // We must fill the ghost cells of these so that the parallel copy works correctly
            state_old[IntVars::xmom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::xmom].FillBoundary(geom[lev].periodicity());
            FPr_u[lev].RegisterCoarseData({&state_old[IntVars::xmom], &state_new[IntVars::xmom]},
                                          {time, time+dt_lev});

            state_old[IntVars::ymom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::ymom].FillBoundary(geom[lev].periodicity());
            FPr_v[lev].RegisterCoarseData({&state_old[IntVars::ymom], &state_new[IntVars::ymom]},
                                          {time, time+dt_lev});

            state_old[IntVars::zmom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::zmom].FillBoundary(geom[lev].periodicity());
            FPr_w[lev].RegisterCoarseData({&state_old[IntVars::zmom], &state_new[IntVars::zmom]},
                                          {time, time+dt_lev});
        }

            //
            // Now create a MultiFab that holds (S_new - S_old) / dt from the coarse level interpolated
            //     on to the coarse/fine boundary at the fine resolution
            //
            Interpolater* mapper_f = &face_cons_linear_interp;

            // PhysBCFunctNoOp null_bc;
            // MultiFab tempx(vars_new[lev+1][Vars::xvel].boxArray(),vars_new[lev+1][Vars::xvel].DistributionMap(),1,0);
            // tempx.setVal(0.0);
            // xmom_crse_rhs[lev+1].setVal(0.0);
            // FPr_u[lev].FillSet(tempx               , time       , null_bc, domain_bcs_type);
            // FPr_u[lev].FillSet(xmom_crse_rhs[lev+1], time+dt_lev, null_bc, domain_bcs_type);
            // MultiFab::Subtract(xmom_crse_rhs[lev+1],tempx,0,0,1,IntVect{0});
            // xmom_crse_rhs[lev+1].mult(1.0/dt_lev,0,1,0);

            // MultiFab tempy(vars_new[lev+1][Vars::yvel].boxArray(),vars_new[lev+1][Vars::yvel].DistributionMap(),1,0);
            // tempy.setVal(0.0);
            // ymom_crse_rhs[lev+1].setVal(0.0);
            // FPr_v[lev].FillSet(tempy               , time       , null_bc, domain_bcs_type);
            // FPr_v[lev].FillSet(ymom_crse_rhs[lev+1], time+dt_lev, null_bc, domain_bcs_type);
            // MultiFab::Subtract(ymom_crse_rhs[lev+1],tempy,0,0,1,IntVect{0});
            // ymom_crse_rhs[lev+1].mult(1.0/dt_lev,0,1,0);

            MultiFab temp_state(zmom_crse_rhs[lev+1].boxArray(),zmom_crse_rhs[lev+1].DistributionMap(),1,0);
            InterpFromCoarseLevel(temp_state,            IntVect{0}, IntVect{0}, state_old[IntVars::zmom], 0, 0, 1,
                                  geom[lev], geom[lev+1], refRatio(lev), mapper_f, domain_bcs_type, BCVars::zvel_bc);
            InterpFromCoarseLevel(zmom_crse_rhs[lev+1],  IntVect{0}, IntVect{0}, state_new[IntVars::zmom], 0, 0, 1,
                                  geom[lev], geom[lev+1], refRatio(lev), mapper_f, domain_bcs_type, BCVars::zvel_bc);
            MultiFab::Subtract(zmom_crse_rhs[lev+1],temp_state,0,0,1,IntVect{0});
            zmom_crse_rhs[lev+1].mult(1.0/dt_lev,0,1,0);
    }

    // ***********************************************************************************************
    // Update the time averaged velocities if they are requested
    // ***********************************************************************************************
    if (solverChoice.time_avg_vel) {
        Time_Avg_Vel_atCC(dt[lev], t_avg_cnt[lev], vel_t_avg[lev].get(), U_new, V_new, W_new);
    }
}
