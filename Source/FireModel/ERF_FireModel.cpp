#include <ERF_FireModel.H>
#include <ERF_FireModelTemplate.H>
#include <ERF_FireAtmInjection.H>
#include <ERF_FirePrerequisites.H>
#include <ERF_FireAtmState.H>
#include <ERF_FireSmoke.H>
#include <ERF_IndexDefines.H>
#include <ERF.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <sstream>
#include <type_traits>
#include <utility>

using namespace amrex;

// The template for new models declares every method FireModelBase requires: a
// pure virtual method added to the interface and missing there makes the
// template abstract. Optional methods (with a default body) are not caught
// here; the interface asks for them to be added to the template by hand.
// Compile time only; no FireModelTemplate is ever created.
static_assert(!std::is_abstract_v<FireModelTemplate>,
              "ERF_FireModelTemplate.H must declare every pure virtual method of FireModelBase");

namespace {

// Where the fire grid sits, printed once at start-up, and the warning when a
// finer level covers part of it: the fire then reads that level's coarser wind,
// and its heat on this level is replaced by the finer level's average-down.
void report_fire_region(const ERF& erf, int lev, const FireGrid& fg)
{
    const Geometry& g = fg.geom;
    amrex::Print() << "[FIRE] Fire grid on level " << lev << " of 0-" << erf.finestLevel()
                   << ": x " << g.ProbLo(0) << " to " << g.ProbHi(0)
                   << " m, y " << g.ProbLo(1) << " to " << g.ProbHi(1) << " m, "
                   << g.Domain().length(0) << " x " << g.Domain().length(1) << " fire cells of "
                   << g.CellSize(0) << " x " << g.CellSize(1) << " m\n";

    if (lev >= erf.finestLevel()) { return; }

    // Share of the fire region that level lev + 1 covers, in columns of level lev.
    const Box region(IntVect(fg.atm_lo[0], fg.atm_lo[1], 0),
                     IntVect(fg.atm_lo[0] + g.Domain().length(0) / fg.C - 1,
                             fg.atm_lo[1] + g.Domain().length(1) / fg.C - 1, 0));
    BoxArray fine = erf.boxArray(lev + 1);
    fine.coarsen(erf.refRatio(lev));
    amrex::Long covered = 0;
    for (int n = 0; n < fine.size(); ++n) {
        Box b = fine[n];
        b.setSmall(2, 0);
        b.setBig(2, 0);
        b &= region;
        if (b.ok()) { covered += b.numPts(); }
    }
    const amrex::Real share = amrex::Real(covered) / amrex::Real(region.numPts());
    amrex::Print() << "[FIRE] WARNING: the fire grid is on level " << lev << " but the finest level is "
                   << erf.finestLevel() << " (erf.fire.anchor_level = " << lev << "), and level " << lev + 1
                   << " covers " << 100.0 * share << " % of the fire grid. The fire reads the wind of level "
                   << lev << " and its heat, moisture and smoke go into level " << lev << " only: with "
                   << "erf.coupling_type = TwoWay the average-down from the finer level replaces the heated cells "
                   << "under it, so that heat is lost, and with OneWay the finer level never sees the fire. "
                   << "Leave erf.fire.anchor_level unset to put the fire grid on the finest level.\n";
}

} // namespace

void
FireModel::ReSize (int nlev)
{
    m_fire.resize(nlev);
    for (auto& f : m_fire) { f = std::make_unique<NullFire>(); }
}

void
FireModel::ReadParameters ()
{
    m_params = FireCouplingParams::read();
}

void
FireModel::SetModel (std::unique_ptr<FireModelBase> model)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(model != nullptr, "[FIRE] FireModel::SetModel needs a model");
    m_pending   = std::move(model);
    m_has_model = true;
}

const MultiFab*
FireModel::Get_Data_Ptr (int lev, const std::string& name) const
{
    return (lev >= 0 && lev < static_cast<int>(m_fire.size())) ? m_fire[lev]->FireDataPtr(name) : nullptr;
}

void
FireModel::Init (const ERF& erf, int lev,
                 const SurfaceLayer* surface_layer,
                 const MultiFab& z_phys_nd,
                 Real stop_time, int max_step, int steps_per_coarse_step,
                 bool restarting)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_has_model, "[FIRE] FireModel::Init without a model (SetModel)");

    // The ERF configuration a fire grid needs; the model checks its own inputs in its Init
    verify_fire_prerequisites(erf, lev, surface_layer, m_params.grid_ratio);

    FireGrid fg = create_fire_grid(erf.boxArray(lev), erf.DistributionMap(lev),
                                   erf.Geom(lev), m_params.grid_ratio);
    fg.lev = lev;
    report_fire_region(erf, lev, fg);
    SetFireGrid(lev, fg, erf.boxArray(lev), erf.DistributionMap(lev), erf.Geom(lev));

    const FireSetup setup{lev, m_fg, erf.Geom(lev), z_phys_nd, surface_layer,
                          stop_time, max_step, steps_per_coarse_step, restarting};
    ActiveModel().Init(setup);

    // The smoke tracer turns the sensible flux into the fuel burnt with the
    // model's heat per kilogram
    if (m_params.smoke_enable && !(ActiveModel().SmokeHeatPerKg() > 0.0_rt)) {
        Abort(std::string("[FIRE] erf.fire.smoke_enable = true, but the fire model (") + ActiveModel().Name()
              + ") reports no heat per kilogram of fuel burnt (FireModelBase::SmokeHeatPerKg)");
    }

    // erf.fire.heat_open_fraction places the heat by the open fraction the
    // model reports; without one the heat is placed as if every column were open
    if (m_params.heat_open_fraction && m_params.injects_flux() &&
        (!ActiveModel().FireDataPtr(FireVar::open_fraction_atm) ||
         !ActiveModel().FireDataPtr(FireVar::roof_height_atm))) {
        amrex::Print() << "[FIRE] WARNING: erf.fire.heat_open_fraction = true, but the fire model ("
                       << ActiveModel().Name() << ") reports no open fraction;"
                       << " the heat is placed as if every column were open\n";
    }
}

void
FireModel::SetFireGrid (int lev, const FireGrid& fg,
                        const BoxArray& ba_atm,
                        const DistributionMapping& dm_atm,
                        const Geometry& geom_atm)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_has_model, "[FIRE] FireModel::SetFireGrid without a model (SetModel)");
    if (static_cast<int>(m_fire.size()) <= lev) {
        const int n0 = static_cast<int>(m_fire.size());
        m_fire.resize(lev + 1);
        for (int l = n0; l <= lev; ++l) { m_fire[l] = std::make_unique<NullFire>(); }
    }
    // Place the installed model on the fire grid's level, once
    if (m_pending) {
        m_fire[lev] = std::move(m_pending);
    } else {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev == m_lev, "[FIRE] FireModel: the fire grid's level cannot change");
    }

    m_lev      = lev;
    m_fg       = fg;
    m_geom_atm = geom_atm;
    const int nz = geom_atm.Domain().length(2);

    // The lagged fluxes live on the k = 0 slabs of the level's boxes, with the
    // level's owners. Ghost columns for a reader that samples the fluxes in the
    // halo around every tile: 5, the most ERF's state carries
    // (ComputeGhostCells() + 1); UpdateFluxBuffer() fills them.
    BoxArray ba_atm_2d = ba_atm;
    ba_atm_2d.coarsen(IntVect(1, 1, nz));
    const IntVect ng_flux(5, 5, 0);
    m_Q_atm_prev     = std::make_unique<MultiFab>(ba_atm_2d, dm_atm, 1, ng_flux);
    m_Q_lat_atm_prev = std::make_unique<MultiFab>(ba_atm_2d, dm_atm, 1, ng_flux);
    m_Q_atm_prev->setVal(0.0_rt);
    m_Q_lat_atm_prev->setVal(0.0_rt);
}

void
FireModel::Advance (int lev, Real time, Real dt,
                    const MultiFab& S_old,
                    const MultiFab& xvel_old, const MultiFab& yvel_old,
                    const MultiFab& zvel_old,
                    const MultiFab& xvel_new, const MultiFab& yvel_new,
                    const MultiFab& zvel_new,
                    const MultiFab& z_phys_cc,
                    const BoxArray& ba2d,
                    const Geometry& geom,
                    const SurfacePrecipAccumulationSources& precip_sources)
{
    if (lev != m_lev) { return; }

    const unsigned inputs = ActiveModel().RequestedInputs();

    // Air temperature and RH at k=0 (a model may drive a fuel moisture model
    // with them), derived from the pre-dycore state. This is used in both
    // lagged and synchronous modes: the moisture ODE uses time-averaged T, and
    // the change in T from fire heating within one atmospheric step is small
    // relative to the moisture time lag. The temperature is theta * Exner(p)
    // (the moisture curves and ignition thresholds take an air temperature,
    // not the k = 0 potential temperature, which is 14-16 K warmer at 840 hPa).
    std::unique_ptr<MultiFab> T_atm_k0, RH_atm_k0;
    if (inputs & FireInput_AirTemp) {
        T_atm_k0 = std::make_unique<MultiFab>(ba2d, S_old.DistributionMap(), 1, 0);
        compute_t_from_conservative(*T_atm_k0, S_old);
    }
    if (inputs & FireInput_RelHum) {
        RH_atm_k0 = std::make_unique<MultiFab>(ba2d, S_old.DistributionMap(), 1, 0);
        compute_rh_from_conservative(*RH_atm_k0, S_old, geom);
    }

    // Surface precipitation accumulation of every column after this step's
    // microphysics, for a model that takes the rain per column; the model
    // differences it against the accumulation it saw at its previous step.
    std::unique_ptr<MultiFab> precip_accum_k0;
    if (inputs & FireInput_SurfacePrecip) {
        precip_accum_k0 = std::make_unique<MultiFab>(ba2d, S_old.DistributionMap(), 1, 0);
        fire_surface_precip_accum_k0(*precip_accum_k0, precip_sources, geom.Domain().smallEnd(2));
    }

    // In passive and lagged modes, fire uses pre-dycore wind (vars_old).
    // In synchronous mode, fire uses post-dycore wind (vars_new) so that
    // the fire spread rate at step n+1 reflects the atmospheric momentum
    // response to fire heating from step n.
    const bool post_dycore = m_params.is_synchronous();
    if (m_params.fire_debug) {
        const char* mode = m_params.is_synchronous() ? "SYNCHRONOUS coupling with POST-dycore"
                         : m_params.is_lagged()      ? "LAGGED coupling with PRE-dycore"
                                                     : "PASSIVE coupling with PRE-dycore";
        amrex::Print() << "[FIRE DEBUG] Fire advance using " << mode << " wind at t="
                       << time << ", dt=" << dt << std::endl;
    }
    FireForcing forcing{time, dt,
                        post_dycore ? xvel_new : xvel_old,
                        post_dycore ? yvel_new : yvel_old,
                        z_phys_cc};
    if (inputs & FireInput_VerticalVel) { forcing.zvel = post_dycore ? &zvel_new : &zvel_old; }
    forcing.T_k0            = T_atm_k0.get();
    forcing.RH_k0           = RH_atm_k0.get();
    forcing.precip_accum_k0 = precip_accum_k0.get();
    ActiveModel().Advance(forcing);

    // Store the fire flux of this step for injection at the next step via the
    // Runge-Kutta source term (AddSources). This applies in both lagged and
    // synchronous modes; in passive mode the buffer keeps its zero unless the
    // smoke tracer needs it.
    UpdateFluxBuffer(geom);
}

void
FireModel::UpdateFluxBuffer (const Geometry& geom_atm)
{
    if (!m_params.injects_flux() && !m_params.smoke_enable) {
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Skipping flux buffer update: coupling_type is passive (injects_flux=false) and smoke_enable=false" << std::endl;
        }
        return;
    }

    if (!m_Q_atm_prev) { return; }
    const MultiFab* sensible = ActiveModel().FireDataPtr(FireVar::sensible_flux);
    const MultiFab* latent   = ActiveModel().FireDataPtr(FireVar::latent_flux);
    if (!sensible) {
        // No flux this step: last step's must not be injected again
        m_Q_atm_prev->setVal(0.0_rt);
        m_Q_lat_atm_prev->setVal(0.0_rt);
        return;
    }

    coarsen_fire_flux_to_atm(*m_Q_atm_prev, *sensible, geom_atm, m_fg);
    if (m_params.fire_debug) {
        amrex::Print() << "[FIRE DEBUG] Sensible flux to the atmosphere: max " << m_Q_atm_prev->max(0)
                       << " W/m2" << std::endl;
    }

    if (m_params.inject_latent && latent) {
        coarsen_fire_flux_to_atm(*m_Q_lat_atm_prev, *latent, geom_atm, m_fg);
    } else {
        if (m_params.fire_debug && m_params.injects_flux()) {
            amrex::Print() << "[FIRE DEBUG] Skipping latent heat injection (inject_latent=false or the model reports no latent flux)" << std::endl;
        }
        m_Q_lat_atm_prev->setVal(0.0_rt);
    }

    // A reader may sample these fluxes in the halo columns around every tile;
    // give them the neighbouring box's value (and the edge column's
    // outside a non-periodic face), so a column sees the same flux whichever box
    // holds it. On a finer level the columns outside the refined region have no
    // fire grid and keep the zero they were allocated with.
    fire_fill_boundary(*m_Q_atm_prev, geom_atm);
    fire_fill_boundary(*m_Q_lat_atm_prev, geom_atm);
}

void
FireModel::AddSources (int lev,
                       MultiFab& cc_source,
                       const MultiFab& S_old,
                       const MultiFab& z_phys_cc,
                       const MultiFab* z_phys_nd,
                       const Geometry& geom,
                       bool has_moisture,
                       int step)
{
    if (lev != m_lev) { return; }

    // Heat and vapour, with the WRF-Fire profile (Coupling/ERF_FireAtmInjection.H)
    if (m_params.injects_flux() && m_Q_atm_prev && m_params.fire_atm_feedback > 0.0_rt) {
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Applying fire coupling to atmosphere (coupling_type="
                           << m_params.coupling_type << ", feedback=" << m_params.fire_atm_feedback
                           << ", max_Q_prev=" << m_Q_atm_prev->max(0) << " W/m2)" << std::endl;
        }

        const MultiFab* Q_lat_ptr = (m_params.inject_latent && has_moisture)
            ? m_Q_lat_atm_prev.get() : nullptr;

        apply_fire_tendency_to_cc_source(
            cc_source,
            *m_Q_atm_prev,
            Q_lat_ptr,
            z_phys_cc,
            z_phys_nd,
            S_old,
            geom,
            m_params.heat_flux_alfg,
            m_params.fire_atm_feedback,
            has_moisture,
            m_params.fire_debug,
            m_params.source_mode == "add",
            ActiveModel().FireDataPtr(FireVar::open_fraction_atm),
            ActiveModel().FireDataPtr(FireVar::roof_height_atm),
            m_params.heat_open_fraction,
            m_params.heat_tendency_density,
            m_params.heat_tendency_exner);
    }

    // The smoke tracer, after the heat: proportional to the fuel burnt, which
    // the model's heat per kilogram recovers from the lagged sensible flux
    if (m_params.smoke_enable && m_Q_atm_prev) {
        inject_smoke_from_fire(
            cc_source,
            *m_Q_atm_prev,
            z_phys_nd,
            geom,
            m_params.smoke_emission_factor,
            ActiveModel().SmokeHeatPerKg(),
            RhoSmoke_comp,
            m_params.fire_debug,
            step);
    }
}

void
FireModel::WritePlotfile (const std::string& prefix, Real time, int step) const
{
    if (m_lev < 0) { return; }
    ActiveModel().WritePlotfile(prefix, time, step);
}

void
FireModel::WriteCheckpoint (const std::string& checkpointname) const
{
    if (m_lev < 0) { return; }

    // The model's own state first, as ReadCheckpoint restores it
    if (ActiveModel().SupportsRestart()) {
        ActiveModel().WriteCheckpoint(checkpointname);
    } else {
        static bool warned = false;
        if (!warned) {
            amrex::Print() << "[FIRE] WARNING: the fire model (" << ActiveModel().Name()
                           << ") does not save its state: a restart from this checkpoint "
                           << "with the fire running will stop\n";
            warned = true;
        }
    }

    // The lagged fluxes injected into the atmosphere on the next step. Without
    // them the first step after a restart injects no fire heat at all.
    if (m_Q_atm_prev) {
        amrex::Print() << "Writing fire atmosphere flux buffer to checkpoint" << std::endl;
        VisMF::Write(*m_Q_atm_prev, MultiFabFileFullPrefix(m_lev, checkpointname, "Level_", "FireQAtmPrev"));
    }
    if (m_Q_lat_atm_prev) {
        VisMF::Write(*m_Q_lat_atm_prev, MultiFabFileFullPrefix(m_lev, checkpointname, "Level_", "FireQLatAtmPrev"));
    }
}

std::string
FireModel::CheckpointMismatch (int chk_lev, const Box& chk_region) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_Q_atm_prev, "[FIRE] FireModel: SetFireGrid() before the checkpoint");
    // The fluxes cover the columns of the level's boxes, which form the fire
    // region (verify_fire_prerequisites, check 7b)
    const Box run_region = m_Q_atm_prev->boxArray().minimalBox();
    if (chk_lev == m_lev && chk_region == run_region) { return {}; }
    std::ostringstream msg;
    msg << "holds the lagged fire fluxes on level " << chk_lev << " over the columns " << chk_region
        << ", but this run's fire grid covers the columns " << run_region << " of level " << m_lev
        << ". Restart with the erf.fire.anchor_level and the refinement the checkpoint was written with.";
    return msg.str();
}

void
FireModel::ReadCheckpoint (const std::string& restart_chkfile, int max_level, int step)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_Q_atm_prev, "[FIRE] FireModel: SetFireGrid() before the checkpoint");
    if (!ActiveModel().SupportsRestart()) {
        Abort(std::string("[FIRE] The fire model (") + ActiveModel().Name() + ") does not support "
              "restart: its state is not in the checkpoint " + restart_chkfile
              + ". Start the run from its initial conditions.");
    }

    // The model first: it may check, naming its own fields, that the
    // checkpoint's fire grid sits on this level and region
    ActiveModel().ReadCheckpoint(restart_chkfile, step);

    // The lagged fluxes, on the level of the fire grid. A checkpoint without them
    // (written with no fire running): the buffers keep their zero, said once.
    // Decided on the I/O rank and broadcast, so every rank takes the same branch
    // into the collective reads below
    int chk_lev = -1;
    if (ParallelDescriptor::IOProcessor()) {
        for (int l = 0; l <= max_level; ++l) {
            if (amrex::FileExists(restart_chkfile + "/Level_" + std::to_string(l) + "/FireQAtmPrev_H")) {
                chk_lev = l;
                break;
            }
        }
    }
    ParallelDescriptor::Bcast(&chk_lev, 1, ParallelDescriptor::IOProcessorNumber());
    if (chk_lev < 0) {
        amrex::Print() << "[FIRE] Checkpoint has no lagged fire fluxes; the first restarted step injects none.\n";
        return;
    }
    // Whatever the model: the fluxes belong to this run's fire grid, on its
    // level and over its region
    const Box chk_region = VisMF(MultiFabFileFullPrefix(chk_lev, restart_chkfile, "Level_", "FireQAtmPrev"))
                               .boxArray().minimalBox();
    if (const std::string why = CheckpointMismatch(chk_lev, chk_region); !why.empty()) {
        amrex::Abort("[FIRE] The checkpoint " + restart_chkfile + " " + why);
    }
    for (const auto& [mf, name] : {std::pair<MultiFab*, const char*>{m_Q_atm_prev.get(), "FireQAtmPrev"},
                                   std::pair<MultiFab*, const char*>{m_Q_lat_atm_prev.get(), "FireQLatAtmPrev"}}) {
        if (mf == nullptr) { continue; }
        // decided on one rank and broadcast, so every rank takes the same
        // path into the collective read below
        const std::string header = restart_chkfile + "/Level_" + std::to_string(m_lev) + "/" + name + "_H";
        int has_field = 0;
        if (ParallelDescriptor::IOProcessor()) { has_field = amrex::FileExists(header) ? 1 : 0; }
        ParallelDescriptor::Bcast(&has_field, 1, ParallelDescriptor::IOProcessorNumber());
        if (has_field == 0) {
            amrex::Print() << "[FIRE] Checkpoint has no " << name << "; keeping the initialized values.\n";
            continue;
        }
        VisMF::Read(*mf, MultiFabFileFullPrefix(m_lev, restart_chkfile, "Level_", name));
        // VisMF restores the valid cells (a checkpoint may hold no ghost
        // columns), so refill the ghosts the way UpdateFluxBuffer() does.
        fire_fill_boundary(*mf, m_geom_atm);
    }
}
