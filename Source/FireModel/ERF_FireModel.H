#ifndef ERF_FIRE_MODEL_H
#define ERF_FIRE_MODEL_H

#include <memory>
#include <string>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Vector.H>
#include <ERF_FireModelBase.H>
#include <ERF_NullFire.H>
#include <ERF_FireCouplingParams.H>
#include <ERF_FireGrid.H>
#include <ERF_SurfacePrecipitation.H>

class ERF;
class SurfaceLayer;

/**
 * @file ERF_FireModel.H
 * @brief The fire model container: one model per level, and the atmosphere
 *        side of the coupling shared by every model
 *
 * ERF holds one FireModel (as it holds one LandSurface), with a NullFire on
 * every level until erf.fire_model selects a model; SetModel<Model>() installs
 * it and Init() places it on the level the fire grid refines
 * (erf.fire.anchor_level). The container:
 * - checks the ERF configuration a fire grid needs and builds that grid;
 * - each step of that level, hands the model the atmospheric state
 *   (FireForcing) and advances it; with lagged or synchronous coupling, or
 *   with the smoke tracer, it averages the fluxes the model reports onto the
 *   atmospheric columns, where they wait for the next step (the one-step
 *   explicit lag of WRF-Fire);
 * - at every Runge-Kutta stage of that next step, injects the lagged heat and
 *   vapour into the cell sources (Coupling/ERF_FireAtmInjection.H) and the
 *   smoke tracer (Features/ERF_FireSmoke.H);
 * - writes the model's state and the lagged fluxes to checkpoints.
 *
 * Every per-level call is a no-op on the levels that hold a NullFire.
 * Coarser levels see the fire through average-down, finer levels not at all.
 */
class FireModel
{
public:
    /// One NullFire per level
    void ReSize (int nlev);

    /// Read and check the erf.fire.* coupling inputs
    void ReadParameters ();

    /// Install the model erf.fire_model selects; Init() places it on the fire grid's level
    template<class Model>
    void SetModel () { SetModel(std::make_unique<Model>()); }
    void SetModel (std::unique_ptr<FireModelBase> model);

    /// A model other than NullFire is installed
    bool HasModel () const { return m_has_model; }

    /// The coupling inputs
    const FireCouplingParams& Params () const { return m_params; }
    /// Test hook: set the coupling inputs without reading them
    void SetParams (const FireCouplingParams& params) { m_params = params; }

    /// The level the fire grid refines (-1 before Init)
    int Level () const { return m_lev; }

    /// The model on level lev (NullFire on every level but the fire grid's)
    FireModelBase&       Model (int lev)       { return *m_fire[lev]; }
    const FireModelBase& Model (int lev) const { return *m_fire[lev]; }

    /// A field of the model on level lev by name (FireVar), or nullptr
    const amrex::MultiFab* Get_Data_Ptr (int lev, const std::string& name) const;

    /**
     * @brief Check the configuration, build the fire grid and set the model up
     *
     * Called once from ERF::InitData_post() after every surface layer exists.
     *
     * @param[in] erf                    ERF instance (public AmrCore accessors only)
     * @param[in] lev                    The level the fire grid refines
     * @param[in] surface_layer          The bottom surface layer, or nullptr (aborts)
     * @param[in] z_phys_nd              Nodal height of level lev [m]
     * @param[in] stop_time              ERF's stop_time [s]
     * @param[in] max_step               ERF's max_step
     * @param[in] steps_per_coarse_step  Steps of level lev per level-0 step
     * @param[in] restarting             ReadCheckpoint() follows
     */
    void Init (const ERF& erf, int lev,
               const SurfaceLayer* surface_layer,
               const amrex::MultiFab& z_phys_nd,
               amrex::Real stop_time, int max_step, int steps_per_coarse_step,
               bool restarting);

    /**
     * @brief Lay the coupling out on a fire grid and place the model on level
     *        lev: what Init() does after its start-up checks and before it sets
     *        the model up
     *
     * Allocates the lagged flux buffers (zero) on the k = 0 slabs of ba_atm,
     * with 5 ghost columns. Public so that a test can drive the coupling with a
     * model of its own without an ERF object.
     */
    void SetFireGrid (int lev, const FireGrid& fg,
                      const amrex::BoxArray& ba_atm,
                      const amrex::DistributionMapping& dm_atm,
                      const amrex::Geometry& geom_atm);

    /// The fire grid
    const FireGrid& GetFireGrid () const { return m_fg; }

    /**
     * @brief Advance the fire one step of level lev and refill the lagged fluxes
     *
     * Called from ERF::Advance() once per step of every level, after the
     * level's dycore; a no-op except on the fire grid's level. The air
     * temperature and humidity at k = 0 come from the state before the dycore;
     * the wind from before the dycore (passive and lagged coupling) or after
     * it (synchronous). Optional inputs are computed only when requested.
     */
    void Advance (int lev, amrex::Real time, amrex::Real dt,
                  const amrex::MultiFab& S_old,
                  const amrex::MultiFab& xvel_old, const amrex::MultiFab& yvel_old,
                  const amrex::MultiFab& zvel_old,
                  const amrex::MultiFab& xvel_new, const amrex::MultiFab& yvel_new,
                  const amrex::MultiFab& zvel_new,
                  const amrex::MultiFab& z_phys_cc,
                  const amrex::BoxArray& ba2d,
                  const amrex::Geometry& geom,
                  const SurfacePrecipAccumulationSources& precip_sources);

    /**
     * @brief Inject the lagged fire heat, vapour and smoke into the cell sources
     *
     * Called from the slow right-hand side at every Runge-Kutta stage of every
     * level, after make_sources has rebuilt cc_source; a no-op except on the
     * fire grid's level. erf.fire.source_mode decides whether the heat replaces
     * or adds to the theta and vapour sources. Heat and vapour only with lagged
     * or synchronous coupling; smoke with erf.fire.smoke_enable.
     */
    void AddSources (int lev,
                     amrex::MultiFab& cc_source,
                     const amrex::MultiFab& S_old,
                     const amrex::MultiFab& z_phys_cc,
                     const amrex::MultiFab* z_phys_nd,
                     const amrex::Geometry& geom,
                     bool has_moisture,
                     int step);

    /// Write the model's fire-grid plotfile
    void WritePlotfile (const std::string& prefix, amrex::Real time, int step) const;

    /**
     * @brief Write the fire to a checkpoint: the model's state
     *        (FireModelBase::WriteCheckpoint), then the lagged fluxes
     *        (FireQAtmPrev and FireQLatAtmPrev) under checkpointname/Level_<lev>/
     */
    void WriteCheckpoint (const std::string& checkpointname) const;

    /**
     * @brief Why a checkpoint's lagged fluxes do not belong to this run's fire
     *        grid, or an empty string when they do
     */
    std::string CheckpointMismatch (int chk_lev, const amrex::Box& chk_region) const;

    /**
     * @brief Restore the model's state, then the lagged fluxes, on a restart
     *
     * Aborts when the model does not support restart, and when the checkpoint
     * holds the lagged fluxes on another level or over another region than this
     * run's fire grid. A checkpoint without them leaves them at zero: the first
     * restarted step injects no fire heat.
     */
    void ReadCheckpoint (const std::string& restart_chkfile, int max_level, int step);

    /**
     * @brief Sensible heat flux [W/m2] of the previous fire step on the
     *        atmospheric columns of the fire grid's level, with 5 ghost columns;
     *        zero outside the fire region. Allocated by SetFireGrid().
     */
    const amrex::MultiFab* Get_Q_atm_prev () const { return m_Q_atm_prev.get(); }
    /// Latent heat flux [W/m2] of the previous fire step, laid out as Get_Q_atm_prev().
    const amrex::MultiFab* Get_Q_lat_atm_prev () const { return m_Q_lat_atm_prev.get(); }

private:
    /// Average the model's fluxes onto the atmospheric columns for the next step.
    void UpdateFluxBuffer (const amrex::Geometry& geom_atm);

    /// The model on the fire grid's level
    FireModelBase&       ActiveModel ()       { return *m_fire[m_lev]; }
    const FireModelBase& ActiveModel () const { return *m_fire[m_lev]; }

    FireCouplingParams                            m_params;
    amrex::Vector<std::unique_ptr<FireModelBase>> m_fire;     ///< one per level: NullFire, the model on m_lev
    std::unique_ptr<FireModelBase>                m_pending;  ///< the model SetModel() installed, until it is placed
    bool                                          m_has_model = false;
    FireGrid                                      m_fg;
    int                                           m_lev = -1;
    amrex::Geometry                               m_geom_atm; ///< geometry of the fire grid's level

    // The model's fluxes averaged onto the atmospheric columns of the fire
    // grid's level (one slab per box of that level), held from the previous
    // step for injection at every stage of the next (one-step explicit lag).
    // Their ghost columns hold the neighbouring box's value for a reader that
    // samples a halo around every tile; outside the refined region of a finer
    // level they hold zero.
    std::unique_ptr<amrex::MultiFab> m_Q_atm_prev;      ///< sensible [W/m2]
    std::unique_ptr<amrex::MultiFab> m_Q_lat_atm_prev;  ///< latent   [W/m2]
};

#endif
