#ifndef ERF_FIRE_MODEL_BASE_H
#define ERF_FIRE_MODEL_BASE_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <string>
#include <ERF_FireGrid.H>

class SurfaceLayer;

/**
 * @file ERF_FireModelBase.H
 * @brief The interface between ERF and a fire model
 *
 * A fire model advances a fire front on the fire grid and reports the heat it
 * releases. Everything that touches the atmosphere belongs to the FireModel
 * container (ERF_FireModel.H): the fire grid's placement, the atmospheric
 * fields handed to the model, the area average of its fluxes onto the
 * atmospheric columns, their one-step lag, their injection into the cell
 * sources and the smoke tracer. A model sees only this interface and the
 * helpers it chooses to use (Tools/).
 *
 * The pattern follows ERF's land surface and wind farm models: erf.fire_model
 * selects a model type, ERF::ReadParameters() installs it with
 * FireModel::SetModel<Model>(), and every level that does not carry the fire
 * grid holds a NullFire. ERF_FireModelTemplate.H declares every method below
 * for a new model.
 *
 * Zero-copy contract. The forcing fields are references to ERF's own
 * MultiFabs of the level the fire grid refines (its index space, box order,
 * DistributionMapping and arena; device memory on GPU builds), and the fields
 * a model reports through FireDataPtr() are the model's own MultiFabs, read
 * in place. Fire box n refines atmospheric box n and has the same owner, so
 * one MFIter serves both and no data moves between ranks.
 */

/// Everything a fire model is given once at start-up
struct FireSetup
{
    int                    lev;                    ///< ERF level the fire grid refines (erf.fire.anchor_level, resolved)
    const FireGrid&        fg;                     ///< the fire grid (ba, dm, geom, C, atm_lo)
    const amrex::Geometry& geom_atm;               ///< geometry of level lev
    const amrex::MultiFab& z_phys_nd;              ///< nodal heights of level lev [m]
    const SurfaceLayer*    surface_layer;          ///< the bottom surface layer
    amrex::Real            stop_time;              ///< ERF stop_time [s] (numeric_limits max when unset)
    int                    max_step;               ///< ERF max_step, counting level-0 steps (negative when unset)
    int                    steps_per_coarse_step;  ///< steps of level lev per level-0 step
    bool                   restarting;             ///< ReadCheckpoint() follows Init()
};

/// Optional atmospheric inputs: the coupling computes only those a model requests
enum FireInput : unsigned {
    FireInput_None          = 0,
    FireInput_AirTemp       = 1u << 0,  ///< FireForcing::T_k0
    FireInput_RelHum        = 1u << 1,  ///< FireForcing::RH_k0
    FireInput_SurfacePrecip = 1u << 2,  ///< FireForcing::precip_accum_k0 (needs a precipitating microphysics)
    FireInput_VerticalVel   = 1u << 3   ///< FireForcing::zvel
};

/**
 * @brief The atmospheric state a fire model reads at the start of its step
 *
 * Every field lives on the level the fire grid refines, in that level's index
 * space and with its ghost cells; the model samples them onto its own grid
 * (FireGrid::atm_lo gives the corner of the fire region). The k = 0 slabs hold
 * one cell in z per box of the level. Valid only during Advance(): a model
 * keeps no pointer to them.
 */
struct FireForcing
{
    amrex::Real            time;        ///< start of the step [s]
    amrex::Real            dt;          ///< step of the fire grid's level [s]
    const amrex::MultiFab& xvel;        ///< x-face velocity [m/s]: before the dycore (passive, lagged) or after it (synchronous)
    const amrex::MultiFab& yvel;        ///< y-face velocity [m/s], the same time level as xvel
    const amrex::MultiFab& z_phys_cc;   ///< cell-centre height [m] of the terrain-following mesh
    // Optional inputs, nullptr unless requested (RequestedInputs())
    const amrex::MultiFab* zvel            = nullptr;  ///< z-face velocity [m/s], the same time level as xvel
    const amrex::MultiFab* T_k0            = nullptr;  ///< air temperature at k = 0 [K] (theta times Exner), before the dycore
    const amrex::MultiFab* RH_k0           = nullptr;  ///< relative humidity at k = 0 [0-1], before the dycore
    const amrex::MultiFab* precip_accum_k0 = nullptr;  ///< surface precipitation accumulated since t = 0 [kg/m2], after this step's microphysics
};

/**
 * @brief Names of the fields a fire model may report through FireDataPtr()
 *
 * The coupling and the optional features look the fields up by these names; a
 * model that has no such field returns nullptr. Flux fields live on the fire
 * grid, per unit ground area, as the mean over the step just taken, and stay
 * valid until the next Advance().
 */
namespace FireVar {
    inline const std::string sensible_flux     = "sensible_flux";      ///< sensible heat into the air [W/m2], fire grid
    inline const std::string latent_flux       = "latent_flux";        ///< latent heat into the air [W/m2], fire grid
    inline const std::string level_set         = "level_set";          ///< level set, < 0 burned, fire grid
    inline const std::string open_fraction_atm = "open_fraction_atm";  ///< fraction of a column open to the air [0-1], level k = 0 slabs, 1 ghost
    inline const std::string roof_height_atm   = "roof_height_atm";    ///< mean roof height above the ground [m], level k = 0 slabs, 1 ghost
}

/**
 * @brief A fire model behind the coupling
 *
 * The FireModel container calls, in this order:
 * - Init() once, with the fire grid built (and ReadCheckpoint() right after it
 *   on a restart);
 * - Advance() once per step of the fire grid's level, after that level's
 *   dycore, then FireDataPtr() for the fluxes when the coupling injects heat
 *   or smoke;
 * - SmokeHeatPerKg() and FireDataPtr() for the building fields at every
 *   Runge-Kutta stage of the next step;
 * - WritePlotfile() on the fire plotfile schedule, WriteCheckpoint() when ERF
 *   writes a checkpoint.
 */
class FireModelBase
{
public:
    virtual ~FireModelBase () = default;

    // ===================== required =====================

    /// Name of the model, for messages
    virtual const char* Name () const = 0;

    /**
     * @brief Set the model up on the fire grid
     *
     * Read the model's own inputs (under its own prefix, erf.fire.<model>.*),
     * check them against the setup and abort naming the input on anything the
     * model cannot run; allocate the model's fields on setup.fg and set the
     * initial fire (fuels, terrain slope, ignition).
     */
    virtual void Init (const FireSetup& setup) = 0;

    /// Advance the fire over [forcing.time, forcing.time + forcing.dt]
    virtual void Advance (const FireForcing& forcing) = 0;

    /// Write the model's plotfile of the fire grid (prefix plus the level-0 step)
    virtual void WritePlotfile (const std::string& prefix, amrex::Real time, int step) const = 0;

    // ===================== optional =====================

    /// A field of the model by name (FireVar), or nullptr when the model has none
    virtual const amrex::MultiFab* FireDataPtr (const std::string& /*name*/) const { return nullptr; }

    /// Bitwise OR of the FireInput flags the model reads
    virtual unsigned RequestedInputs () const { return FireInput_None; }

    /// Heat released per kilogram of fuel burned [J/kg], matching the sensible
    /// flux last reported; the smoke tracer needs it (0 = none)
    virtual amrex::Real SmokeHeatPerKg () const { return amrex::Real(0.0); }

    /// The model writes and reads its state in checkpoints; without it a restart
    /// with the fire running stops before reading the checkpoint
    virtual bool SupportsRestart () const { return false; }

    /// Everything the model carries from one step to the next, under
    /// checkpointname/Level_<lev>/ (names must not start with FireQ)
    virtual void WriteCheckpoint (const std::string& /*checkpointname*/) const {}

    /// Read back what WriteCheckpoint() wrote; step is the fire grid level's step
    /// count at the checkpoint. Abort if the checkpoint's fire grid is not this run's.
    virtual void ReadCheckpoint (const std::string& /*restart_chkfile*/, int /*step*/) {}
};

#endif
