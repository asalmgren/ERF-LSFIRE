#ifndef ERF_FIRE_MODEL_TEMPLATE_H
#define ERF_FIRE_MODEL_TEMPLATE_H

#include <ERF_FireModelBase.H>

/**
 * @file ERF_FireModelTemplate.H
 * @brief A starting point for a new fire model: declarations only
 *
 * This class declares every FireModelBase method a new fire model provides and
 * says what ERF hands it and expects back. Nothing here has a body, no input
 * selects it and nothing creates one, so it adds no code to ERF. The stand-in
 * model of the FireModel unit tests (Tests/Unit/Fire/ERF_GTestFireCoupler.cpp)
 * is the smallest model that works.
 *
 * To add a model (as for a land surface or wind farm model):
 *  1. Create Source/FireModel/<Model>/, copy this file there, rename the class
 *     and the include guard, and write the bodies in a .cpp file there. Delete
 *     the declarations of the optional methods the model does not provide
 *     (their defaults in FireModelBase then apply); a declaration left without
 *     a body fails to link.
 *  2. Read the model's inputs under its own prefix (erf.fire.<model>.*), so they
 *     never mix with the coupling's keys (FireCouplingParams).
 *  3. Add the model to FireType (ERF_DataStruct.H) and install it in
 *     ERF::ReadParameters() with fire_model.SetModel<Model>() under the model's
 *     build define; add a Make.package and the CMake sources, as for the other
 *     model directories.
 *  4. Run the FireModel unit tests with the new model in place of their
 *     stand-in, and a restart of a coupled run against a straight run (when the
 *     model supports restart): the atmosphere and fire plotfiles must agree at
 *     zero tolerance.
 *
 * What the container already does for every model: it builds the fire grid
 * over erf.fire.anchor_level, hands the model the atmosphere's state each
 * step, averages the model's fluxes onto the atmospheric columns, holds them
 * one step, injects heat, vapour and smoke, and saves those lagged fluxes in
 * the checkpoint. A model owns only its fire physics and its own state.
 *
 * The call order is the one documented on FireModelBase (ERF_FireModelBase.H).
 * ERF_FireModel.cpp checks at compile time that this class declares every
 * method FireModelBase requires (the pure virtual ones), and override checks
 * the signature of every method declared here. The optional methods are listed
 * by hand: a new optional method in FireModelBase is added here too.
 */
class FireModelTemplate : public FireModelBase
{
public:
    /// The model's name, for messages
    const char* Name () const override;

    /**
     * Read the model's own inputs and check them against the setup, aborting
     * naming the input on anything the model cannot run (a missing fuel file, a
     * grid ratio it does not support). Allocate the model's fields on the fire
     * grid setup.fg (fg.ba, fg.dm, fg.geom; fg.atm_lo is the fire region's corner
     * in the level's index space) and set the initial fire: fuels, terrain slope,
     * ignition. setup.stop_time and setup.max_step give the run's end.
     */
    void Init (const FireSetup& setup) override;

    /**
     * Advance the fire over one step of the fire grid's level (forcing.dt from
     * forcing.time). forcing holds the face velocities, the cell-centre heights
     * and the optional inputs RequestedInputs() asks for; sample them onto the
     * fire grid here and keep no pointer to them after this call.
     */
    void Advance (const FireForcing& forcing) override;

    /// Write the model's plotfile of the fire grid, named prefix plus the step.
    void WritePlotfile (const std::string& prefix, amrex::Real time, int step) const override;

    /**
     * Optional, nullptr unless declared: the model's fields by name (FireVar).
     * FireVar::sensible_flux and FireVar::latent_flux are the step's mean heat
     * fluxes into the air on the fire grid [W/m2] (nullptr for none); they must
     * stay valid until the next Advance(). FireVar::open_fraction_atm and
     * FireVar::roof_height_atm serve erf.fire.heat_open_fraction, on the level's
     * k = 0 slabs with one ghost cell.
     */
    const amrex::MultiFab* FireDataPtr (const std::string& name) const override;

    /// Optional, FireInput_None unless declared: the optional inputs Advance() reads.
    unsigned RequestedInputs () const override;

    /**
     * Optional, 0 unless declared: heat released per kilogram of fuel burned
     * [J/kg], matching the sensible flux last reported; the smoke tracer needs it.
     */
    amrex::Real SmokeHeatPerKg () const override;

    /// Optional, false unless declared: the model saves and restores its state.
    bool SupportsRestart () const override;

    /**
     * Write everything the model carries from one step to the next under
     * checkpointname/Level_<lev>/ (lev the fire grid's level). The container
     * writes its lagged fluxes after this; use names that do not start with FireQ.
     */
    void WriteCheckpoint (const std::string& checkpointname) const override;

    /**
     * Read back what WriteCheckpoint() wrote, on a restart, after Init() and
     * before the first step; step is the fire grid level's step count at the
     * checkpoint. Abort if the checkpoint's fire grid is not this run's.
     */
    void ReadCheckpoint (const std::string& restart_chkfile, int step) override;
};

#endif
