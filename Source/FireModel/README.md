# Fire model interface

`Source/FireModel` couples a wildland fire model to ERF. It follows the pattern
of ERF's land surface and wind farm models: an input selects the model, ERF
holds a container with one model per level, and a model plugs in by deriving
from a base class. No fire model is built into this directory; the container
and its coupling to the atmosphere are what every model shares.

## Layout

```
Source/FireModel/
├── ERF_FireModel.H/.cpp       container: one model per level + the coupling state
├── ERF_FireModelBase.H        interface every model derives from
├── ERF_FireModelTemplate.H    declarations-only starting point for a new model
├── Null/       ERF_NullFire.H                 no fire (and every level off the fire grid)
├── Coupling/   atmosphere side, used by the container
│   ├── ERF_FireGrid.H/.cpp            the fire grid over one ERF level
│   ├── ERF_FireCouplingParams.H/.cpp  erf.fire.* keys, range-checked
│   ├── ERF_FirePrerequisites.H/.cpp   start-up checks of the ERF configuration
│   └── ERF_FireAtmInjection.H         flux coarsening + vertical heat profile
├── Tools/      helpers a model may call
│   ├── ERF_FireWindExtract.H/.cpp     wind sampled onto the fire grid
│   ├── ERF_FireTerrain.H/.cpp         ground height of the atmospheric columns
│   └── ERF_FireAtmState.H             optional inputs: T, RH, surface precipitation
├── Features/   optional, switched by inputs
│   └── ERF_FireSmoke.H/.cpp           smoke tracer
└── <Model>/    one directory per fire model
```

Dependencies: the container uses `Coupling/`, `Tools/` and `Features/`; a model
sees only `ERF_FireModelBase.H` and the `Tools/` it chooses.

## How ERF drives it

ERF holds `FireModel fire_model;` (as it holds `LandSurface lsm;`).

| ERF call site | Container call | What happens |
|---|---|---|
| `ERF::initializeFire` (constructor) | `ReSize`, `ReadParameters`, `SetModel<Model>()` | `NullFire` on every level; `erf.fire.*` checked; the model `erf.fire_model` selects is installed |
| `ERF::InitData_post` | `Init`, `ReadCheckpoint` on restart | ERF checks, fire grid built on `erf.fire.anchor_level`, model placed on that level and `Init(FireSetup)` called |
| `ERF::Advance`, after the dycore | `Advance(lev, ...)` | model advanced once per step of its level; its fluxes averaged onto the atmospheric columns and held for the next step |
| slow RHS, every Runge-Kutta stage | `AddSources(lev, ...)` | held heat and vapour injected with an exponential profile; smoke injected |
| plotfile / checkpoint output | `WritePlotfile`, `WriteCheckpoint` | model's fire-grid plotfile; model state + held fluxes |

Every per-level call is a no-op on levels holding a `NullFire`.

## The model interface (`FireModelBase`)

Required:

| Method | Purpose |
|---|---|
| `Name()` | name for messages |
| `Init(const FireSetup&)` | read own inputs (`erf.fire.<model>.*`), check them, allocate fields on the fire grid, set fuels / terrain / ignition |
| `Advance(const FireForcing&)` | advance the fire over `[time, time + dt]` |
| `WritePlotfile(prefix, time, step)` | the model's fire-grid plotfile |

Optional (defaults in the base class):

| Method | Default | Purpose |
|---|---|---|
| `FireDataPtr(name)` | `nullptr` | the model's fields by name (`FireVar`): `sensible_flux`, `latent_flux` [W/m²], `level_set`, `open_fraction_atm`, `roof_height_atm` |
| `RequestedInputs()` | none | optional forcing: `FireInput_AirTemp`, `_RelHum`, `_SurfacePrecip`, `_VerticalVel` |
| `SmokeHeatPerKg()` | 0 | heat per kg of fuel burnt, needed by the smoke tracer |
| `SupportsRestart()` | false | without it a restart with the fire running stops |
| `WriteCheckpoint` / `ReadCheckpoint` | no-op | the model's own state |

What the model receives:

- `FireSetup` (once): ERF level, fire grid (`ba`, `dm`, `geom`, ratio `C`, region corner `atm_lo`), level geometry, nodal heights, surface layer, `stop_time`, `max_step`, steps per level-0 step, restart flag.
- `FireForcing` (each step): `time`, `dt`, x/y face velocities (before the dycore for `passive`/`lagged` coupling, after it for `synchronous`), cell-centre heights, and the requested optional fields.

## Zero-copy contract

- Forcing fields are references to ERF's own MultiFabs of the fire grid's level (its index space, box order, DistributionMapping and arena; device memory on GPU builds), valid only during `Advance()`.
- Fields a model reports through `FireDataPtr()` are its own MultiFabs, read in place; fluxes stay valid until the next `Advance()`.
- Fire box *n* refines atmospheric box *n* and has the same owner, so one `MFIter` serves both and no data moves between ranks.

## Inputs

- `erf.fire_model = None | <model>` selects the model (`FireType` in `ERF_DataStruct.H`).
- `erf.fire.*` are the coupling keys (grid ratio, anchor level, coupling type, heat injection, smoke); see `Docs/sphinx_doc/Inputs.rst`.
- `erf.fire.<model>.*` are each model's own keys.

## Adding a fire model

1. Create `Source/FireModel/<Model>/`, copy `ERF_FireModelTemplate.H` there, rename the class and implement the methods.
2. Read the model's inputs under `erf.fire.<model>.*`; use `Tools/` for wind and terrain.
3. Add the model to `FireType` and a `case` in `ERF::initializeFire` (`fire_model.SetModel<Model>()`) under the model's build define.
4. Add its sources to `CMake/BuildERFExe.cmake` and its `Make.package` to `Exec/Make.ERF`; an external library also needs its path.
5. Run the `FireModel` unit tests (`Tests/Unit/Fire/ERF_GTestFireModel.cpp`, whose stand-in is the smallest working model) with the new model.

No file in `Coupling/`, `Tools/` or `Features/` changes when a model is added.

The physics of the coupling (fire grid, coupling modes, heat profile, smoke)
is described in `Docs/sphinx_doc/theory/FireCoupling.rst`.
