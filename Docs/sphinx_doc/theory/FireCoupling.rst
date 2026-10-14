.. role:: cpp(code)
   :language: c++

.. _sec:FireCoupling:

Fire-Atmosphere Coupling
========================

A wildland fire and the atmosphere drive each other: the wind spreads the
fire, and the fire's heat drives the wind. ERF provides the atmosphere side of
that coupling: where the fire grid sits, what the fire is handed each step,
and how its heat, water vapour and smoke enter the atmosphere. The fire
itself, the model that spreads the front and computes the heat it releases,
is a *fire model* that plugs into this coupling, the way land surface and
wind farm models plug into ERF. No fire model is built into this repository;
the coupling is the common part every model shares.

Building with the fire coupling
-------------------------------

The fire model container and its coupling are compiled with the CMake option
``ERF_ENABLE_FIRE=ON`` (``USE_FIRE = TRUE`` with GNU make); it is off by
default. With it on, the state carries one more passive scalar, the smoke
tracer (``smoke`` in the plotfile), whether or not a fire runs. Every
checkpoint, from any build, records the number of passive scalars in a small
file ``ScalarLayout``, so that a build configured differently refuses to
restart from it instead of reading the components shifted; a checkpoint
without the file is read as one scalar, so older checkpoints restart in builds
without the fire coupling.

:cpp:`erf.fire_model` selects the fire model (``None`` by default); a model
other than ``None`` must be built into the executable, and a build without
``ERF_ENABLE_FIRE`` stops on a fire deck. The keys that selected the fire
before :cpp:`erf.fire_model`, :cpp:`erf.fire.enable` and
:cpp:`erf.fire.engine`, stop the run with a message pointing to it.

The fire model container and fire models
----------------------------------------

The code follows the pattern of ERF's land surface and wind farm models, so
that more than one fire model can be coupled to ERF the same way
(``Source/FireModel``):

* **The container** (``FireModel`` in ``ERF_FireModel.H``, held by ERF as it
  holds its ``LandSurface``) keeps one model per level: a ``NullFire``
  (``Null/ERF_NullFire.H``) on every level, and the selected model on the
  level the fire grid refines, :cpp:`erf.fire.anchor_level`. It also holds
  everything on the atmosphere side (``Coupling/``): it checks that the mesh
  can carry a fire grid, builds that grid and hands the model the atmospheric
  state at the start of each of its steps. With ``"lagged"`` or
  ``"synchronous"`` coupling it averages the heat the model reports onto the
  atmospheric columns, holds it for one step, and injects it at every
  Runge-Kutta stage of the next step; with :cpp:`erf.fire.smoke_enable` it
  injects the smoke tracer the same way (``Features/``), also with
  ``"passive"`` coupling, which injects no heat. It saves the heat it holds in
  the checkpoint, so a restarted run injects the same heat in its first step as
  an uninterrupted one.
* **The fire model** (derived from ``FireModelBase`` in
  ``ERF_FireModelBase.H``) is the fire physics: it advances the front on the
  fire grid and reports, after each step, the sensible and latent heat flux
  into the air on the fire grid [W/m²], averaged over its step. The model
  reports its fields by name (``FireVar`` in ``ERF_FireModelBase.H``: the two
  fluxes, a level set, and the open fraction and roof height of
  :cpp:`erf.fire.heat_open_fraction`); a field it does not have is simply
  absent. It writes its own fire plotfile and, when it supports restart, its
  own checkpoint state; without restart support a restart with the fire
  running stops before reading the checkpoint.
* **Tools** (``Tools/``) are helpers a model may use: a wind sampler
  (``ERF_FireWindExtract.H``), the ground height of each atmospheric column
  (``ERF_FireTerrain.H``) and the optional atmospheric inputs
  (``ERF_FireAtmState.H``).

The fire grid refines the columns of its level by :cpp:`erf.fire.grid_ratio`
in x and y. Each fire box is an atmospheric box refined by the ratio, owned
by the same rank, and every x and y box length on that level must be
divisible by the ratio; the run stops at start-up naming the first box that
is not, with the inputs to change. The level must keep its boxes for the whole
run (on a level above 0, :cpp:`erf.regrid_int` must be off), its boxes must
form one rectangle, no box may be split in z, and the bottom boundary must be
a surface layer.

The model reads, at the start of each step, the face velocities of the fire
grid's level (before the dynamical core with ``"passive"`` and ``"lagged"``
coupling, after it with ``"synchronous"``) and the cell-centre heights. The
container computes the optional inputs only for a model that requests them
(``FireModelBase::RequestedInputs``): the vertical velocity, the air
temperature and relative humidity of the lowest cell before the dynamical
core, and the surface precipitation of every column from the microphysics
accumulators (the run then needs a microphysics scheme with precipitation).
These fields are ERF's own, read in place (no copy is made), and valid only
during the model's step. What the model does with them, including the height
it samples the wind at and any correction it applies to that wind, is the
model's choice.

Adding a fire model
~~~~~~~~~~~~~~~~~~~

``Source/FireModel/ERF_FireModelTemplate.H`` is the starting point for a new
fire model. It declares every method a fire model provides, each with a
comment on what ERF hands it, what ERF expects back and when it is called. It
holds declarations only: nothing in it has a body, no input selects it and
nothing creates one, so it adds no code to a run. ``ERF_FireModel.cpp``
checks at compile time that it declares every method ``FireModelBase``
requires (the pure virtual ones), and ``override`` checks the signature of
every method it declares; the optional methods are kept in it by hand.

To add a model:

#. Create ``Source/FireModel/<Model>/``, copy the template there, rename the
   class, and write the method bodies in a ``.cpp`` file there. Delete the
   optional methods the model does not provide; their defaults then apply.
#. Read the model's inputs under its own prefix (``erf.fire.<model>.*``), so
   they never mix with the coupling's keys or another model's.
#. Add the model to ``FireType`` (``ERF_DataStruct.H``) and install it in
   ``ERF::initializeFire`` with ``fire_model.SetModel<Model>()`` under the
   model's build define. Add the folder's sources to
   ``CMake/BuildERFExe.cmake`` (inside ``if(ERF_ENABLE_FIRE)``) and its
   ``Make.package`` to ``Exec/Make.ERF`` under ``USE_FIRE``.
#. Run the ``FireModel`` unit tests (``Tests/Unit/Fire/ERF_GTestFireModel.cpp``,
   whose stand-in model is the smallest model that works) with the new
   model, and a restart of a coupled run against an uninterrupted one: the
   atmosphere and fire plotfiles must agree at zero tolerance.

The model owns its fire physics and its own state, including its checkpoint
and plotfile; the container does the rest, as described on this page.

Coupling modes
--------------

:cpp:`erf.fire.coupling_type` selects which wind the fire sees and whether its
heat returns to the atmosphere:

- ``"passive"``: the fire spreads on the wind from the start of the step and
  its heat flux is never injected. Use it for fast screening runs and for
  tests where the atmosphere should not respond.
- ``"lagged"`` (default): the fire spreads on the wind from the start of the
  step; the flux it produces is stored and injected during the next step's
  Runge-Kutta stages. This one-step explicit lag is the WRF-SFIRE approach
  (Mandel et al., 2011).
- ``"synchronous"``: the fire spreads on the wind *after* the dynamical core
  has advanced, so the spread at step :math:`n+1` already reflects the
  momentum response to the heat injected from step :math:`n`; the flux still
  enters at the next step.

The fire runs once per step of its level in all three modes.

Injection into the atmosphere
-----------------------------

The fire-grid fluxes are area-averaged onto the atmospheric columns of the
fire grid's level and distributed vertically with an exponential profile,

.. math::

   H(z) = \frac{Q}{c_p}\, e^{-z/\alpha_g}, \qquad
   \frac{\partial (\rho\theta)}{\partial t} = -\frac{\partial H}{\partial z},

where :math:`z` is height above the first cell centre of the column (the
profile's datum, half a cell above the terrain, so the lowest cell takes the
share :math:`1 - e^{-\Delta z_1/\alpha_g}` of the flux, with
:math:`\Delta z_1` the distance from the first cell centre to the second, and
nothing is lost below it) and
:math:`\alpha_g` = :cpp:`erf.fire.heat_flux_alfg` (default 45 m) is the
e-folding depth: 37% of the flux remains at 45 m and less than 1% at 225 m.
Each cell's tendency is the flux difference between the profile at its own
centre and at the next one, divided by the cell's own thickness (face to
face, from the nodal heights on a terrain-fitted or stretched grid; only a
degenerate cell of zero thickness is floored), so the
column receives :math:`(1 - e^{-z_{top}/\alpha_g})\,Q` on any vertical
spacing, where :math:`z_{top}` is the height of the top cell centre above the
first one plus the mean vertical spacing of the level (the profile closes one
mean spacing above the top centre). On a uniform grid, heights above the first
centre equal the heights of the cell faces above the ground, so each cell takes
the profile between its own faces; on a grid stretched by a constant ratio
:math:`r` the weights (apart from the top cell's) are those of the face-based profile with the e-folding
depth :math:`2\alpha_g/(1+r)`, which places slightly more of the heat in the
lowest cells (about 1 % of the flux moves down at :math:`r = 1.05`). The latent
flux is distributed the same way into the vapour equation (divided by the
latent heat of vaporisation) when :cpp:`erf.fire.inject_latent` is true and a
moisture model is active.

The tendency above is a :math:`\rho\theta` budget: heating :math:`q` per unit
volume gives :math:`\partial\theta/\partial t = q / (\rho c_p \Pi)` with the
Exner function :math:`\Pi = (p/p_0)^{R_d/c_p}`, so the enthalpy
:math:`c_p \Pi\, \partial(\rho\theta)/\partial t` the column receives is
:math:`\Pi` times the flux, 4.5 % short at 850 hPa.
:cpp:`erf.fire.heat_tendency_exner = true` divides the tendency by
:math:`\Pi` of each cell so that the enthalpy integrates to the flux; the
default omits the factor, as WRF-SFIRE does. :cpp:`erf.fire.heat_tendency_density = true`
multiplies the tendency by the air density, so the column receives
:math:`\rho` times the flux (kept for comparison with runs made in that
form); set together, the two factors stack (:math:`\rho / \Pi`). With
:cpp:`erf.fire.fire_debug` the coupling prints the column-integrated heating
against the flux handed in, whose ratio is :math:`1 - e^{-z_{top}/\alpha_g}`
with neither factor (the factors multiply it by the profile-weighted mean of
:math:`\rho`, of :math:`1/\Pi`, or of :math:`\rho/\Pi` with both, essentially
their values near the surface).

:cpp:`erf.fire.fire_atm_feedback` multiplies both fluxes before injection;
zero gives one-way coupling with the fire still responding to the wind, one
gives the full feedback. At zero the coupling writes nothing into the source,
so with ``source_mode = "overwrite"`` the other theta sources survive there,
while any positive value replaces them.

The cell source is rebuilt by the atmosphere at every Runge-Kutta stage and
the fire tendency is applied after that rebuild. :cpp:`erf.fire.source_mode`
selects how: ``"overwrite"`` (default) replaces the potential-temperature and
vapour slots, which keeps the injected energy independent of the number of
stages but discards any other source already in those slots, such as the
immersed-forcing scalar relaxation, radiative heating or Rayleigh damping;
``"add"`` accumulates the fire tendency into the freshly rebuilt source, once
per stage, so nothing is lost and nothing is double counted. Use ``"add"``
whenever another source writes into those slots: immersed forcing applied on
the slow step (``erf.immersed_forcing_substep = false``, the anelastic
default), radiative heating, Rayleigh damping. With immersed forcing on the
acoustic substeps (the compressible default) the scalar relaxation never
meets the slow-step source and the two modes are identical, as they are with
no other source. Because the default would silently drop such a source, a run
that leaves :cpp:`erf.fire.source_mode` unset while the fire injects heat and
the atmosphere has radiation, :cpp:`erf.four_stream_radiation`,
:cpp:`erf.rayleigh_damp_T`, custom theta forcing or subsidence (or custom
moisture forcing with latent heat on), or slow-step immersed forcing stops at
start-up and names them; setting the key to either value is the explicit
choice.

The heat held between steps carries five ghost columns, holding the
neighbouring box's value (or the edge column's outside a non-periodic face),
for a reader that samples a halo around each tile.

Levels
~~~~~~

The sensible and latent tendencies and the smoke source go into the source of
the level the fire grid refines (:cpp:`erf.fire.anchor_level`, the finest
level by default) at every Runge-Kutta stage of that level's steps. Outside
the refined region, which has no fire grid, the held heat is zero. Coarser
levels receive the heat through average-down (``erf.coupling_type =
"TwoWay"``), which keeps the volume integral, so the heat on level 0 is the
heat injected on the fire's level. With the fire grid on a coarser level than
the finest, the same average-down replaces the heated cells under the finer
level with unheated ones and that heat is lost, and with ``"OneWay"`` the
finer level never sees it; the coupling warns at start-up whenever the level
is not the finest.

A restart that re-makes the level-0 grids (with
:cpp:`erf.regrid_level_0_on_restart`, or on more ranks than the checkpoint has
level-0 boxes) is not supported with a fire running and stops before
reading the checkpoint.

Heat placement around buildings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A model that knows where buildings stand can report, per atmospheric column,
the open fraction :math:`f` of its footprint and the mean roof height :math:`H`
of the rest, above the ground. With :cpp:`erf.fire.heat_open_fraction` the
coupling places the heat accordingly: a cell's share of the exponential profile
is scaled by :math:`f` where the cell lies below the roof, by one above it, and
linearly across the cell that holds the roof (the cell's faces measured from
the ground), and the column is renormalised over the shares that remain,

.. math::

   w_k = \frac{\bigl(e^{-z_k/\alpha_g} - e^{-z_{k+1}/\alpha_g}\bigr)\, s_k}
              {\sum_j \bigl(e^{-z_j/\alpha_g} - e^{-z_{j+1}/\alpha_g}\bigr)\, s_j}
         \,\bigl(1 - e^{-z_{top}/\alpha_g}\bigr),

so the energy injected per column is what the plain profile injects (exactly,
with neither tendency factor; the factors weigh the lifted heat by the
:math:`\rho` or :math:`1/\Pi` of the cells it moves to). With :math:`B` the
plain profile's share in the cells below the roof and
:math:`L = 1 - e^{-z_{top}/\alpha_g}`, the share placed below the roof falls
from :math:`B` to :math:`f B L / (L - (1 - f) B)`, more than :math:`f B`
because the renormalisation also raises the cells below the roof; the rest is
lifted above it. A model
that reports no open fraction leaves the plain profile in place, with a
warning at start-up.

Smoke tracer
------------

With :cpp:`erf.fire.smoke_enable` (which needs :cpp:`erf.transport_scalar`)
the smoke tracer is fed at the surface with

.. math::

   F_{smoke} = \epsilon\, \frac{Q}{H_c} \quad [\mathrm{kg\,m^{-2}\,s^{-1}}],

where :math:`\epsilon` = :cpp:`erf.fire.smoke_emission_factor` (kg of smoke
per kg of fuel, default 0.02) and :math:`H_c` is the heat the model releases
per kilogram of fuel burnt (``FireModelBase::SmokeHeatPerKg``). The emission
uses the same held, area-averaged heat flux as the energy injection and is
placed in the lowest atmospheric cell. The tracer is advected and diffused
like any other scalar and is plotted as ``smoke``; its integral is printed
with the other integrated quantities. All passive scalars share one boundary
condition, so on a face that ingests the ``scalar`` boundary plane
(:cpp:`erf.bndry_input_var_names`) the smoke takes the base scalar's plane
values too.

Inputs
------

The coupling's inputs are listed under :ref:`Inputs <sec:InputsFire>`; an
model's own inputs are its own.

References
----------

- Mandel, J., Beezley, J. D. and Kochanski, A. K. (2011). Coupled atmosphere-wildland fire modeling with WRF 3.3 and SFIRE 2011. Geoscientific Model Development, 4, 591-610.
