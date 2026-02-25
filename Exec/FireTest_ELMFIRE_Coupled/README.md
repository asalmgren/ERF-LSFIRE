# ELMFIRE Atmospheric Coupling Test Case

## Overview

This test case demonstrates ELMFIRE fire spread with atmospheric coupling (two-way interaction between fire and atmosphere). Fire extracts wind from ERF atmosphere and returns heat flux to drive atmospheric circulation.

**NOTE**: This test case requires Phase 4 implementation (atmospheric coupling in ELMFIRE library). It will abort with "not yet implemented" message until Phase 4 is complete.

## Test Configuration

- **Domain**: 3 km × 3 km horizontal, 2 km vertical
- **Resolution**: 100 × 100 × 64 cells (30 m horizontal, ~31 m vertical)
- **Duration**: 3600 s (1 hour)
- **Fuel Model**: 102 (GR2 - moderate grass)
- **Initial Wind**: 5 m/s eastward (uniform)
- **Ignition**: Center of domain (1500 m, 1500 m), 50 m radius circle, at t = 600 s
- **Terrain**: Flat
- **Turbulence**: Smagorinsky LES (Cs = 0.17)

## Running the Test

### Build

```bash
cd /path/to/ERF/mybuild
cmake .. -DERF_ENABLE_FIRE=ON
make -j8 erf_elmfire_coupled
```

### Run

```bash
cd Exec/FireTest_ELMFIRE_Coupled
./erf_elmfire_coupled inputs_elmfire_coupled
```

### MPI Run

```bash
mpirun -np 4 ./erf_elmfire_coupled inputs_elmfire_coupled
```

## Expected Output (after Phase 4 implementation)

- **Plot files**: `plt00000` through `plt00120` (every 30 seconds)
- **Fire area**: ~400 m² burned by end of simulation (larger than constant wind due to feedback)
- **Fire shape**: Elliptical with dynamic wind interaction
- **Atmospheric response**: Convective plume above fire, wind acceleration near fire front

## Plot Variables

### Fire Variables
- `fire_phi`: Level set function (φ < 0 = burned)
- `fire_heat_flux`: Sensible heat flux from fire (W/m²)
- `fire_vel_u`: Fire spread velocity, x-component (m/s)
- `fire_vel_v`: Fire spread velocity, y-component (m/s)

### Atmospheric Variables
- `density`: Air density (kg/m³)
- `x_velocity`, `y_velocity`, `z_velocity`: Wind components (m/s)
- `temp`: Temperature (K)
- `theta`: Potential temperature (K)

## Atmospheric Coupling Mechanism

### Phase 4 Implementation Details

1. **Wind Extraction** (Fire → Atmosphere):
   - ELMFIRE extracts surface wind (`U`, `V` at `k=0`) from ERF MultiFabs
   - Wind field is spatially varying and time-dependent
   - Fire spread responds to atmospheric wind

2. **Heat Flux Application** (Atmosphere → Fire):
   - ELMFIRE computes sensible heat flux from burning rate
   - Typical values: 10⁵ - 10⁶ W/m² in active fire regions
   - ERF applies heat flux as surface boundary condition
   - Drives convective plume and atmospheric circulation

3. **Feedback Loop**:
   - Fire-induced convection modifies surface wind
   - Modified wind affects fire spread
   - Positive feedback can accelerate fire spread

## Configuration Parameters

Key parameters in `inputs_elmfire_coupled`:

```
# Fire model selection
erf.fire_model = "ELMFIRE"

# ELMFIRE configuration
fire.mode = "atmospheric"         # Coupled mode
fire.uniform_fuel_id = 102        # GR2 moderate grass
fire.ignition_x = 1500.0          # m
fire.ignition_y = 1500.0          # m
fire.ignition_radius = 50.0       # m
fire.ignition_time = 600.0        # s (after atmosphere spins up)

# Turbulence (important for coupled mode)
erf.les_type = "Smagorinsky"
erf.Cs = 0.17
```

## Physics Comparison: Constant vs Coupled

| Aspect | Constant Wind | Atmospheric Coupling |
|--------|--------------|---------------------|
| Fire wind | Fixed 4.89 m/s | Variable, from ERF |
| Heat flux | Not applied to atmosphere | Applied as BC |
| Convective plume | No | Yes |
| Wind feedback | No | Yes |
| Fire spread rate | ~0.3 m/s (head) | ~0.35 m/s (faster due to feedback) |
| Simulation cost | Lower | Higher (atmosphere + turbulence) |

## Verification (after Phase 4)

Expected fire-atmosphere characteristics:
- **Convective plume**: Vertical velocity > 5 m/s above fire
- **Temperature rise**: +10-30 K in plume
- **Wind acceleration**: 20-40% increase near fire front
- **Fire spread enhancement**: 10-20% faster than constant wind case

## Troubleshooting

### "not yet implemented" error
- This is expected until Phase 4 is complete
- Use `FireTest_ELMFIRE_Constant` for testing until then

### Simulation crashes (after Phase 4)
- Reduce `erf.fixed_dt` (try 0.5 s)
- Reduce `fire.cfl` (try 0.5)
- Increase `fire.band_thickness` (try 10)
- Check CFL condition for atmospheric solver

### Fire doesn't respond to atmosphere
- Verify `fire.mode = "atmospheric"`
- Check that ERF is passing wind MultiFabs to fire model
- Verify surface wind extraction is working (add debug output)

## Development Status

- **Phase 1**: ✅ ELMFIRE library structure
- **Phase 2**: ✅ ERF integration (FireELMFIRE class)
- **Phase 3**: ✅ Test case setup (this directory)
- **Phase 4**: ⏳ Atmospheric coupling implementation (in progress)

## Phase 4 Implementation Checklist

- [ ] Extract surface wind from ERF MultiFabs (`U`, `V` at `k=0`)
- [ ] Update ELMFIRE node weather data with spatially-varying wind
- [ ] Calculate sensible heat flux from burning rate
- [ ] Return heat flux MultiFab to ERF
- [ ] Apply heat flux as surface boundary condition in ERF
- [ ] Test and validate fire-atmosphere feedback

## Notes

- Ignition delayed to t = 600 s to allow atmosphere to spin up turbulence
- Checkpoint files enabled (every 10 minutes) for long simulations
- LES turbulence recommended for realistic atmospheric response
- Periodic lateral boundaries to avoid edge effects

## References

- Clark et al. (1996): A coupled atmosphere-fire model
- Coen (2013): WRF-Fire: Coupled atmosphere-fire model
- Mandel et al. (2011): Coupled atmosphere-fire modeling with WRF-Fire
- Phase 4 implementation plan: See coupling plan documentation
