# ELMFIRE Constant Wind Test Case

## Overview

This test case demonstrates ELMFIRE fire spread with constant wind (offline mode, no atmospheric coupling). The atmosphere is quiescent, and fire wind is specified via configuration parameters.

## Test Configuration

- **Domain**: 3 km × 3 km horizontal, 100 × 100 cells (30 m resolution)
- **Duration**: 3600 s (1 hour)
- **Fuel Model**: 102 (GR2 - moderate grass)
- **Wind**: 4.89 m/s eastward (constant)
- **Ignition**: Center of domain (1500 m, 1500 m), 50 m radius circle
- **Terrain**: Flat

## Running the Test

### Build

```bash
cd /path/to/ERF/mybuild
cmake .. -DERF_ENABLE_FIRE=ON
make -j8 erf_elmfire_constant
```

### Run

```bash
cd Exec/FireTest_ELMFIRE_Constant
./erf_elmfire_constant inputs_elmfire_constant
```

### MPI Run

```bash
mpirun -np 4 ./erf_elmfire_constant inputs_elmfire_constant
```

## Expected Output

- **Plot files**: `plt00000` through `plt00120` (every 30 seconds)
- **Fire area**: ~350 m² burned by end of simulation
- **Fire shape**: Elliptical, elongated in wind direction (eastward)

## Plot Variables

- `fire_phi`: Level set function (φ < 0 = burned)
- `fire_heat_flux`: Sensible heat flux from fire (W/m²)
- `fire_vel_u`: Fire spread velocity, x-component (m/s)
- `fire_vel_v`: Fire spread velocity, y-component (m/s)

## Visualization

### Using ParaView

1. Open ParaView
2. File → Open → Select all `plt*` directories
3. Click "Apply"
4. Select variable: `fire_phi`
5. Add filter: Contour (value = 0 for fire perimeter)
6. Play animation to see fire spread

### Using Python

```python
import yt

# Load dataset
ds = yt.load("plt00120/Header")

# Plot fire level set
slc = yt.SlicePlot(ds, 'z', 'fire_phi', origin='native')
slc.set_zlim('fire_phi', -1000, 1000)
slc.set_cmap('fire_phi', 'RdBu')
slc.save()
```

## Configuration Parameters

Key parameters in `inputs_elmfire_constant`:

```
# Fire model selection
erf.fire_model = "ELMFIRE"

# ELMFIRE configuration
fire.mode = "constant_wind"       # Offline mode
fire.constant_wind_u = 4.89       # m/s (eastward)
fire.constant_wind_v = 0.0        # m/s
fire.uniform_fuel_id = 102        # GR2 moderate grass
fire.ignition_x = 1500.0          # m
fire.ignition_y = 1500.0          # m
fire.ignition_radius = 50.0       # m
fire.ignition_time = 0.0          # s
```

## Verification

Expected fire spread characteristics:
- **Head fire ROS**: ~0.3-0.4 m/s (downwind)
- **Flank fire ROS**: ~0.15-0.20 m/s (crosswind)
- **Back fire ROS**: ~0.05-0.08 m/s (upwind)
- **Length-to-width ratio**: ~2.5:1

## Troubleshooting

### Fire doesn't ignite
- Check `fire.ignition_time` is within simulation time
- Check ignition location is inside domain
- Verify `fire.uniform_fuel_id` matches fuel database

### Simulation crashes
- Reduce `fire.cfl` (try 0.5 instead of 0.9)
- Increase `fire.band_thickness` (try 10 instead of 5)
- Check that ELMFIRE library was linked correctly

### No plot output
- Check `erf.plot_int_1` is reasonable (30 s recommended)
- Verify `erf.plot_file_1 = plt` is set
- Ensure `fire_phi` is in `erf.plot_vars_1`

## Notes

- This test case runs in **offline mode**: fire wind is constant and does not interact with atmosphere
- For **coupled mode** with fire-atmosphere interaction, see `FireTest_ELMFIRE_Coupled`
- Fuel moisture is constant (3% for dead fuels, 30-60% for live fuels)
- Terrain slope is flat (0°)

## References

- ELMFIRE user guide: See ELMFIRE-CPP_MPI/README.md
- Rothermel (1972): A mathematical model for predicting fire spread
- Anderson (2015): Aids to determining fuel models for estimating fire behavior
