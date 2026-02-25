#ifndef ERF_FIRE_LEVELSET_H
#define ERF_FIRE_LEVELSET_H

#include "ERF_NullFire.H"
#include <AMReX_TableData.H>

namespace FireVar {
    enum {
        phi = 0,        // Level set function (distance to fire front)
        fuel_frac,      // Fuel fraction (0=burned, 1=unburned)
        fuel_time,      // Time since ignition
        heat_flux,      // Heat flux to atmosphere
        NumVars
    };
}

/**
 * Level-set based fire spread model
 * Based on: Mell et al. (2007), Int. J. Wildland Fire
 */
class FireLevelSet : public NullFire {
public:
    FireLevelSet() : m_wind_u(nullptr), m_wind_v(nullptr) {}  // ← Initialize pointers
    ~FireLevelSet() override = default;

    void Define(SolverChoice& sc) override;
    
    void Init(int lev,
             const amrex::MultiFab& cons_in,
             const amrex::Geometry& geom,
             const amrex::Real& dt_advance) override;
    
    void Advance(const amrex::Real& dt,
                const amrex::MultiFab* U,
                const amrex::MultiFab* V,
                const amrex::MultiFab* terrain) override;

    void GetHeatFlux(amrex::MultiFab& heat_flux) override;
    void GetMoistureFlux(amrex::MultiFab& moisture_flux) override;
    
    bool is_active() const override { return m_active; }
    amrex::Real get_fire_area() const override;

    void GetPlotVarNames(amrex::Vector<std::string>& names) const override;
    
    void GetPlotVar(const std::string& name,
                   amrex::MultiFab& mf,
                   int lev) const override;

private:
    // Core computation methods
    void ComputeTerrainSlopes(const amrex::MultiFab* terrain, const amrex::Real* dx);
    void ComputeSpreadRate(const amrex::MultiFab& wind_speed,
                          const amrex::MultiFab& wind_dir);
    void AdvectLevelSet(const amrex::Real& dt);
    void ReinitializeLevelSet();
    void BurnFuel(const amrex::Real& dt);
    void ComputeHeatRelease();
    
    // Helper methods
    void SetIgnition(const amrex::Vector<amrex::Real>& location,
                    const amrex::Real& radius);
    
    // REMOVED TagRefinement - will add in Phase 1b with proper AMR support
    
    // Data members
    bool m_active = false;
    int m_lev = 0;
    amrex::Real m_dt = 0.0;
    
    // Fire state variables
    amrex::MultiFab m_fire_data;  // phi, fuel_frac, etc.
    
    // Fuel parameters (Rothermel model)
    amrex::Real m_fuel_density = 512.0;      // kg/m³ (ovendry)
    amrex::Real m_fuel_moisture = 0.10;      // fraction
    amrex::Real m_fuel_heat = 18.6e6;        // J/kg
    amrex::Real m_fuel_depth = 0.6;          // m
    amrex::Real m_fuel_savr = 3500.0;        // surface-to-volume ratio (1/m)
    
    // Spread parameters
    amrex::Real m_ros_base = 0.1;            // base rate of spread (m/s)
    amrex::Real m_wind_factor = 0.04;        // wind influence
    amrex::Real m_slope_factor = 0.2;        // slope influence
    
    // Refinement parameters
    amrex::Real m_refine_distance = 50.0;    // refine within 50m of fire front
    
    // Ignition parameters (from SolverChoice)
    amrex::Real m_ignition_time = 0.0;       // When to ignite (seconds)
    amrex::Real m_ignition_x = 0.0;          // Ignition x-coordinate (m)
    amrex::Real m_ignition_y = 0.0;          // Ignition y-coordinate (m)
    amrex::Real m_ignition_radius = 10.0;    // Ignition radius (m)

    // Wind source control
    std::string m_met_source = "Atmosphere";     // "Constant" or "Atmosphere"
    amrex::Real m_constant_u = 0.0;             // For idealized tests
    amrex::Real m_constant_v = 0.0;             // For idealized tests
    
    // Geometry
    const amrex::Geometry* m_geom = nullptr;

    // Wind field pointers
    const amrex::MultiFab* m_wind_u = nullptr;
    const amrex::MultiFab* m_wind_v = nullptr;
    
    // Terrain
    const amrex::MultiFab* m_terrain = nullptr;
};

#endif