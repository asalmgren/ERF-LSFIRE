#ifndef ERF_FIRE_DATA_GATHERER_H
#define ERF_FIRE_DATA_GATHERER_H

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <map>
#include <vector>
#include <set>

/**
 * @brief Helper class to centralize MPI communication for fire module
 * 
 * Handles gathering atmospheric/terrain data to Rank 0 and broadcasting results
 */
class FireDataGatherer {
public:
    FireDataGatherer(const amrex::Geometry* geom);
    ~FireDataGatherer() = default;
    
    // Static data gathering (once during Init)
    void GatherTerrainStatic(const amrex::MultiFab* slope_x, 
                             const amrex::MultiFab* slope_y);
    void GatherFuelModelsStatic(const amrex::iMultiFab* luindex, 
                                int default_fuel_model);
    
    // Dynamic data gathering (every timestep)
    void GatherWindsDynamic(const amrex::MultiFab* U, 
                           const amrex::MultiFab* V,
                           const std::vector<std::pair<int,int>>& fire_cells,
                           int search_radius = 10);
    
    void SetConstantWind(amrex::Real u, amrex::Real v);
    
    // Data access (Rank 0 only)
    std::pair<amrex::Real, amrex::Real> GetWind(int i, int j) const;
    std::pair<amrex::Real, amrex::Real> GetSlope(int i, int j) const;
    int GetFuelModel(int i, int j) const;
    bool HasData(int i, int j) const;
    
    void PrintStatistics() const;

private:
    const amrex::Geometry* m_geom;
    
    std::map<std::pair<int,int>, std::pair<amrex::Real, amrex::Real>> m_wind_map;
    std::map<std::pair<int,int>, std::pair<amrex::Real, amrex::Real>> m_slope_map;
    std::map<std::pair<int,int>, int> m_fuel_map;
    
    amrex::Real m_constant_u = 0.0;
    amrex::Real m_constant_v = 0.0;
    bool m_use_constant_wind = false;
    
    std::vector<std::pair<int,int>> BuildGatherList(
        const std::vector<std::pair<int,int>>& fire_cells,
        int search_radius) const;
};

#endif