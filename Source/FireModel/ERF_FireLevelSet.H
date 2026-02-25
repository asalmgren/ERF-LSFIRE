#ifndef ERF_FIRE_LEVELSET_H
#define ERF_FIRE_LEVELSET_H

#include "ERF_NullFire.H"
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <vector>
#include <map>
#include <string>
#include <memory>

#include "ERF_FireDataGatherer.H"  

namespace FireVar {
    enum {
        phi = 0,
        fuel_frac,
        fuel_time,
        heat_flux,
        burned_area,
        NumVars
    };
}

// ========================================
// FUEL MODEL DATA (from CSV database)
// ========================================
struct FuelModelData {
    int id;
    std::string name;
    bool has_live_fuels;
    amrex::Real w_1hr;
    amrex::Real w_10hr;
    amrex::Real w_100hr;
    amrex::Real w_lh;
    amrex::Real w_lw;
    amrex::Real sigma_1hr;
    amrex::Real sigma_10hr;
    amrex::Real sigma_100hr;
    amrex::Real depth;
    amrex::Real mex_dead;
    amrex::Real hoc;
    
    FuelModelData() : 
        id(0), name(""), has_live_fuels(false),
        w_1hr(0), w_10hr(0), w_100hr(0), w_lh(0), w_lw(0),
        sigma_1hr(2000), sigma_10hr(109), sigma_100hr(30),
        depth(1.0), mex_dead(0.15), hoc(8000) {}
};

// ========================================
// FUEL MODEL TABLE
// ========================================
struct FuelModelTable {
    amrex::Real wprimenumer[4];
    amrex::Real f[6];
    amrex::Real feps[6];
    amrex::Real fmex[4];
    amrex::Real mex_dead;
    amrex::Real mex_live;
    amrex::Real r_mprimedenome14sum_mex_dead;
    amrex::Real rhob;
    amrex::Real sigma_bar;
    amrex::Real beta;
    amrex::Real beta_opt;
    amrex::Real eta_s;
    amrex::Real xi;
    amrex::Real tr;
    amrex::Real gp_wnd_emd_es_hoc;
    amrex::Real gp_wnl_eml_es_hoc;
    amrex::Real f_dead;
    amrex::Real f_live;
    
    FuelModelTable() : 
        mex_dead(0.3), mex_live(1.5), 
        r_mprimedenome14sum_mex_dead(0),
        rhob(0), sigma_bar(2000.0), beta(0.001), beta_opt(0.02), eta_s(0.41),
        xi(0.4), tr(0.2),
        gp_wnd_emd_es_hoc(0), gp_wnl_eml_es_hoc(0),
        f_dead(0.5), f_live(0.5) {
        for (int i=0; i<4; i++) {
            wprimenumer[i] = 0;
            fmex[i] = 0;
        }
        for (int i=0; i<6; i++) {
            f[i] = 0;
            feps[i] = 0;
        }
    }
};

// ========================================
// FIRE NODE (Each rank's local fire cells)
// ========================================
struct FireNode {
    int i, j;
    amrex::Real time_added;
    amrex::Real phip_old;
    amrex::Real phip_rk1;  // ✅ Added for RK2 integration
    amrex::Real velocity;
    amrex::Real vs0;
    amrex::Real velocity_dms;
    amrex::Real velocity_dms_surface;
    amrex::Real vback;
    amrex::Real ir;
    amrex::Real hpua_surface;
    amrex::Real flin_dms_surface;
    amrex::Real dphidx, dphidy;  // ✅ Added raw gradients
    amrex::Real gradmag;          // ✅ Added gradient magnitude
    amrex::Real dphidx_limited;
    amrex::Real dphidy_limited;
    amrex::Real dphidx_upwind, dphidy_upwind;  // ✅ Added for upwind schemes
    amrex::Real normx, normy;
    amrex::Real normx_dms, normy_dms;
    amrex::Real ux, uy;
    amrex::Real low;
    amrex::Real b_coeff;
    amrex::Real phiwterm;
    amrex::Real phisterm;
    amrex::Real phisx, phisy;
    amrex::Real phiw_surface;
    amrex::Real phis_surface;
    amrex::Real uxousx, uyousy;
    amrex::Real ws20, wd20, u_wind, v_wind;
    amrex::Real m1, m10, m100;
    amrex::Real mlh, mlw;
    amrex::Real slope_x, slope_y;
    amrex::Real tanslp2;
    amrex::Real wsmf;
    amrex::Real adj;
    amrex::Real curvature;  // ✅ Added for advanced models
    int ifbfm;
    bool burned;
    bool front;
    
    FireNode() : i(0), j(0), time_added(0), phip_old(0), phip_rk1(0),
                 velocity(0), vs0(0), velocity_dms(0), velocity_dms_surface(0), vback(0),
                 ir(0), hpua_surface(0), flin_dms_surface(0),
                 dphidx(0), dphidy(0), gradmag(1.0),
                 dphidx_limited(0), dphidy_limited(0),
                 dphidx_upwind(0), dphidy_upwind(0),
                 normx(0), normy(1),
                 normx_dms(1), normy_dms(0),
                 ux(0), uy(0),
                 low(1.0),
                 b_coeff(2.0), phiwterm(1.0), phisterm(5.275),
                 phisx(0), phisy(0),
                 phiw_surface(0), phis_surface(0),
                 uxousx(1), uyousy(1),
                 ws20(0), wd20(0), u_wind(0), v_wind(0),
                 m1(0.06), m10(0.08), m100(0.10),
                 mlh(0.5), mlw(1.0),
                 slope_x(0), slope_y(0), tanslp2(0), 
                 wsmf(0), adj(1.0), curvature(0),
                 ifbfm(1), burned(false), front(false) {}
};

// ========================================
// MAIN CLASS: DISTRIBUTED FIRE MODEL
// ========================================
class FireLevelSet : public NullFire {
public:
    FireLevelSet() = default;
    ~FireLevelSet() override;
    // ============================================
    // PUBLIC INTERFACE (ERF integration)
    // ============================================
    void Define(SolverChoice& sc) override;
    void Init(int lev, const amrex::MultiFab& cons_in, 
              const amrex::Geometry& geom, const amrex::Real& dt_advance) override;
    void Advance(const amrex::Real& dt, const amrex::MultiFab* U, 
                 const amrex::MultiFab* V, const amrex::MultiFab* terrain,
                 const amrex::iMultiFab* luindex) override;
    void GetHeatFlux(amrex::MultiFab& heat_flux) override;
    void GetMoistureFlux(amrex::MultiFab& moisture_flux) override;
    amrex::Real get_fire_area() const override;
    void GetPlotVarNames(amrex::Vector<std::string>& names) const override;
    void GetPlotVar(const std::string& name, amrex::MultiFab& mf, int lev) const override;

private:
    // ============================================
    // FUEL DATABASE (Static, shared by all ranks)
    // ============================================
    static std::map<int, FuelModelData> m_fuel_database;
    static bool m_fuel_database_loaded;
    
    void LoadFuelDatabase();
    void UpdateFuelModelsFromLUINDEX(const amrex::iMultiFab* luindex);
    int MapLUINDEXtoFuelModel(int lu_index);    
    FuelModelTable GetFuelModelTable(int fuel_model_id, amrex::Real mlh);
    FuelModelTable GetHardcodedFuelModel(int fuel_model_id, amrex::Real mlh);
    int m_default_fuel_model = 101;
    
    // ============================================
    // REINITIALIZATION (Distributed with gather-broadcast)
    // ============================================
    void ReinitializeSignedDistanceDistributed();

    // ============================================
    // NARROW BAND MANAGEMENT (Each rank manages local cells)
    // ============================================
    void TagBandLocal(int i_center, int j_center, amrex::Real time);
    void UntagCellsLocal(amrex::Real current_time);
    void UpdateBurnedCellsLocal(amrex::Real current_time);
    void BurnFuelLocal(const amrex::Real& dt);
    
    // ============================================
    // EXTERNAL DATA (Wind, terrain, ignition)
    // ============================================
    void UpdateNodeWeatherDataLocal();
    void ComputeTerrainSlopes(const amrex::MultiFab* terrain, const amrex::Real* dx);
    void SetIgnition(const amrex::Vector<amrex::Real>& location, const amrex::Real& radius);
    
    // ============================================
    // OUTPUT (Each rank contributes)
    // ============================================
    void ComputeHeatReleaseLocal();
    bool IsNonburnable(int i, int j);
    
    // ============================================
    // PHYSICS (Rothermel model, each rank local)
    // ============================================
    void SurfaceSpreadRateLocal();
    void UxAndUyEllipticalLocal(int step, amrex::Real time);
    
    // ============================================
    // NUMERICS (Level set integration, local)
    // ============================================
    void CalcNormalVectorsLocal(int step);
    void ApplyFluxLimiterLocal();
    void Rk2IntegrateLocal(const amrex::Real& dt, int step);  
    amrex::Real ComputeAdaptiveTimestep(amrex::Real max_velocity);
    amrex::Real GetMaxFireVelocityLocal();
    
    // ✅ Additional numerics utilities
    void ComputeGradientsLocal();
    void ComputeCurvatureLocal();
    void ComputeUpwindGradientsLocal();
    void FastMarchingLocal();
    amrex::Real ComputeDistanceToPoint(int i, int j, int i_ref, int j_ref) const;
    void CheckGradientQuality();
    amrex::Real ComputeMaxCFL(amrex::Real dt) const;
    void CheckMassConservation();

    // ============================================
    // DATA MEMBERS
    // ============================================
    
    amrex::Real m_current_time = 0.0;
    int m_lev = 0;
    
    // ✅ DISTRIBUTED FIRE LISTS (Each rank has its own)
    std::vector<FireNode> m_local_tagged;   // Active fire cells owned by this rank
    std::vector<FireNode> m_local_burned;   // Burned cells owned by this rank
    std::vector<FireNode> m_local_front; // front cells owned by this rank

    // ✅ Data gatherer for wind/terrain communication
    std::unique_ptr<FireDataGatherer> m_data_gatherer;
    
    // ============================================
    // PARAMETERS (All ranks have same values)
    // ============================================
    int m_bandthickness = 10;
    amrex::Real m_target_cfl = 0.8;
    amrex::Real m_simulation_dt_max = 30.0;
    int m_untag_interval = 5;

    int m_verbosity = 2;
    
    amrex::Real m_ignition_time = 0.0;
    amrex::Real m_ignition_x = 0.0;
    amrex::Real m_ignition_y = 0.0;
    amrex::Real m_ignition_radius = 500.0;

    // ============================================
    // MULTIFAB DATA (Distributed with ghost cells)
    // ============================================
    amrex::MultiFab m_fire_data;  // φ, fuel_frac, heat_flux, etc. (with ghost cells!)
    std::unique_ptr<amrex::MultiFab> m_slope_x;
    std::unique_ptr<amrex::MultiFab> m_slope_y;
    bool m_needs_initial_tagging = false;

    // ============================================
    // WEATHER DATA
    // ============================================
    std::string m_met_source = "Constant";
    amrex::Real m_constant_u = 5.0;
    amrex::Real m_constant_v = 0.0;
    amrex::Real m_fuel_moisture = 0.06;
    
    // ============================================
    // GEOMETRY & EXTERNAL DATA POINTERS
    // ============================================
    const amrex::Geometry* m_geom = nullptr;
    const amrex::MultiFab* m_wind_u = nullptr;
    const amrex::MultiFab* m_wind_v = nullptr;
    const amrex::MultiFab* m_terrain = nullptr;
    const amrex::iMultiFab* m_luindex = nullptr;
};

#endif