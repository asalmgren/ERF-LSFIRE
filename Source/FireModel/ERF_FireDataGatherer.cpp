#include "ERF_FireDataGatherer.H"

using namespace amrex;

FireDataGatherer::FireDataGatherer(const Geometry* geom) 
    : m_geom(geom) 
{
    m_wind_map.clear();
    m_slope_map.clear();
    m_fuel_map.clear();
}

void FireDataGatherer::GatherTerrainStatic(const MultiFab* slope_x, 
                                           const MultiFab* slope_y) {
    if (!slope_x || !slope_y) {
        if (ParallelDescriptor::MyProc() == 0) {
            const Box& domain = m_geom->Domain();
            for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); j++) {
                for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); i++) {
                    m_slope_map[{i, j}] = {0.0, 0.0};
                }
            }
        }
        return;
    }
    
    BL_PROFILE("FireDataGatherer::GatherTerrainStatic");
    
#ifdef AMREX_USE_MPI
    int nprocs = ParallelDescriptor::NProcs();
    int myproc = ParallelDescriptor::MyProc();
    
    std::vector<Real> local_data;
    
    for (MFIter mfi(*slope_x); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> sx_arr = slope_x->const_array(mfi);
        Array4<Real const> sy_arr = slope_y->const_array(mfi);
        
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++) {
                local_data.push_back(static_cast<Real>(i));
                local_data.push_back(static_cast<Real>(j));
                local_data.push_back(sx_arr(i, j, 0));
                local_data.push_back(sy_arr(i, j, 0));
            }
        }
    }
    
    int local_count = local_data.size() / 4;
    
    std::vector<int> counts(nprocs, 0);
    MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 
               0, ParallelDescriptor::Communicator());
    
    std::vector<int> displs;
    int total_count = 0;
    if (myproc == 0) {
        displs.resize(nprocs);
        for (int p = 0; p < nprocs; p++) {
            displs[p] = total_count;
            total_count += counts[p];
        }
    }
    
    std::vector<Real> global_data;
    if (myproc == 0) {
        global_data.resize(total_count * 4);
    }
    
    std::vector<int> counts_bytes(nprocs), displs_bytes(nprocs);
    for (int p = 0; p < nprocs; p++) {
        counts_bytes[p] = counts[p] * 4;
        displs_bytes[p] = (myproc == 0) ? displs[p] * 4 : 0;
    }
    
    MPI_Gatherv(local_data.data(), local_count * 4, ParallelDescriptor::Mpi_typemap<Real>::type(),
                global_data.data(), counts_bytes.data(), displs_bytes.data(), 
                ParallelDescriptor::Mpi_typemap<Real>::type(), 0, ParallelDescriptor::Communicator());
    
    if (myproc == 0) {
        m_slope_map.clear();
        
        for (int k = 0; k < total_count; k++) {
            int i = static_cast<int>(global_data[4*k]);
            int j = static_cast<int>(global_data[4*k + 1]);
            Real sx = global_data[4*k + 2];
            Real sy = global_data[4*k + 3];
            
            m_slope_map[{i, j}] = {sx, sy};
        }
        
        if (ParallelDescriptor::IOProcessor()) {
            Print() << "🗻  Gathered " << m_slope_map.size() 
                    << " terrain slopes to Rank 0\n";
        }
    }
#else
    m_slope_map.clear();
    for (MFIter mfi(*slope_x); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> sx_arr = slope_x->const_array(mfi);
        Array4<Real const> sy_arr = slope_y->const_array(mfi);
        
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++) {
                m_slope_map[{i, j}] = {sx_arr(i, j, 0), sy_arr(i, j, 0)};
            }
        }
    }
#endif
    
    ParallelDescriptor::Barrier();
}

void FireDataGatherer::GatherFuelModelsStatic(const iMultiFab* luindex,
                                              int default_fuel_model) {
    if (!luindex) {
        if (ParallelDescriptor::MyProc() == 0) {
            const Box& domain = m_geom->Domain();
            for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); j++) {
                for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); i++) {
                    m_fuel_map[{i, j}] = default_fuel_model;
                }
            }
            
            if (ParallelDescriptor::IOProcessor()) {
                Print() << "🌲 Using default fuel model " << default_fuel_model 
                        << " everywhere\n";
            }
        }
        return;
    }
    
    BL_PROFILE("FireDataGatherer::GatherFuelModelsStatic");
    
#ifdef AMREX_USE_MPI
    int nprocs = ParallelDescriptor::NProcs();
    int myproc = ParallelDescriptor::MyProc();
    
    std::vector<int> local_data;
    
    for (MFIter mfi(*luindex); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<int const> lu_arr = luindex->const_array(mfi);
        
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++) {
                int lu_index = lu_arr(i, j, 0);
                int fuel_id = (lu_index > 0) ? lu_index : default_fuel_model;
                
                local_data.push_back(i);
                local_data.push_back(j);
                local_data.push_back(fuel_id);
            }
        }
    }
    
    int local_count = local_data.size() / 3;
    
    std::vector<int> counts(nprocs, 0);
    MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 
               0, ParallelDescriptor::Communicator());
    
    std::vector<int> displs;
    int total_count = 0;
    if (myproc == 0) {
        displs.resize(nprocs);
        for (int p = 0; p < nprocs; p++) {
            displs[p] = total_count;
            total_count += counts[p];
        }
    }
    
    std::vector<int> global_data;
    if (myproc == 0) {
        global_data.resize(total_count * 3);
    }
    
    std::vector<int> counts_bytes(nprocs), displs_bytes(nprocs);
    for (int p = 0; p < nprocs; p++) {
        counts_bytes[p] = counts[p] * 3;
        displs_bytes[p] = (myproc == 0) ? displs[p] * 3 : 0;
    }
    
    MPI_Gatherv(local_data.data(), local_count * 3, MPI_INT,
                global_data.data(), counts_bytes.data(), displs_bytes.data(), 
                MPI_INT, 0, ParallelDescriptor::Communicator());
    
    if (myproc == 0) {
        m_fuel_map.clear();
        
        for (int k = 0; k < total_count; k++) {
            int i = global_data[3*k];
            int j = global_data[3*k + 1];
            int fuel_id = global_data[3*k + 2];
            
            m_fuel_map[{i, j}] = fuel_id;
        }
        
        if (ParallelDescriptor::IOProcessor()) {
            Print() << "🌲 Gathered " << m_fuel_map.size() 
                    << " fuel model assignments to Rank 0\n";
        }
    }
#else
    m_fuel_map.clear();
    for (MFIter mfi(*luindex); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<int const> lu_arr = luindex->const_array(mfi);
        
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++) {
                int lu_index = lu_arr(i, j, 0);
                int fuel_id = (lu_index > 0) ? lu_index : default_fuel_model;
                m_fuel_map[{i, j}] = fuel_id;
            }
        }
    }
#endif
    
    ParallelDescriptor::Barrier();
}

std::vector<std::pair<int,int>> FireDataGatherer::BuildGatherList(
    const std::vector<std::pair<int,int>>& fire_cells,
    int search_radius) const {
    
    const Box& domain = m_geom->Domain();
    std::set<std::pair<int,int>> cells_to_gather;
    
    for (const auto& fc : fire_cells) {
        for (int j = fc.second - search_radius; j <= fc.second + search_radius; j++) {
            for (int i = fc.first - search_radius; i <= fc.first + search_radius; i++) {
                if (domain.contains(IntVect(i, j, 0))) {
                    cells_to_gather.insert({i, j});
                }
            }
        }
    }
    
    return std::vector<std::pair<int,int>>(cells_to_gather.begin(), cells_to_gather.end());
}

void FireDataGatherer::GatherWindsDynamic(const MultiFab* U, const MultiFab* V,
                                          const std::vector<std::pair<int,int>>& fire_cells,
                                          int search_radius) {
    if (!U || !V || m_use_constant_wind) return;
    
    BL_PROFILE("FireDataGatherer::GatherWindsDynamic");
    
    auto cells_to_gather = BuildGatherList(fire_cells, search_radius);
    
#ifdef AMREX_USE_MPI
    int nprocs = ParallelDescriptor::NProcs();
    int myproc = ParallelDescriptor::MyProc();
    
    std::vector<Real> local_data;
    
    for (MFIter mfi(*U); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> u_arr = U->const_array(mfi);
        Array4<Real const> v_arr = V->const_array(mfi);
        
        for (const auto& ij : cells_to_gather) {
            int i = ij.first;
            int j = ij.second;
            
            if (bx.contains(IntVect(i, j, 0))) {
                local_data.push_back(static_cast<Real>(i));
                local_data.push_back(static_cast<Real>(j));
                local_data.push_back(u_arr(i, j, 0));
                local_data.push_back(v_arr(i, j, 0));
            }
        }
    }
    
    int local_count = local_data.size() / 4;
    
    std::vector<int> counts(nprocs, 0);
    MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 
               0, ParallelDescriptor::Communicator());
    
    std::vector<int> displs;
    int total_count = 0;
    if (myproc == 0) {
        displs.resize(nprocs);
        for (int p = 0; p < nprocs; p++) {
            displs[p] = total_count;
            total_count += counts[p];
        }
    }
    
    std::vector<Real> global_data;
    if (myproc == 0) {
        global_data.resize(total_count * 4);
    }
    
    std::vector<int> counts_bytes(nprocs), displs_bytes(nprocs);
    for (int p = 0; p < nprocs; p++) {
        counts_bytes[p] = counts[p] * 4;
        displs_bytes[p] = (myproc == 0) ? displs[p] * 4 : 0;
    }
    
    MPI_Gatherv(local_data.data(), local_count * 4, ParallelDescriptor::Mpi_typemap<Real>::type(),
                global_data.data(), counts_bytes.data(), displs_bytes.data(), 
                ParallelDescriptor::Mpi_typemap<Real>::type(), 0, ParallelDescriptor::Communicator());
    
    if (myproc == 0) {
        m_wind_map.clear();
        
        for (int k = 0; k < total_count; k++) {
            int i = static_cast<int>(global_data[4*k]);
            int j = static_cast<int>(global_data[4*k + 1]);
            Real u = global_data[4*k + 2];
            Real v = global_data[4*k + 3];
            
            m_wind_map[{i, j}] = {u, v};
        }
    }
#else
    m_wind_map.clear();
    for (MFIter mfi(*U); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> u_arr = U->const_array(mfi);
        Array4<Real const> v_arr = V->const_array(mfi);
        
        for (const auto& ij : cells_to_gather) {
            int i = ij.first;
            int j = ij.second;
            
            if (bx.contains(IntVect(i, j, 0))) {
                m_wind_map[{i, j}] = {u_arr(i, j, 0), v_arr(i, j, 0)};
            }
        }
    }
#endif
}

std::pair<Real, Real> FireDataGatherer::GetWind(int i, int j) const {
    if (m_use_constant_wind) {
        return {m_constant_u, m_constant_v};
    }
    
    auto it = m_wind_map.find({i, j});
    if (it != m_wind_map.end()) {
        return it->second;
    }
    return {m_constant_u, m_constant_v};
}

std::pair<Real, Real> FireDataGatherer::GetSlope(int i, int j) const {
    auto it = m_slope_map.find({i, j});
    if (it != m_slope_map.end()) {
        return it->second;
    }
    return {0.0, 0.0};
}

int FireDataGatherer::GetFuelModel(int i, int j) const {
    auto it = m_fuel_map.find({i, j});
    if (it != m_fuel_map.end()) {
        return it->second;
    }
    return 101;
}

bool FireDataGatherer::HasData(int i, int j) const {
    return m_wind_map.find({i, j}) != m_wind_map.end();
}

void FireDataGatherer::SetConstantWind(Real u, Real v) {
    m_constant_u = u;
    m_constant_v = v;
    m_use_constant_wind = true;
}

void FireDataGatherer::PrintStatistics() const {
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "\n=== FireDataGatherer Statistics ===\n";
        Print() << "  Wind data points: " << m_wind_map.size() << "\n";
        Print() << "  Slope data points: " << m_slope_map.size() << "\n";
        Print() << "  Fuel model assignments: " << m_fuel_map.size() << "\n";
        Print() << "====================================\n\n";
    }
}