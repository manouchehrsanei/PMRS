//
//  TPMRSSimulationData.h
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#ifndef TPMRSSimulationData_h
#define TPMRSSimulationData_h

#include <stdio.h>
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include "tinystr.h"
#include "tinyxml.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <tuple>
#include "pzerror.h"
#include "TPZGmshReader.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckgeom.h"
#include "TPMRSUndrainedParameters.h"
#include "TPMRSPoroMechParameters.h"
#include "TPMRSPhiParameters.h"
#include "TPMRSKappaParameters.h"
#include "TPMRSPlasticityParameters.h"


class TPMRSSimulationData
{
    
protected:
    
    /// Brief Time step size
    REAL m_dt;
    
    /// Brief Number of time steps
    int m_n_steps;
    
    /// Brief Store time values to be reported
    TPZStack< REAL , 500 > m_reporting_times;
    
    /// Brief Current time value
    REAL m_time;
    
    /// Brief Number of iteration
    int m_n_iteraions;
    
    //// Brief Residue overal tolerance
    REAL m_epsilon_res;
    
    /// Brief Correction overal tolerance
    REAL m_epsilon_cor;
    
    /// Brief Number of iteration for fss scheme
    int m_n_fss_iterations;
    
    /// Brief Number of iteration for enforced fss scheme
    int m_n_enf_fss_iterations;
    
    /// Brief Number of thread
    int m_n_threads;
    
    /// Brief Directive that states the use of dual (true) or pirmal (false) formulation for monophacic flow
    bool m_is_dual_formulation_Q;
    
    /// Brief Directive that states if the last memory solution is being transferred to the current memory solution
    bool m_transfer_current_to_last_solution_Q;
    
    /// Brief Spatial refinemenet level
    int m_h_level;
    
    /// Brief Polynomial order for elasticity component
    int m_elasticity_order;
    
    /// Brief Polynomial order for diffusion component
    int m_diffusion_order;
    
    /// Brief Physical dimension of the domain
    int m_dimesion;
    
    /// Brief Name for the Gmsh geometry file being used
    std::string m_geometry_file;
    
    /// Brief Neopz geometry description
    TPZGeoMesh * m_geometry;
    
    /// Brief Name for the vtk files being postprocessed
    std::string m_vtk_file;
    
    /// Brief Number of vtk resolution during postprocessing
    int m_vtk_resolution;
    
    /// Brief Number of geomechanics outputs
    int m_n_outputs_geo;
    
    /// Brief Number of reservoir outputs
    int m_n_outputs_res;
    
    /// Brief Vector that storage scalar names for reservoir postprocessing
    TPZManVector<std::string,50> m_s_names_res;
    
    /// Brief Vector that storage vectors names for reservoir postprocessing
    TPZManVector<std::string,50> m_v_names_res;
    
    /// Brief Vector that storage scalar names for geomechanics postprocessing
    TPZManVector<std::string,50> m_s_names_geo;
    
    /// Brief Vector that storage vectors names for geomechanics postprocessing
    TPZManVector<std::string,50> m_v_names_geo;
    
    /// Brief Vector that storage tensor names for geomechanics postprocessing
    TPZManVector<std::string,50> m_t_names_geo;
    
    /// Brief Gravity field
    TPZManVector<REAL,3> m_g;
    
    /// Brief Integer that define the number of regions presented in the geometry
    int m_n_regions;
    
    /// Brief Material and boundaries identifiers sorted per region
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12> m_mat_ids;
    
    /// Brief Material properties sorted per region
    TPZManVector<std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters>,12> m_mat_props;
    
    /// Controled by the kernel
    
    /// Brief Initial state directive
    bool m_is_initial_state_Q;
    
    /// Brief Current time directive
    bool m_is_current_state_Q;
    
    /// Brief Use for Crank-Nicolson method directive
    bool m_is_crank_nicolson_Q;
        
    /// Brief Map that storage the boundary condition of Geomechanic Simulator identifier with the numerical values provided
    std::map<int, std::vector<REAL> > m_bc_id_to_values_geo_un;
    
    /// Brief Map that storage the provided bc identifiers with the type of boundary condition of Geomechanic Simulator
    std::map<int, std::string> m_bc_id_to_type_geo_un;
    
    /// Brief Map that storage all the boundary conditions of Geomechanic Simulator supported
    std::map< std::string,std::pair<int,std::vector<std::string> > >  m_condition_type_to_index_value_names_geo;
    
    /// Brief Map that storage the boundary condition of Geomechanic Simulator identifier with the numerical values provided
    std::map<int, std::vector<REAL> > m_bc_id_to_values_geo;
    
    /// Brief Map that storage the provided bc identifiers with the type of boundary condition of Geomechanic Simulator
    std::map<int, std::string> m_bc_id_to_type_geo;
    
    /// Brief Map that storage all the boundary conditions of Reservoir Simulator supported
    std::map< std::string,std::pair<int,std::vector<std::string> > >  m_condition_type_to_index_value_names_reser;
    
    /// Brief Map that storage the boundary condition of Reservoir Simulator identifier with the numerical values provided
    std::map<int, std::vector<REAL> > m_bc_id_to_values_reser;
    
    /// Brief Map that storage the provided bc identifiers with the type of boundary condition of Reservoir Simulator
    std::map<int, std::string> m_bc_id_to_type_reser;
    
    /// Brief Directive that states if the current solution must be accepted inside the memory
    bool m_must_accept_solution_Q;
    
    /// Brief update pressure from undrain condition
    bool m_update_pressure_from_undrain_solution_Q;
    
    
public:
    
    
    /// Brief default constructor
    TPMRSSimulationData();
    
    /// Brief Copy constructor
    TPMRSSimulationData(const TPMRSSimulationData & other)
    {
        m_dt                                      = other.m_dt;
        m_n_steps                                 = other.m_n_steps;
        m_reporting_times                         = other.m_reporting_times;
        m_time                                    = other.m_time;
        m_n_iteraions                             = other.m_n_iteraions;
        m_epsilon_res                             = other.m_epsilon_res;
        m_epsilon_cor                             = other.m_epsilon_cor;
        m_n_fss_iterations                        = other.m_n_fss_iterations;
        m_n_enf_fss_iterations                    = other.m_n_enf_fss_iterations;
        m_n_threads                               = other.m_n_threads;
        m_is_dual_formulation_Q                   = other.m_is_dual_formulation_Q;
        m_transfer_current_to_last_solution_Q     = other.m_transfer_current_to_last_solution_Q;
        m_h_level                                 = other.m_h_level;
        m_elasticity_order                        = other.m_elasticity_order;
        m_diffusion_order                         = other.m_diffusion_order;
        m_dimesion                                = other.m_dimesion;
        m_geometry_file                           = other.m_geometry_file;
        m_geometry                                = other.m_geometry;
        m_vtk_file                                = other.m_vtk_file;
        m_vtk_resolution                          = other.m_vtk_resolution;
        m_n_outputs_geo                           = other.m_n_outputs_geo;
        m_n_outputs_res                           = other.m_n_outputs_res;
        m_s_names_res                             = other.m_s_names_res;
        m_s_names_geo                             = other.m_s_names_geo;
        m_v_names_res                             = other.m_v_names_res;
        m_v_names_geo                             = other.m_v_names_geo;
        m_t_names_geo                             = other.m_t_names_geo;
        m_g                                       = other.m_g;
        m_n_regions                               = other.m_n_regions;
        m_mat_ids                                 = other.m_mat_ids;
        m_mat_props                               = other.m_mat_props;
        m_is_initial_state_Q                      = other.m_is_initial_state_Q;
        m_is_current_state_Q                      = other.m_is_current_state_Q;
        m_is_crank_nicolson_Q                     = other.m_is_crank_nicolson_Q;
        m_update_pressure_from_undrain_solution_Q = other.m_update_pressure_from_undrain_solution_Q;
    }
    
    /// Brief Assignement constructor
    TPMRSSimulationData &operator=(const TPMRSSimulationData &other)
    {
        if (this != & other) /// prevent self-assignment
        {
            m_dt                                      = other.m_dt;
            m_n_steps                                 = other.m_n_steps;
            m_reporting_times                         = other.m_reporting_times;
            m_time                                    = other.m_time;
            m_n_iteraions                             = other.m_n_iteraions;
            m_epsilon_res                             = other.m_epsilon_res;
            m_epsilon_cor                             = other.m_epsilon_cor;
            m_n_fss_iterations                        = other.m_n_fss_iterations;
            m_n_enf_fss_iterations                    = other.m_n_enf_fss_iterations;
            m_n_threads                               = other.m_n_threads;
            m_is_dual_formulation_Q                   = other.m_is_dual_formulation_Q;
            m_transfer_current_to_last_solution_Q     = other.m_transfer_current_to_last_solution_Q;
            m_h_level                                 = other.m_h_level;
            m_elasticity_order                        = other.m_elasticity_order;
            m_diffusion_order                         = other.m_diffusion_order;
            m_dimesion                                = other.m_dimesion;
            m_geometry_file                           = other.m_geometry_file;
            m_geometry                                = other.m_geometry;
            m_vtk_file                                = other.m_vtk_file;
            m_vtk_resolution                          = other.m_vtk_resolution;
            m_n_outputs_geo                           = other.m_n_outputs_geo;
            m_n_outputs_res                           = other.m_n_outputs_res;
            m_s_names_res                             = other.m_s_names_res;
            m_s_names_geo                             = other.m_s_names_geo;
            m_v_names_res                             = other.m_v_names_res;
            m_v_names_geo                             = other.m_v_names_geo;
            m_t_names_geo                             = other.m_t_names_geo;
            m_g                                       = other.m_g;
            m_n_regions                               = other.m_n_regions;
            m_mat_ids                                 = other.m_mat_ids;
            m_mat_props                               = other.m_mat_props;
            m_is_initial_state_Q                      = other.m_is_initial_state_Q;
            m_is_current_state_Q                      = other.m_is_current_state_Q;
            m_is_crank_nicolson_Q                     = other.m_is_crank_nicolson_Q;
            m_update_pressure_from_undrain_solution_Q = other.m_update_pressure_from_undrain_solution_Q;
        }

        return *this;
    }
    
    /// Brief destructor
    ~TPMRSSimulationData();
    
    /// Brief Read the xml input file
    void ReadSimulationFile(char *simulation_file);
    
    /// Brief Set the update pressure from undrain condition
    void SetupdatePressureFromUndrainSolutionQ(bool update_pressure_from_undrain_solution_Q) { m_update_pressure_from_undrain_solution_Q = update_pressure_from_undrain_solution_Q; }
    
    /// Brief Set the directive that states if the current memory solution is being transferred to the last memory solution
    void SetTransferCurrentToLastQ(bool transfer_current_to_last_solution_Q) { m_transfer_current_to_last_solution_Q = transfer_current_to_last_solution_Q; }

    /// Brief Get the update pressure from undrain condition
    bool GetupdatePressureFromUndrainSolutionQ() { return m_update_pressure_from_undrain_solution_Q; }
    
    /// Brief Get the directive that states if the current memory solution is being transferred to the last memory solution
    bool GetTransferCurrentToLastQ() { return m_transfer_current_to_last_solution_Q; }
    
    /// Brief Set initial state
    void SetInitialStateQ(bool state) { m_is_initial_state_Q = state; }
    
    /// Brief Get initial state
    bool IsInitialStateQ() {return m_is_initial_state_Q;}
    
    /// Brief Set current time state
    void SetCurrentStateQ(bool state) { m_is_current_state_Q = state; }
    
    /// Brief Get current time state
    bool IsCurrentStateQ() {return m_is_current_state_Q;}

    /// Brief Setup for reporting times and time step size
    void SetTimeControls(int n_times, REAL dt, bool crank_nicolson_Q);
    
    /// Brief Set the directive that states if the current solution must be accepted inside the memory
    bool Set_must_accept_solution_Q(bool must_accept_solution_Q){
        m_must_accept_solution_Q = must_accept_solution_Q;}
    
    /// Brief Set the the use of dual (true) or pirmal (false) formulation for monophacic flow
    bool Set_is_dual_formulation_Q(bool is_dual_formulation_Q){
        m_is_dual_formulation_Q = is_dual_formulation_Q;}
    
    /// Brief Setup for Newton method controls
    void SetNumericControls(int n_iterations, REAL epsilon_res, REAL epsilon_cor);
    
    /// Brief Setup for fixed stress split schemes
    void SetFixedStressSplitSchemes(int n_fss_iterations, int n_enf_fss_iterations);
    
    /// Brief Get time values being reported
    TPZStack< REAL , 500 > ReportingTimes(){
        return m_reporting_times;
    }
    
    /// Brief Time step size
    REAL dt() { return m_dt; }
    
    /// Brief Set the current time value
    void SetTime(REAL time) { m_time = time; }
    
    /// Brief Get the current time value
    REAL t() { return m_time; }
    
    /// Brief Get the number of time steps
    int n_steps() { return m_n_steps; }
    
    /// Brief Get the number of iterations steps
    int n_iterations() { return m_n_iteraions; }
    
    /// Brief Get the residue overal tolerance
    REAL epsilon_res() { return m_epsilon_res; }
    
    /// Brief Get the correction overal tolerance
    REAL epsilon_cor() { return m_epsilon_cor; }
    
    /// Brief Get the maximum number of fixed stress split scheme
    int n_fss_iterations() { return m_n_fss_iterations; }
    
    /// Brief Get the number of enforced for fixed stress split scheme
    int n_enf_fss_iterations() { return m_n_enf_fss_iterations; }
    
    /// Brief Get the number of threads
    int n_threads() { return m_n_threads; }
    
    /// Brief Get Name for the vtk files being postprocessed
    std::string name_vtk_file() { return m_vtk_file; }
    
    /// Brief Get Number of vtk resolution during postprocessing
    int n_div() { return m_vtk_resolution; }
    
    /// Brief Get Number of geomechanics output
    int num_outputs_geo() { return m_n_outputs_geo; }
    
    /// Brief Get Number of reservoir output
    int num_outputs_res() { return m_n_outputs_res; }
    
    /// Brief Get Vector that storage scalar names for reservoir postprocessing
    TPZManVector<std::string,50> s_names_res() { return m_s_names_res; }
    
    /// Brief Get Vector that storage vector names for reservoir postprocessing
    TPZManVector<std::string,50> v_names_res() { return m_v_names_res; }
    
    /// Brief Get Vector that storage scalar names for geomechanics postprocessing
    TPZManVector<std::string,50> s_names_geo() { return m_s_names_geo; }
    
    /// Brief Get Vector that storage vector names for geomechanics postprocessing
    TPZManVector<std::string,50> v_names_geo() { return m_v_names_geo; }
    
    /// Brief Get Vector that storage tensor names for geomechanics postprocessing
    TPZManVector<std::string,50> t_names_geo() { return m_t_names_geo; }
    
    /// Brief Get the gravity field
    TPZVec<REAL> & Gravity()
    {
        return m_g;
    }
    
    /// Brief Get the neopz geometry description
    TPZGeoMesh * Geometry()
    {
        return m_geometry;
    }
    
    /// Brief dimension of the model
    int Dimension() const {return m_dimesion;}

    /// Brief Get the number of regions presented in the geometry
    int NumberOfRegions() { return m_n_regions; }
    
    /// Brief Get the material and boundaries identifiers sorted per region
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12> & MaterialIds() { return m_mat_ids; }
    
    /// Brief Get the material properties sorted per region
    TPZManVector<std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters>,12> & MaterialProps() { return m_mat_props; }
    
    /// Brief Get the physical dimension of the domain
    int Dimension() { return m_dimesion; }
    
    /// Brief Get the spatial refinemenet level
    int HLevel() { return m_h_level; }
    
    /// Brief Get the polynomial order for elasticity component
    int ElasticityOrder() { return m_elasticity_order; }
    
    /// Brief Get the polynomial order for diffusion component
    int DiffusionOrder() { return m_diffusion_order; }
    
    /// Brief Print the all members
    void Print();
    
    /// Brief Print the geometry member
    void PrintGeometry();
    
    /// Brief Get the map that storage all the boundary conditions of Reservoir Simulator supported
    std::map< std::string,std::pair<int,std::vector<std::string> > > & ConditionTypeToBCIndexReservoirs() { return m_condition_type_to_index_value_names_reser; }
    
    /// Brief Get the map that storage the type of boundary condition of Reservoir Simulator with the numerical values provided
    std::map< int , std::vector<REAL> > & BCIdToBCValuesReservoirs() { return m_bc_id_to_values_reser; }
    
    /// Brief Get the map that storage the provided bc identifiers with the type of boundary condition of Reservoir Simulator
    std::map<int, std::string> & BCIdToConditionTypeReservoirs() { return m_bc_id_to_type_reser; }
    
    /// Brief Get the map that storage the type of boundary condition of Geomechanic Simulator with the numerical values provided
    std::map< int , std::vector<REAL> > & BCIdToBCValuesGeomechanicsUndrained() { return m_bc_id_to_values_geo_un; }
    
    /// Brief Get the map that storage the provided bc identifiers with the type of boundary condition of Geomechanic Simulator
    std::map<int, std::string> & BCIdToConditionTypeGeomechanicsUndrained() { return m_bc_id_to_type_geo_un; }
    
    /// Brief Get the map that storage all the boundary conditions of Geomechanic Simulator supported
    std::map< std::string,std::pair<int,std::vector<std::string> > > & ConditionTypeToBCIndexGeomechanics() { return m_condition_type_to_index_value_names_geo; }
    
    /// Brief Get the map that storage the type of boundary condition of Geomechanic Simulator with the numerical values provided
    std::map< int , std::vector<REAL> > & BCIdToBCValuesGeomechanics() { return m_bc_id_to_values_geo; }
    
    /// Brief Get the map that storage the provided bc identifiers with the type of boundary condition of Geomechanic Simulator
    std::map<int, std::string> & BCIdToConditionTypeGeomechanics() { return m_bc_id_to_type_geo; }
    
    /// Brief Get the directive that states if the current solution must be accepted inside the memory
    bool Get_must_accept_solution_Q() { return m_must_accept_solution_Q; }
    
    /// Brief Get the the use of dual (true) or pirmal (false) formulation for monophacic flow
    bool Get_is_dual_formulation_Q() { return m_is_dual_formulation_Q; }
    
    //// Brief Get crank nicolson directive for time derivative (false Euler method)
    bool Get_is_crank_nicolson_Q() { return m_is_crank_nicolson_Q; }
    
    
private:
    
    /// Brief Read the Gmsh file and set the geometry member
    void ReadGeometry();
    
    /// Brief Fillup the map that storage all the boundary conditions of Reservoir Simulator supported
    void LoadBoundaryConditionsReservoirs();
    
    /// Brief Fillup the map that storage all the boundary conditions of Geomechanic Simulator supported
    void LoadBoundaryConditionsGeomechanics();

    /// Brief Apply uniform refinements
    void UniformRefinement();
    
    // @TODO:: MS, please implement and comment this function
    void ReadRegionsAndMaterials();
    
    // @TODO:: MS, please implement and comment this function
    void ReadBCForGeomechanicSimulator();
    
    // @TODO:: MS, please implement and comment this function
    void ReadBCForReservoirSimulator();
    
};


#endif /* TPMRSSimulationData_h */
