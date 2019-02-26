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
#include "TPMRSInterpolator.h"


class TPMRSSimulationData
{
    
protected:
    
    /// Time step size
    REAL m_dt;
    
    /// Number of time steps
    int m_n_steps;
    
    /// Store time values to be reported
    TPZStack< REAL , 500 > m_reporting_times;
    
    /// Current time value
    REAL m_time;
    
    /// Number of iteration
    int m_n_iteraions;
    
    //// Residue overal tolerance
    REAL m_epsilon_res;
    
    /// Correction overal tolerance
    REAL m_epsilon_cor;
    
    /// Number of iteration for fss scheme
    int m_n_fss_iterations;
    
    /// Number of iteration for enforced fss scheme
    int m_n_enf_fss_iterations;
    
    /// The maximum substeps of plasticity
    REAL m_max_plastic_strain;
    
    /// Name of nonlinear acceleration method
    std::string m_n_nonlinear_acceleration;
    
    /// Number of thread
    int m_n_threads;
    
    /// value of scale factor
    REAL m_scale_factor;
    
    /// Directive that states the use of dual (true) or pirmal (false) formulation for monophacic flow
    bool m_is_dual_formulation_Q;
    
    /// Directive that states if the current memory solution is being transferred to the last memory solution
    bool m_transfer_current_to_last_solution_Q;
    
    /// Spatial refinemenet level
    int m_h_level;
    
    /// Polynomial order for elasticity component
    int m_elasticity_order;
    
    /// Polynomial order for diffusion component
    int m_diffusion_order;
    
    /// Physical dimension of the domain
    int m_dimesion;
    
    /// Name for the Gmsh geometry file being used
    std::string m_geometry_file;
    
    /// Neopz geometry description
    TPZGeoMesh * m_geometry;
    
    /// Name for the vtk files being postprocessed
    std::string m_vtk_file;
    
    /// Number of vtk resolution during postprocessing
    int m_vtk_resolution;
    
    /// Number of geomechanics outputs
    int m_n_outputs_geo;
    
    /// Number of reservoir outputs
    int m_n_outputs_res;
    
    /// Vector that storage scalar names for reservoir postprocessing
    TPZManVector<std::string,50> m_s_names_res;
    
    /// Vector that storage vectors names for reservoir postprocessing
    TPZManVector<std::string,50> m_v_names_res;
    
    /// Vector that storage scalar names for geomechanics postprocessing
    TPZManVector<std::string,50> m_s_names_geo;
    
    /// Vector that storage vectors names for geomechanics postprocessing
    TPZManVector<std::string,50> m_v_names_geo;
    
    /// Vector that storage tensor names for geomechanics postprocessing
    TPZManVector<std::string,50> m_t_names_geo;
    
    /// Gravity field
    TPZManVector<REAL,3> m_g;
    
    /// Integer that define the number of regions presented in the geometry
    int m_n_regions;
    
    /// Material and boundaries identifiers sorted per region
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12> m_mat_ids;
    
    /// Material properties sorted per region
    TPZManVector<std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters>,12> m_mat_props;
    
    /// Controled by the kernel
    
    /// Initial state directive
    bool m_is_initial_state_Q;
    
    /// Current time directive
    bool m_is_current_state_Q;
    
    /// Use for Crank-Nicolson method directive
    bool m_is_crank_nicolson_Q;
        
    /// Map that storage the boundary condition of Geomechanic Simulator identifier with the numerical values provided (initial)
    std::map<int, TPMRSInterpolator > m_bc_id_to_values_geo_initial;
    
    /// Map that storage the provided bc identifiers with the type of boundary condition of Geomechanic Simulator (initial)
    std::map<int, std::string> m_bc_id_to_type_geo_initial;
    
    /// Map that storage all the boundary conditions of Geomechanic Simulator supported
    std::map< std::string,std::pair<int,std::vector<std::string> > >  m_condition_type_to_index_value_names_geo;
    
    /// Map that storage the boundary condition of Geomechanic Simulator identifier with the numerical values provided
    std::map<int, TPMRSInterpolator > m_bc_id_to_values_geo;
    
    /// Map that storage the provided bc identifiers with the type of boundary condition of Geomechanic Simulator
    std::map<int, std::string> m_bc_id_to_type_geo;
    
    /// Map that storage all the boundary conditions of Reservoir Simulator supported
    std::map< std::string,std::pair<int,std::vector<std::string> > >  m_condition_type_to_index_value_names_res;
    
    /// Map that storage the boundary condition of Reservoir Simulator identifier with the numerical values provided (initial)
    std::map<int, TPMRSInterpolator > m_bc_id_to_values_res_initial;
    
    /// Map that storage the provided bc identifiers with the type of boundary condition of Reservoir Simulator (initial)
    std::map<int, std::string> m_bc_id_to_type_res_initial;
    
    /// Map that storage the boundary condition of Reservoir Simulator identifier with the numerical values provided
    std::map<int, TPMRSInterpolator > m_bc_id_to_values_res;
    
    /// Map that storage the provided bc identifiers with the type of boundary condition of Reservoir Simulator
    std::map<int, std::string> m_bc_id_to_type_res;
    
    /// Directive that states if the current solution must be accepted inside the memory
    bool m_must_accept_solution_Q;
    
    /// update pressure from undrain condition
    bool m_update_pressure_from_undrain_solution_Q;
    
    /// Use substeps during solving geomecanics
    bool m_must_use_sub_stepping_Q;
    
    /// Number of level of substeps during solving geomecanics
    int m_n_sub_step_level;
    
    /// String that stands for the used acceleration method (None,Aitken,..,ect)
    std::string m_nonlinear_acceleration;
    
    
public:
    
    
    /// default constructor
    TPMRSSimulationData();
    
    /// Copy constructor
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
        m_n_nonlinear_acceleration                = other.m_n_nonlinear_acceleration;
        m_max_plastic_strain                      = other.m_max_plastic_strain;
        m_n_threads                               = other.m_n_threads;
        m_scale_factor                            = other.m_scale_factor;
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
        m_must_use_sub_stepping_Q                 = other.m_must_use_sub_stepping_Q;
        m_n_sub_step_level                        = other.m_n_sub_step_level;
    }
    
    /// Assignement constructor
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
            m_n_nonlinear_acceleration                = other.m_n_nonlinear_acceleration;
            m_max_plastic_strain                      = other.m_max_plastic_strain;
            m_n_threads                               = other.m_n_threads;
            m_scale_factor                            = other.m_scale_factor;
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
            m_must_use_sub_stepping_Q                 = other.m_must_use_sub_stepping_Q;
            m_n_sub_step_level                        = other.m_n_sub_step_level;
        }

        return *this;
    }
    
    /// destructor
    ~TPMRSSimulationData();
    
    /// Read the xml input file
    void ReadSimulationFile(char *simulation_file);
    
    /// Set the update pressure from undrain condition
    void SetupdatePressureFromUndrainSolutionQ(bool update_pressure_from_undrain_solution_Q) { m_update_pressure_from_undrain_solution_Q = update_pressure_from_undrain_solution_Q; }
    
    /// Get the update pressure from undrain condition
    bool GetupdatePressureFromUndrainSolutionQ() { return m_update_pressure_from_undrain_solution_Q; }
    
    /// Set the directive that states if the current memory solution is being transferred to the last memory solution
    void SetTransferCurrentToLastQ(bool transfer_current_to_last_solution_Q) { m_transfer_current_to_last_solution_Q = transfer_current_to_last_solution_Q; }

    
    /// Get the directive that states if the current memory solution is being transferred to the last memory solution
    bool GetTransferCurrentToLastQ() { return m_transfer_current_to_last_solution_Q; }
    
    /// Set initial state
    void SetInitialStateQ(bool state) { m_is_initial_state_Q = state; }
    
    /// Get initial state
    bool IsInitialStateQ() {return m_is_initial_state_Q;}
    
    /// Set current time state
    void SetCurrentStateQ(bool state) { m_is_current_state_Q = state; }
    
    /// Get current time state
    bool IsCurrentStateQ() {return m_is_current_state_Q;}

    /// Setup for reporting times and time step size
    void SetTimeControls(int n_times, REAL dt, bool crank_nicolson_Q);
    
    /// Set the directive that states if the current solution must be accepted inside the memory
    void Set_must_accept_solution_Q(bool must_accept_solution_Q){
        m_must_accept_solution_Q = must_accept_solution_Q;
    }
    
    /// Set the use of dual (true) or primal (false) formulation for monophacic flow
    void Set_is_dual_formulation_Q(bool is_dual_formulation_Q){
        m_is_dual_formulation_Q = is_dual_formulation_Q;
    }
    
    /// Setup for Newton method controls
    void SetNumericControls(int n_iterations, REAL epsilon_res, REAL epsilon_cor);
    
    /// Setup for fixed stress split schemes
    void SetFixedStressSplitSchemes(int n_fss_iterations, int n_enf_fss_iterations);
    
    /// Get time values being reported
    TPZStack< REAL , 500 > ReportingTimes(){
        return m_reporting_times;
    }
    
    /// Set time step size
    void Setdt(REAL & dt) { m_dt = dt; }
    
    /// Get time step size
    REAL dt() { return m_dt; }
    
    /// Set the current time value
    void SetTime(REAL time) { m_time = time; }
    
    /// Get the current time value
    REAL t() { return m_time; }
    
    /// Get the number of time steps
    int n_steps() { return m_n_steps; }
    
    /// Get the number of iterations steps
    int n_iterations() { return m_n_iteraions; }
    
    /// Get the residue overal tolerance
    REAL epsilon_res() { return m_epsilon_res; }
    
    /// Get the correction overal tolerance
    REAL epsilon_cor() { return m_epsilon_cor; }
    
    /// Get the maximum number of fixed stress split scheme
    int n_fss_iterations() { return m_n_fss_iterations; }
    
    /// Get the number of enforced for fixed stress split scheme
    int n_enf_fss_iterations() { return m_n_enf_fss_iterations; }
    
    /// Get Name for the nonlinear acceleration
    std::string name_nonlinear_acceleration() { return m_n_nonlinear_acceleration; }
    
    /// Get the maximum substeps of plasticity
    REAL Get_max_plastic_strain() { return m_max_plastic_strain; }
    
    /// Get the number of threads
    int n_threads() { return m_n_threads; }
    
    /// Get the value of scale factor
    REAL scale_factor_val() { return m_scale_factor; }
    
    /// Get Name for the vtk files being postprocessed
    std::string name_vtk_file() { return m_vtk_file; }
    
    /// Get Number of vtk resolution during postprocessing
    int n_div() { return m_vtk_resolution; }
    
    /// Get Number of geomechanics output
    int num_outputs_geo() { return m_n_outputs_geo; }
    
    /// Get Number of reservoir output
    int num_outputs_res() { return m_n_outputs_res; }
    
    /// Get Vector that storage scalar names for reservoir postprocessing
    TPZManVector<std::string,50> s_names_res() { return m_s_names_res; }
    
    /// Get Vector that storage vector names for reservoir postprocessing
    TPZManVector<std::string,50> v_names_res() { return m_v_names_res; }
    
    /// Get Vector that storage scalar names for geomechanics postprocessing
    TPZManVector<std::string,50> s_names_geo() { return m_s_names_geo; }
    
    /// Get Vector that storage vector names for geomechanics postprocessing
    TPZManVector<std::string,50> v_names_geo() { return m_v_names_geo; }
    
    /// Get Vector that storage tensor names for geomechanics postprocessing
    TPZManVector<std::string,50> t_names_geo() { return m_t_names_geo; }
    
    /// Get the gravity field
    TPZVec<REAL> & Gravity()
    {
        return m_g;
    }
    
    /// Get the neopz geometry description
    TPZGeoMesh * Geometry()
    {
        return m_geometry;
    }
    
    /// dimension of the model
    int Dimension() const {return m_dimesion;}

    /// Get the number of regions presented in the geometry
    int NumberOfRegions() { return m_n_regions; }
    
    /// Get the material and boundaries identifiers sorted per region
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12> & MaterialIds() { return m_mat_ids; }
    
    /// Get the material properties sorted per region
    TPZManVector<std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters>,12> & MaterialProps() { return m_mat_props; }
    
    /// Get the physical dimension of the domain
    int Dimension() { return m_dimesion; }
    
    /// Get the spatial refinemenet level
    int HLevel() { return m_h_level; }
    
    /// Get the polynomial order for elasticity component
    int ElasticityOrder() { return m_elasticity_order; }
    
    /// Get the polynomial order for diffusion component
    int DiffusionOrder() { return m_diffusion_order; }
    
    /// Print the all members
    void Print();
    
    /// Print the geometry member
    void PrintGeometry();
    
    /// Get the map that storage all the boundary conditions of Reservoir Simulator supported
    std::map< std::string,std::pair<int,std::vector<std::string> > > & ConditionTypeToBCIndexReservoir() { return m_condition_type_to_index_value_names_res; }
    
    /// Get the map that storage the type of boundary condition of Reservoir Simulator with the numerical values provided
    std::map< int , TPMRSInterpolator > & BCIdToBCValuesReservoir() { return m_bc_id_to_values_res; }
    
    /// Get the map that storage the provided bc identifiers with the type of boundary condition of Reservoir Simulator
    std::map<int, std::string> & BCIdToConditionTypeReservoir() { return m_bc_id_to_type_res; }
    
    /// Get the map that storage the type of boundary condition of Geomechanic Simulator with the numerical values provided (initial)
    std::map< int , TPMRSInterpolator > & BCIdToBCValuesGeomechanicsInitial() { return m_bc_id_to_values_geo_initial; }
    
    /// Get the map that storage the provided bc identifiers with the type of boundary condition of Geomechanic Simulator (initial)
    std::map<int, std::string> & BCIdToConditionTypeGeomechanicsInitial() { return m_bc_id_to_type_geo_initial; }
    
    /// Get the map that storage the type of boundary condition of Reservoir Simulator with the numerical values provided
    std::map< int , TPMRSInterpolator > & BCIdToBCValuesReservoirInitial() { return m_bc_id_to_values_res_initial; }
    
    /// Get the map that storage the provided bc identifiers with the type of boundary condition of Reservoir Simulator
    std::map<int, std::string> & BCIdToConditionTypeReservoirInitial() { return m_bc_id_to_type_res_initial; }
    
    /// Get the map that storage all the boundary conditions of Geomechanic Simulator supported
    std::map< std::string,std::pair<int,std::vector<std::string> > > & ConditionTypeToBCIndexGeomechanics() { return m_condition_type_to_index_value_names_geo; }
    
    /// Get the map that storage the type of boundary condition of Geomechanic Simulator with the numerical values provided
    std::map<int, TPMRSInterpolator > & BCIdToBCValuesGeomechanics() { return m_bc_id_to_values_geo; }
    
    /// Get the map that storage the provided bc identifiers with the type of boundary condition of Geomechanic Simulator
    std::map<int, std::string> & BCIdToConditionTypeGeomechanics() { return m_bc_id_to_type_geo; }
    
    /// Get the directive that states if the current solution must be accepted inside the memory
    bool Get_must_accept_solution_Q() { return m_must_accept_solution_Q; }
    
    /// Get the use of dual (true) or primal (false) formulation for monophacic flow
    bool Get_is_dual_formulation_Q() { return m_is_dual_formulation_Q; }
    
    //// Get crank nicolson directive for time derivative (false Euler method)
    bool Get_is_crank_nicolson_Q() { return m_is_crank_nicolson_Q; }
    
    //// Set directive for using substeps during solving geomecanics
    void Set_must_use_sub_stepping_Q(bool must_use_sub_stepping_Q) {
        m_must_use_sub_stepping_Q = must_use_sub_stepping_Q;
    }
    
    //// Get directive for using substeps during solving geomecanics
    bool Get_must_use_sub_stepping_Q() { return m_must_use_sub_stepping_Q; }
    
    //// Set level of substeps during solving geomecanics
    void Set_n_sub_step_level(int n_sub_step_level) {
        m_n_sub_step_level = n_sub_step_level;
    }
    
    //// Get level of substeps during solving geomecanics
    int Get_n_sub_step_level() { return m_n_sub_step_level; }
    
private:
    
    /// Read the Gmsh file and set the geometry member
    void ReadGeometry();
    
    /// Fillup the map that storage all the boundary conditions of Reservoir Simulator supported
    void LoadBoundaryConditionsReservoirs();
    
    /// Fillup the map that storage all the boundary conditions of Geomechanic Simulator supported
    void LoadBoundaryConditionsGeomechanics();

    /// Apply uniform refinements
    void UniformRefinement();
    
    /// @TODO:: MS, please implement and comment this function
    void ReadRegionsAndMaterials();
    
    /// @TODO:: MS, please implement and comment this function
    void ReadBCForGeomechanicSimulator();
    
    /// @TODO:: MS, please implement and comment this function
    void ReadBCForReservoirSimulator();
    
};


#endif /* TPMRSSimulationData_h */
