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
    
    /// Name of nonlinear Newton methods
    std::string m_nonlinear_Newton_method;
    
    /// Use Secant method for reservoir
    bool m_is_secant_reservoir_Q;
    
    /// Use Secant method for geomechanics
    bool m_is_secant_geomechanics_Q;
    
    /// Number of update jacobian reservoir
    int m_n_update_jac_res;
    
    /// Number of update jacobian geomechanics
    int m_n_update_jac_geo;
    
    
    
    /// Number of iteration for fss scheme
    int m_n_fss_iterations;
    
    /// Number of iteration for enforced fss scheme
    int m_n_enf_fss_iterations;
    
    /// The maximum substeps of plasticity
    REAL m_max_plastic_strain;
    
    /// Directive that states the reset undarined respose data
    bool m_reset_undarined_respose_data_Q;
    
    /// Name of nonlinear acceleration method
    std::string m_n_nonlinear_acceleration;
    
    /// Number of states of acceleration method
    int m_n_state;
    
    REAL m_max_theta_value;
    
    /// Number of thread
    int m_n_threads;
    
    /// value of scale factor
    REAL m_scale_factor;
    
    /// Directive that states the use of fully coupled solver or SFI
    bool m_is_fully_coupled_Q;
    
    /// Directive that states the use of dual (true) or primal (false) formulation for monophacic flow
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
    
    /// Directive that states to print the initial data
    bool m_draw_initial_data_Q;
    
    /// Directive that states to print geometry or not
    bool m_is_draw_geometry_Q;
    
    /// Directive that states to show summary of performance
    bool m_performance_summary_Q;
    
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
    
    
    /// Map that storage all the boundary conditions of Fully coupled Simulator supported
    std::map< std::string,std::pair<int,std::vector<std::string> > >  m_condition_type_to_index_value_names_fc;
    
    /// Map that storage the boundary condition of Fully coupled Simulator identifier with the numerical values provided
    std::map<int, TPMRSInterpolator > m_bc_id_to_values_fc;
    
    /// Map that storage the provided bc identifiers with the type of boundary condition of Fully coupled Simulator
    std::map<int, std::string> m_bc_id_to_type_fc;
    
    
    /// Directive that states if the current solution must be accepted inside the memory
    bool m_must_accept_solution_Q;
    
    /// update pressure from undrain condition
    bool m_update_pressure_from_undrain_solution_Q;
    
    /// Use substeps during solving geomecanics
    bool m_must_use_sub_stepping_Q;
    
    /// Number of level of substeps during solving geomecanics
    int m_n_sub_step_level;
    
    /// String that stands for the used acceleration method (None,FDM,SDM)
    std::string m_nonlinear_acceleration;
    
public:
    
    
    /// default constructor
    TPMRSSimulationData();
    
    /// Copy constructor
    TPMRSSimulationData(const TPMRSSimulationData & other);
    
    /// Assignement constructor
    TPMRSSimulationData &operator=(const TPMRSSimulationData &other);
    
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

    
    /// Set the use of fully coupled solver or SFI
    void Set_is_fully_coupled_Q(bool is_fully_coupled_Q){
        m_is_fully_coupled_Q = is_fully_coupled_Q;
    }
    
    /// Set the use of dual (true) or primal (false) formulation for monophacic flow
    void Set_is_dual_formulation_Q(bool is_dual_formulation_Q){
        m_is_dual_formulation_Q = is_dual_formulation_Q;
    }
    
    /// Setup for Newton method controls
    void SetNumericControls(int n_iterations, REAL epsilon_res, REAL epsilon_cor);
    
    /// Setup for fixed stress split schemes
    void SetFixedStressSplitSchemes(int n_fss_iterations, int n_enf_fss_iterations);
    
    /// Set time values being reported
    void SetReportingTimes(TPZStack< REAL , 500 > & reporting_times){
        m_reporting_times = reporting_times;
    }
    
    /// Get time values being reported
    TPZStack< REAL , 500 > & ReportingTimes(){
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
    
    /// Set the number of time steps
    void Set_n_steps(int nstep){m_n_steps = nstep;}
    
    /// Get the number of time steps
    int n_steps() { return m_n_steps; }
    
    /// Get the number of iterations steps
    int n_iterations() { return m_n_iteraions; }
    
    /// Get the residue overal tolerance
    REAL epsilon_res() { return m_epsilon_res; }
    
    /// Get the correction overal tolerance
    REAL epsilon_cor() { return m_epsilon_cor; }
    
    /// Get Name for the nonlinear Newton method
    std::string name_nonlinear_Newton_method() {return m_nonlinear_Newton_method;}
    
    //// Get Secant method for reservoir
    bool Get_is_secant_reservoir_Q() { return m_is_secant_reservoir_Q; }
    
    //// Set Secant method for reservoir
    void Set_is_secant_reservoir_Q(bool is_secant_reservoir_Q){
        m_is_secant_reservoir_Q = is_secant_reservoir_Q;
    }
    
    //// Get Secant method for geomechanics
    bool Get_is_secant_geomechanics_Q() { return m_is_secant_geomechanics_Q; }
    
    //// Set Secant method for geomechanics
    void Set_is_secant_geomechanics_Q(bool is_secant_geomechanics_Q){
        m_is_secant_geomechanics_Q = is_secant_geomechanics_Q;
    }
    
    /// Setup for Secant numerical method controls
    void SetSecantMethod(bool is_secant_reservoir_Q, bool is_secant_geomechanics_Q);
    
    /// Get the number of update jacobian reservoir
    int Get_n_update_jac_res() { return m_n_update_jac_res; }
    
    /// Get the number of update jacobian geomechanics
    int Get_n_update_jac_geo() { return m_n_update_jac_geo; }
    
    
    /// Setup for update Jacobian numerical method controls
    void SetUpdateJacobianMethod(int num_update_jac_res,int num_update_jac_geo);
    
    
    /// Get the maximum number of fixed stress split scheme
    int n_fss_iterations() { return m_n_fss_iterations; }
    
    /// Get the number of enforced for fixed stress split scheme
    int n_enf_fss_iterations() { return m_n_enf_fss_iterations; }
    
    /// Get Name for the nonlinear acceleration
    std::string name_nonlinear_acceleration() { return m_n_nonlinear_acceleration; }
    
    /// Get number of state for acceleration methods
    int n_state_acceleration() {return m_n_state; }
    
    /// Set the maximum theta value for adaptive acceleration
    void Set_max_theta_value (REAL maximum_theta_value) {
        m_max_theta_value = maximum_theta_value;
    }
    
    /// Get the maximum theta value for adaptive acceleration
    REAL Get_max_theta_value() { return m_max_theta_value; }
    
    /// Get the maximum substeps of plasticity
    REAL Get_max_plastic_strain() { return m_max_plastic_strain; }
    
    /// Get the ask to reset undarined respose data
    bool Get_reset_undarined_respose_data_Q() { return m_reset_undarined_respose_data_Q; }
    
    /// Set the ask to reset undarined respose data
    void Set_reset_undarined_respose_data_Q(bool is_reset_undarined_respose_dataQ){
        m_reset_undarined_respose_data_Q = is_reset_undarined_respose_dataQ;
    }
    
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
    
    
    /// Set the ask to print the initial data
    void Set_is_draw_initial_data_Q(bool is_draw_initial_data_Q){
        m_draw_initial_data_Q = is_draw_initial_data_Q;
    }
    
    /// Get the ask to print the initial data
    bool Get_is_draw_initial_data_Q() { return m_draw_initial_data_Q; }
    
    
    /// Set the ask to print the geometry or not
    void Set_is_draw_geometry_Q(bool is_draw_geometry_Q){
        m_is_draw_geometry_Q = is_draw_geometry_Q;
    }
    
    /// Get the ask to print the geometry or not
    bool Get_is_draw_geometry_Q() { return m_is_draw_geometry_Q; }
    
    
    /// Set the ask to show the performance
    void Set_is_performance_summary_Q(bool is_performance_summary_Q){
        m_performance_summary_Q = is_performance_summary_Q;
    }
    
    /// Get the ask to show the performance
    bool Get_is_performance_summary_Q() { return m_performance_summary_Q; }
    
    
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
    
    /// Get the map that storage all the boundary conditions of Fully coupled Simulator supported
    std::map< std::string,std::pair<int,std::vector<std::string> > > & ConditionTypeToBCIndexFullyCoupled() { return m_condition_type_to_index_value_names_fc; }
    
    /// Get the map that storage the type of boundary condition of Fully coupled Simulator with the numerical values provided
    std::map<int, TPMRSInterpolator > & BCIdToBCValuesFullyCoupled() { return m_bc_id_to_values_fc; }
    
    /// Get the map that storage the provided bc identifiers with the type of boundary condition of Fully coupled Simulator
    std::map<int, std::string> & BCIdToConditionTypeFullyCoupled() { return m_bc_id_to_type_fc; }
    
    /// Get the directive that states if the current solution must be accepted inside the memory
    bool Get_must_accept_solution_Q() { return m_must_accept_solution_Q; }
    
    /// Get the use of fully coupled solver or SFI
    bool Get_is_fully_coupled_Q() {return m_is_fully_coupled_Q;}
    
    /// Get the use of dual (true) or primal (false) formulation for monophacic flow
    bool Get_is_dual_formulation_Q() { return m_is_dual_formulation_Q; }
    
    //// Get crank nicolson directive for time derivative (false Euler method)
    bool Get_is_crank_nicolson_Q() { return m_is_crank_nicolson_Q; }
    
    //// Set crank nicolson directive for time derivative (false Euler method)
    void Set_is_crank_nicolson_Q(bool is_crank_nicolson_Q){
        m_is_current_state_Q = is_crank_nicolson_Q;
    }
    
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
    
    /// Fillup the map that storage all the boundary conditions of Fully coupled Simulator supported
    void LoadBoundaryConditionsFC();

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
