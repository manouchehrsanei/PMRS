
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>
#include <tuple>
#include <cmath>
#include <iostream>
#include <math.h>

/// Geometry
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"
#include "tpzhierarquicalgrid.h"

/// Computational mesh
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzintel.h"
#include "pzlog.h"

/// Materials
#include "pzl2projection.h"
#include "pzbndcond.h"
#include "TPZBndCondWithMem.h"

/// Monophasic
#include "TPMRSMonoPhasic_impl.h"
#include "TPMRSMonoPhasicCG_impl.h"
#include "TPMRSMonoPhasicAnalysis.h"
#include "TPMRSMonoPhasicMemory.h"


/// Geomechanics
#include "TPMRSElastoPlastic_impl.h"
#include "TPMRSGeomechanicAnalysis.h"
#include "TPMRSElastoPlasticMemory.h"

/// Segregated solver
#include "TPMRSMemory.h"
#include "TPMRSSegregatedAnalysis.h"


/// Elasticity
#include "TPZElasticCriterion.h"

/// Plasticity
#include "TPZPlasticStepPV.h"
#include "TPZSandlerExtended.h"
#include "TPZYCMohrCoulombPV.h"

#include "pzporoelastoplasticmem.h"
#include "TPZSandlerDimaggio.h"

/// Runge-Kutta solver
#include "TPMRSRKSolver_impl.h"

/// Analysis
#include "pzpostprocanalysis.h"
#include "pzanalysis.h"
#include "TPMRSFullyCoupledAnalysis.h"


/// Matrix
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzfstrmatrix.h"

/// Simulation data structure
#include "TPMRSSimulationData.h"
#include "TPMRSUndrainedParameters.h"
#include "TPMRSPoroMechParameters.h"
#include "TPMRSPhiParameters.h"
#include "TPMRSKappaParameters.h"
#include "TPMRSPlasticityParameters.h"

/// Poroelastoplastic material by means of monolithic approach
#include "TPMRSPoroElastoPlastic_impl.h"

#ifdef LOG4CXX
static LoggerPtr log_data(Logger::getLogger("pz.PMRS"));
#endif

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


/// PoroElastic Full Coupling
// H1 mesh for Deformation
TPZCompMesh * CMesh_Deformation(TPMRSSimulationData * sim_data);
// H1 mesh for Pore Pressure
TPZCompMesh * CMesh_PorePressure(TPMRSSimulationData * sim_data);
// Multiphysics Coupling
TPZCompMesh * CMesh_FullyCoupled(TPZManVector<TPZCompMesh * , 2 > & mesh_vector, TPMRSSimulationData * sim_data);

/// Monophasic Reservoir Simulator
// Hdiv mesh
TPZCompMesh * CMesh_Flux(TPMRSSimulationData * sim_data);
// L2 mesh
TPZCompMesh * CMesh_PorePressure_disc(TPMRSSimulationData * sim_data);
// Mixed mesh
TPZCompMesh * CMesh_Mixed(TPZManVector<TPZCompMesh * , 2 > & mesh_vector, TPMRSSimulationData * sim_data);
TPZMaterial * ConfigurateAndInsertVolumetricMaterialsRes(bool IsMixedQ, int index, int matid, TPMRSSimulationData * sim_data, TPZCompMesh * cmesh);

// CG mesh for pressure
TPZCompMesh * CMesh_Primal(TPMRSSimulationData * sim_data);
// Configurate and insert volumetric materials

// Geomechanic Simulator
// H1 mesh for displacements
TPZCompMesh * CMesh_Geomechanics(TPMRSSimulationData * sim_data);
// Configurate and insert volumetric materials
TPZMaterial * ConfigurateAndInsertVolumetricMaterialsGeo(int index, int matid, TPMRSSimulationData * sim_data, TPZCompMesh * cmesh);

// Method that creates a fully coupled solver
TPMRSFullyCoupledAnalysis * CreateFCSolver(TPMRSSimulationData * sim_data);
// Configurate and insert volumetric materials
TPZMaterial * ConfigurateAndInsertVolumetricMaterialsFC(int index, int matid, TPMRSSimulationData * sim_data, TPZCompMesh * cmesh);

/// Restructuring implementation of Reservoir Geomechanics Simulator

// Method that makes use of TPMRSMonophasic and TPMRSElastoPlastic with a common memory
TPMRSSegregatedAnalysis * CreateSFISolver(TPMRSSimulationData * sim_data);


/// Deprecated.

// Method that makes use of TPMRSMonophasic for parabolic solutions
void RuningMonophasic(TPMRSSimulationData * sim_data);

// Method that makes use of TPMRSElastoPlastic
void RuningGeomechanics(TPMRSSimulationData * sim_data);

/// Shear-enhanced compaction and strain localization:
// Inelastic deformation and constitutive modeling of four porous sandstones
// This function generate the data associated to the Figure 2a. Darley Dale Sandstone Data for Cap Model
/// Begin
void LEDSPorosityReductionPlot();

// Method that apply stress
void Apply_Stress(TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> &LEDS, TPZFMatrix<REAL> &De, TPZFMatrix<REAL> &De_inv, TPZTensor<REAL> &sigma, TPZTensor<REAL> &epsilon);

// Read experimental duplet
TPZFMatrix<REAL> Read_Duplet(int n_data, std::string file);
/// End

void RunRKApproximation(TPMRSSimulationData * sim_data);

int main(int argc, char *argv[])
{
    
    bool RK_approximation_Q = true;

#ifdef LOG4CXX
    if(log_data->isInfoEnabled())
    {
        std::stringstream sout;
        sout << " Defining the case MS -> Review... " << std::endl;
        LOGPZ_DEBUG(log_data,sout.str())
    }
#endif
    
    //    Reading arguments
    char *simulation_file = NULL;
    {
        using namespace std;
        if (argc != 2)
        {
            cout << "Size: " << argc << " Number of Arguments " << endl;
            cout << "Usage: " << argv[0] << " my_input_file.xml " << endl;
            cout <<    "Program must stop: not xml file found \n"    << endl;
            DebugStop();
        }
        
        if (argc == 2)
        {
            cout << "Control File used : " << argv[1]  << "\n"  << endl;
            simulation_file        = argv[1];
        }
    }
    
    // Simulation data to be configurated
    TPMRSSimulationData * sim_data = new TPMRSSimulationData;
    sim_data->ReadSimulationFile(simulation_file);
    
    if (RK_approximation_Q) {
        RunRKApproximation(sim_data);
        return 0;
    }
    
    bool is_fully_coupled_Q = sim_data->Get_is_fully_coupled_Q();
    
#ifdef USING_BOOST
    boost::posix_time::ptime int_case_t1 = boost::posix_time::microsec_clock::local_time();
#endif

    if (is_fully_coupled_Q) {
        
        TPMRSSegregatedAnalysis * SFI_analysis = CreateSFISolver(sim_data);
        TPMRSFullyCoupledAnalysis * FC_analysis = CreateFCSolver(sim_data);
        // Liking the memory to FC solver
        SFI_analysis->ApplyMemoryLink(SFI_analysis->GetGeomechanicsSolver()->Mesh(), FC_analysis->Mesh());
        
        REAL t_0 = 0;
        SFI_analysis->ConfigureGeomechanicsBC(t_0,true);
        SFI_analysis->ConfigureReservoirBC(t_0,true);
        
        /// vtk file
        std::string name = sim_data->name_vtk_file();
        std::string file = name + "_fc.vtk";
        { /// Initial states and postprocess them.
            /// Compute initial response
            SFI_analysis->ExecuteStaticSolution();
            FC_analysis->PostProcessTimeStep(file);
            
            /// Compute undrained response
            SFI_analysis->ExecuteUndrainedStaticSolution();
            FC_analysis->PostProcessTimeStep(file);
        }
        // Load initial conditions in FC for dof
        FC_analysis->Meshvec()[0]->Solution() = SFI_analysis->GetGeomechanicsSolver()->Solution();
        FC_analysis->Meshvec()[1]->Solution() = SFI_analysis->GetReservoirSolver()->Solution();
        TPZBuildMultiphysicsMesh::TransferFromMeshes(FC_analysis->Meshvec(), FC_analysis->Mesh());
        FC_analysis->Solution()=FC_analysis->Mesh()->Solution();
        
        { ///  Printing FC bc conditions
            
            std::cout << "Begining:: Printing FC bc conditions " << std::endl;
            for(auto i :sim_data->BCIdToConditionTypeFullyCoupled()){
                std::cout << "BC id = " << i.first << std::endl;
                std::cout << "BC type = " << i.second << std::endl;
                std::cout << "BC index = " << sim_data->ConditionTypeToBCIndexFullyCoupled()[i.second].first << std::endl;
            }
            std::cout << "Ending:: Printing FC bc conditions " << std::endl;
        }
        
        FC_analysis->ExecuteTimeEvolution();
        
    }
    else
    {

        TPMRSSegregatedAnalysis * SFI_analysis = CreateSFISolver(sim_data);
        REAL t_0 = 0;
        SFI_analysis->ConfigureGeomechanicsBC(t_0,true);
        SFI_analysis->ConfigureReservoirBC(t_0,true);
        SFI_analysis->ExecuteStaticSolution();
        SFI_analysis->ExecuteUndrainedStaticSolution();
        SFI_analysis->ExecuteTimeEvolution();
        
        /// Writing summaries
        TPZFMatrix<REAL> iterations = SFI_analysis->IterationsSummary();
        TPZFMatrix<REAL> residuals  = SFI_analysis->ResidualsSummary();
        TPZFMatrix<REAL> cpu_time   = SFI_analysis->TimeSummary();
        
        if (sim_data->Get_is_performance_summary_Q()) {
            std::string name = sim_data->name_vtk_file();
            std::string file_performance = name + "_performance_summary.txt";
            std::ofstream summary_file(file_performance.c_str());
            iterations.Print("iterations = ",summary_file,EMathematicaInput);
            residuals.Print("residuals = ",summary_file,EMathematicaInput);
            cpu_time.Print("time = ",summary_file,EMathematicaInput);
            summary_file.flush();
        }
    }
    
#ifdef USING_BOOST
    boost::posix_time::ptime int_case_t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    REAL case_solving_time = boost::numeric_cast<double>((int_case_t2-int_case_t1).total_milliseconds());
    std::cout << "Case closed in :" << setw(10) <<  case_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
    std::cout << std::endl;
#endif
    

    
	return EXIT_SUCCESS;
}

void RunRKApproximation(TPMRSSimulationData * sim_data){
    
    std::ofstream rk_file_data("rk_data.txt");
    /// Discretization
    int n_steps = 1000;
    REAL rw = 0.1;
    REAL re = 10.0;
    
    TPZTensor<REAL> sigma_0;
    REAL p_0,u_r,sigma_r,sigma_t,sigma_z,p_r,q_r;
    REAL eps_p_r,eps_p_t,eps_p_z;
    sigma_0.Zero();
    
    
    sigma_0.XX() = -10.0;
    sigma_0.YY() = -10.0;
    sigma_0.ZZ() = -4.0;
    p_0 = 30.0;
    
    int A = 30;
    switch (A) {
        case 0:
        {
            
            u_r        = -0.000495691007;
            sigma_r    = 0.0;
            sigma_t    = -30.85560036;
            sigma_z    = -6.25865984;
            p_r        = 20.0;
            q_r        = -2.171439886;
            eps_p_r    = 0.001573640038;
            eps_p_t    = -0.0005245369975;
            eps_p_z    = 0.0;
            
        }
            break;
        case 10:
        {
            
//{-0.0004978550132,19.99970055,-1.964660048,-0.4517470002,-30.87969971,-6.266290188,0.00152619998,-0.0005087259924,-1.550029968*10^(-18)}
            
            u_r        = -0.0004949550132;
            sigma_r    = 0.0;
            sigma_t    = -30.87969971;
            sigma_z    = -6.266290188;
            p_r        = 20.0;
            q_r        = -1.964660048;
            eps_p_r    = 0.00152619998;
            eps_p_t    = -0.0005087259924;
            eps_p_z    = 0.0;
            
        }
            break;
        case 20:
        {
            
//{-0.0004967530258,19.99939919,-1.785099983,-0.4637269974,-30.90069962,-6.272890091,0.001483049942,-0.0004943400272,3.880409978*10^(-19)}
            
            u_r        = -0.0004943130258;
            sigma_r    = 0.0;
            sigma_t    = -30.90069962;
            sigma_z    = -6.272890091;
            p_r        = 20.0;
            q_r        = -1.785099983;
            eps_p_r    = 0.001483049942;
            eps_p_t    = -0.0004943400272;
            eps_p_z    = 0.0;
        }
            break;
        case 30:
        {
            
//{-0.0004957389901,19.99909973,-1.628520012,-0.4737080038,-30.91860008,-6.278460026,0.001444009948,-0.0004813270061,4.340599762*10^(-19)}
            
            u_r        = -0.0004938589901;
            sigma_r    = 0.0;
            sigma_t    = -30.91860008;
            sigma_z    = -6.278460026;
            p_r        = 20.0;
            q_r        = -1.628520012;
            eps_p_r    = 0.001444009948;
            eps_p_t    = -0.0004813270061;
            eps_p_z    = 0.0;
        }
        break;
        default:
        {
            DebugStop();
        }
            break;
    }

    /// Initial data at re
    std::vector<REAL> y_0;
    TPZTensor<REAL> sigma,eps, eps_p;
    sigma.Zero();
    
    y_0.push_back(u_r);
    y_0.push_back(sigma_r);
    y_0.push_back(p_r);
    y_0.push_back(q_r);
    
    sigma.XX() = sigma_r;
    sigma.YY() = sigma_t;
    sigma.ZZ() = sigma_z;
    
    eps_p.XX() = eps_p_r;
    eps_p.YY() = eps_p_t;
    eps_p.ZZ() = eps_p_z;
    
    /// Configuring the elastoplastic integrator
    TPMRSMemory default_memory;
    
    {
        if (sim_data->MaterialProps().size()!=1) {
            DebugStop();
        }
        
        int index = 0;
        
        std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters> chunk = sim_data->MaterialProps()[index];
        
        TPMRSUndrainedParameters udrained_parameters(std::get<0>(chunk));
        TPMRSPoroMechParameters poro_parameters(std::get<1>(chunk));
        TPMRSKappaParameters kappa_parameters(std::get<3>(chunk));
        std::vector<REAL> undrained_pars = udrained_parameters.GetParameters();
        std::vector<REAL> poroperm_pars  = poro_parameters.GetParameters();
        std::vector<REAL> kappa_param    = kappa_parameters.GetParameters();
        
        REAL phi_0   = undrained_pars[2];
        REAL kappa_0 = undrained_pars[3];
        REAL E       = poroperm_pars[0];
        REAL nu      = poroperm_pars[1];
        REAL Ks      = poroperm_pars[2];
        REAL c_f     = poroperm_pars[3];
        REAL eta     = poroperm_pars[4];
        REAL rho     = poroperm_pars[5];
        REAL Kdr     = E/(3.0*(1.0-2.0*nu));
//        REAL alpha   = 1.0 - (Kdr/Ks);
    
        TPZElasticResponse ER;
        ER.SetEngineeringData(E, nu);
        
        // Plastic corrector
        TPMRSPlasticityParameters plasticity_parameters(std::get<4>(chunk));
        std::vector<REAL> p_pars = plasticity_parameters.GetParameters();
        
        if (p_pars.size() == 0) {
            // Elastic material
            TPZElasticCriterion Elastic;
            Elastic.SetElasticResponse(ER);
            
            Elastic.ApplyLoad(sigma_0, eps);
            eps.ZZ() = 0.0;
            eps.Zero();

            
            default_memory.SetKs(Ks);
            default_memory.SetKdr(Kdr);
            default_memory.Setphi_0(phi_0);
            default_memory.Setphi_n(phi_0);
            default_memory.Setkappa_0(kappa_0);
            default_memory.Setkappa_n(kappa_0);
            default_memory.SetSigma_0(sigma_0);
            default_memory.SetSigma(sigma_0);
            default_memory.SetSigma_n(sigma_0);
            default_memory.Setp_0(p_0);
            default_memory.GetPlasticState_0().m_eps_t = eps;
            default_memory.GetPlasticState().m_eps_t = eps;
            default_memory.GetPlasticState_n().m_eps_t = eps;
            
            /// Configuring the solver
            TPMRSRKSolver<TPZElasticCriterion,TPMRSMemory> RKSolver;
            
            /// Configuring the permeability model
            if (kappa_param.size() == 0) {
                RKSolver.SetKappaParameters(kappa_parameters);
                
            }else{
                // Permeability model
                switch (kappa_parameters.GetModel()) {
                    case kappa_parameters.k_petunin: {
                        RKSolver.SetKappaParameters(kappa_parameters);
                    }
                        break;
                    case kappa_parameters.k_davies: {
                        RKSolver.SetKappaParameters(kappa_parameters);
                    }
                        break;
                    case kappa_parameters.k_costa: {
                        RKSolver.SetKappaParameters(kappa_parameters);
                    }
                        break;
                    case kappa_parameters.k_nelson: {
                        RKSolver.SetKappaParameters(kappa_parameters);
                    }
                        break;
                    case kappa_parameters.k_bayles: {
                        RKSolver.SetKappaParameters(kappa_parameters);
                    }
                        break;
                    default:{
                        std::cout << "Permeability model is not implemented. " << std::endl;
                        DebugStop();
                    }
                        break;
                }
            }
            
            RKSolver.SetPlasticIntegrator(Elastic);
            RKSolver.SetDefaultMemory(default_memory);
            RKSolver.SetInitialData(y_0);
            DebugStop();
            RKSolver.SetElastoPlasticInitialData(sigma,sigma,sigma);
            RKSolver.SetFluidData(eta, c_f, rho);
            RKSolver.SetDiscretization(rw, re, n_steps);
            RKSolver.SetGrainBulkModulus(Ks);
            RKSolver.SetFourthOrderApproximation();
            RKSolver.Synchronize();
            RKSolver.ExecuteRKApproximation();
            RKSolver.PrintRKApproximation(rk_file_data);
            RKSolver.PrintSecondaryVariables(rk_file_data);
            
        }else{
            // Elastoplastic material
            switch (plasticity_parameters.GetModel()) {
                case plasticity_parameters.ep_mc: {
                    // Mohr Coulomb data
                    REAL cohesion    = p_pars[0];
                    REAL phi         = (p_pars[1]*M_PI/180);
                    REAL psi         = phi;
                    
                    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
                    LEMC.SetElasticResponse(ER);
                    LEMC.fYC.SetUp(phi, psi, cohesion, ER);
                    
                    LEMC.ApplyLoad(sigma_0, eps);
                    eps.ZZ() = 0.0;
//                    eps.Zero();
                    
                    
                    /// computing elastoplastic state at initial point
                    TPZTensor<REAL> eps_e;
                    ER.ComputeStrain(sigma, eps_e);
                    
                    default_memory.SetKs(Ks);
                    default_memory.SetKdr(Kdr);
                    default_memory.Setphi_0(phi_0);
                    default_memory.Setphi_n(phi_0);
                    default_memory.Setkappa_0(kappa_0);
                    default_memory.Setkappa_n(kappa_0);
                    default_memory.SetSigma_0(sigma_0);
                    default_memory.SetSigma(sigma_0);
                    default_memory.SetSigma_n(sigma_0);
                    default_memory.Setp_0(p_0);
                    default_memory.GetPlasticState_0().m_eps_t = eps;
                    default_memory.GetPlasticState().m_eps_t = eps;
                    default_memory.GetPlasticState_n().m_eps_t = eps;
                    
                    /// Configuring the solver
                    TPMRSRKSolver<TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>,TPMRSMemory> RKSolver;
                    
                    /// Configuring the permeability model
                    if (kappa_param.size() == 0) {
                        RKSolver.SetKappaParameters(kappa_parameters);
                        
                    }else{
                        // Permeability models
                        switch (kappa_parameters.GetModel()) {
                            case kappa_parameters.k_petunin: {
                                RKSolver.SetKappaParameters(kappa_parameters);
                            }
                                break;
                            case kappa_parameters.k_davies: {
                                RKSolver.SetKappaParameters(kappa_parameters);
                            }
                                break;
                            case kappa_parameters.k_costa: {
                                RKSolver.SetKappaParameters(kappa_parameters);
                            }
                                break;
                            case kappa_parameters.k_nelson: {
                                RKSolver.SetKappaParameters(kappa_parameters);
                            }
                                break;
                            case kappa_parameters.k_bayles: {
                                RKSolver.SetKappaParameters(kappa_parameters);
                            }
                                break;
                            default:{
                                std::cout << "Permeability model is not implemented. " << std::endl;
                                DebugStop();
                            }
                                break;
                        }
                    }
                    
                    RKSolver.SetPlasticIntegrator(LEMC);
                    RKSolver.SetDefaultMemory(default_memory);
                    RKSolver.SetInitialData(y_0);
                    RKSolver.SetElastoPlasticInitialData(eps_e, eps_p, sigma);
                    RKSolver.SetFluidData(eta, c_f, rho);
//                    RKSolver.SetDefineDataAtRe();
                    RKSolver.SetDiscretization(rw, re, n_steps);
                    RKSolver.SetGrainBulkModulus(Ks);
                    RKSolver.SetFourthOrderApproximation();
                    RKSolver.Synchronize();
                    RKSolver.ExecuteRKApproximation();
                    RKSolver.PrintRKApproximation(rk_file_data);
                    RKSolver.PrintSecondaryVariables(rk_file_data);
                    
                }
                    break;
                case plasticity_parameters.ep_ds: {
                    // Dimaggio Sandler data
                    
                    STATE G   = E / (2.0 * (1.0 + nu));
                    STATE K   = E / (3.0 * (1.0 - 2 * nu));
                    
                    REAL A    = p_pars[0];
                    REAL B    = p_pars[1];
                    REAL C    = p_pars[2];
                    REAL D    = p_pars[3];
                    REAL R    = p_pars[4];
                    REAL W    = p_pars[5];
                    REAL X_0  = p_pars[6];
                    REAL phi = 0, psi = 1.0, N = 0;
                    
                    REAL Pc = X_0/3.0;
                    TPZTensor<REAL> sigma;
                    sigma.Zero();
                    
                    sigma.XX() = Pc;
                    sigma.YY() = Pc;
                    sigma.ZZ() = Pc;
                    
                    
                    TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
                    LEDS.SetElasticResponse(ER);
                    LEDS.fYC.SetUp(A, B, C, D, K, G, W, R, phi, N, psi);
                    
                    
                    // Initial damage data
                    REAL k_0;
                    LEDS.InitialDamage(sigma, k_0);
                    LEDS.fN.m_hardening = k_0;
                    LEDS.fYC.SetInitialDamage(k_0);
                    LEDS.fN.m_eps_t.Zero();
                    LEDS.fN.m_eps_p.Zero();
                    
                    LEDS.ApplyLoad(sigma, eps);
                    default_memory.SetKs(Ks);
                    default_memory.SetKdr(Kdr);
                    default_memory.Setphi_0(phi_0);
                    default_memory.Setphi_n(phi_0);
                    default_memory.Setkappa_0(kappa_0);
                    default_memory.Setkappa_n(kappa_0);
                    default_memory.SetSigma_0(sigma_0);
                    default_memory.SetSigma(sigma_0);
                    default_memory.SetSigma_n(sigma_0);
                    default_memory.GetPlasticState_0().m_eps_t = eps;
                    default_memory.GetPlasticState().m_eps_t = eps;
                    default_memory.GetPlasticState_n().m_eps_t = eps;
                    
                    
                    /// Configuring the solver
                    TPMRSRKSolver<TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>,TPMRSMemory> RKSolver;
                    
                    /// Configuring the permeability model
                    if (kappa_param.size() == 0) {
                        RKSolver.SetKappaParameters(kappa_parameters);
                        
                    }else{
                        // Permeability models
                        switch (kappa_parameters.GetModel()) {
                            case kappa_parameters.k_petunin: {
                                RKSolver.SetKappaParameters(kappa_parameters);
                            }
                                break;
                            case kappa_parameters.k_davies: {
                                RKSolver.SetKappaParameters(kappa_parameters);
                            }
                                break;
                            case kappa_parameters.k_costa: {
                                RKSolver.SetKappaParameters(kappa_parameters);
                            }
                                break;
                            case kappa_parameters.k_nelson: {
                                RKSolver.SetKappaParameters(kappa_parameters);
                            }
                                break;
                            case kappa_parameters.k_bayles: {
                                RKSolver.SetKappaParameters(kappa_parameters);
                            }
                                break;
                            default:{
                                std::cout << "Permeability model is not implemented. " << std::endl;
                                DebugStop();
                            }
                                break;
                        }
                    }
                    
                    RKSolver.SetPlasticIntegrator(LEDS);
                    RKSolver.SetDefaultMemory(default_memory);
                    RKSolver.SetInitialData(y_0);
                    RKSolver.SetFluidData(eta, c_f, rho);
                    RKSolver.SetDiscretization(rw, re, n_steps);
                    RKSolver.SetGrainBulkModulus(Ks);
                    RKSolver.SetFourthOrderApproximation();
                    RKSolver.Synchronize();
                    RKSolver.ExecuteRKApproximation();
                    RKSolver.PrintRKApproximation(rk_file_data);
                    RKSolver.PrintSecondaryVariables(rk_file_data);

                
                }
                    break;
                default:{
                    std::cout << "Material not implemented. " << std::endl;
                    DebugStop();
                }
                    break;
            }
        }
    }
    
    
}


TPMRSFullyCoupledAnalysis * CreateFCSolver(TPMRSSimulationData * sim_data)
{
    
    // Create multiphysisc mesh
    TPZManVector<TPZCompMesh * , 2 > mesh_vector(2);
    TPZCompMesh * cmesh_poro_perm_coupling = CMesh_FullyCoupled(mesh_vector,sim_data);
    
    // The initial condition is set up to zero for Deformation and Pore Pressure
    // Create and run the Transient analysis
    
    bool mustOptimizeBandwidth = true;
    int number_threads = sim_data->n_threads();
    
    TPMRSFullyCoupledAnalysis * fc_analysis = new TPMRSFullyCoupledAnalysis;
    fc_analysis->SetCompMesh(cmesh_poro_perm_coupling,mustOptimizeBandwidth);
    fc_analysis->SetSimulationData(sim_data);
    fc_analysis->SetMeshvec(mesh_vector);
    fc_analysis->AdjustVectors();
    
#ifdef USING_MKL
    TPZSpStructMatrix struct_mat(cmesh_poro_perm_coupling); // NonSymm Pardiso MKL flag
#else
    
    TPZSkylineNSymStructMatrix struct_mat(cmesh_poro_perm_coupling);
#endif
    
    TPZStepSolver<STATE> step;
    struct_mat.SetNumThreads(number_threads);
    step.SetDirect(ELU);
    fc_analysis->SetSolver(step);
    fc_analysis->SetStructuralMatrix(struct_mat);
    fc_analysis->ConfiguratePostProcessor();
    
    return fc_analysis;
}


TPMRSSegregatedAnalysis * CreateSFISolver(TPMRSSimulationData * sim_data){
    
    // The Geomechanics Simulator cmesh
    TPZCompMesh * cmesh_geomechanic = CMesh_Geomechanics(sim_data);
    
    // The Reservoir Simulator cmesh
    TPZCompMesh * cmesh_res = NULL;
    TPZManVector<TPZCompMesh * , 2 > mesh_vector(2);
    if (sim_data->Get_is_dual_formulation_Q()) {
        cmesh_res = CMesh_Mixed(mesh_vector,sim_data);
    }else{
        cmesh_res = CMesh_Primal(sim_data);
    }
 
    TPMRSSegregatedAnalysis * sfi_analysis = new TPMRSSegregatedAnalysis;
    sfi_analysis->ConfigurateAnalysis(ECholesky, ELU, sim_data, cmesh_geomechanic, cmesh_res, mesh_vector);
    return sfi_analysis;
}


// Method that makes use of
void RuningGeomechanics(TPMRSSimulationData * sim_data){
    
    TPZCompMesh * cmesh_geomechanic = CMesh_Geomechanics(sim_data);
    bool mustOptimizeBandwidth = false;
    TPMRSGeomechanicAnalysis * analysis = new TPMRSGeomechanicAnalysis;
    analysis->SetCompMesh(cmesh_geomechanic,mustOptimizeBandwidth);
    analysis->ConfigurateAnalysis(ECholesky, sim_data);
    
    int n_time_steps = sim_data->ReportingTimes().size();
    std::string file_name("Geomechanic.vtk");
    for (int it = 0; it < n_time_steps; it++) {
        analysis->ExecuteOneTimeStep(true);
        analysis->PostProcessTimeStep(file_name);
    }
}


TPZCompMesh * CMesh_Geomechanics(TPMRSSimulationData * sim_data){
    
    // Getting mesh dimension
    int dim = sim_data->Dimension();

    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = sim_data->MaterialIds();
    
    std::map<int, std::string>::iterator it_bc_id_to_type;
    std::map< std::string,std::pair<int,std::vector<std::string> > >::iterator  it_condition_type_to_index_value_names;
    std::map<int , TPMRSInterpolator >::iterator it_bc_id_to_values;
    
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        
        TPZMaterial  * material = ConfigurateAndInsertVolumetricMaterialsGeo(iregion,matid,sim_data,cmesh);
        
        // Inserting boundary conditions
        int n_bc = material_ids[iregion].second.first.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second.first [ibc];
            
            it_bc_id_to_type = sim_data->BCIdToConditionTypeGeomechanics().find(bc_id);
            it_bc_id_to_values = sim_data->BCIdToBCValuesGeomechanics().find(bc_id);
            it_condition_type_to_index_value_names = sim_data->ConditionTypeToBCIndexGeomechanics().find(it_bc_id_to_type->second);
            
            int bc_index = it_condition_type_to_index_value_names->second.first;
            int n_bc_values = it_bc_id_to_values->second.n_functions();
            TPZFMatrix<STATE> val1(0,0,0.), val2(n_bc_values,1,0.);
            for (int i = 0; i < n_bc_values; i++) {
                REAL value = 0.0; // Values are currently interpolated using time functions
                val2(i,0) = value;
            }
            
            TPZBndCondWithMem<TPMRSElastoPlasticMemory> * bc = new  TPZBndCondWithMem<TPMRSElastoPlasticMemory>(material, bc_id, bc_index, val1, val2);
            cmesh->InsertMaterialObject(bc);
        }
    }
    
    // Setting H1 approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(sim_data->ElasticityOrder());
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("Cmesh_Geomechanics.txt");
    cmesh->Print(out);
#endif
    return cmesh;
    
}

TPZMaterial * ConfigurateAndInsertVolumetricMaterialsGeo(int index, int matid, TPMRSSimulationData * sim_data, TPZCompMesh * cmesh){
    
    int dim = sim_data->Dimension();
    
    std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters> chunk =    sim_data->MaterialProps()[index];
    
    // Elastic predictor
    TPMRSPoroMechParameters poro_parameters(std::get<1>(chunk));
    std::vector<REAL> e_pars = poro_parameters.GetParameters();
    REAL E  = e_pars[0];
    REAL nu = e_pars[1];
    
    TPZElasticResponse ER;
    ER.SetEngineeringData(E, nu);

    // Plastic corrector
    TPMRSPlasticityParameters plasticity_parameters(std::get<4>(chunk));
    std::vector<REAL> p_pars = plasticity_parameters.GetParameters();
    
    if (p_pars.size() == 0) {
        // Elastic material
        TPZElasticCriterion Elastic;
        Elastic.SetElasticResponse(ER);
        
        TPMRSElastoPlastic<TPZElasticCriterion, TPMRSMemory> * material = new TPMRSElastoPlastic <TPZElasticCriterion, TPMRSMemory>(matid);
        material->SetDimension(dim);
        material->SetPlasticIntegrator(Elastic);
        material->SetSimulationData(sim_data);
        cmesh->InsertMaterialObject(material);
        return material;
        
    }else{
        // Elastoplastic material
        switch (plasticity_parameters.GetModel()) {
            case plasticity_parameters.ep_mc: {
                // Mohr Coulomb data
                REAL cohesion    = p_pars[0];
                REAL phi         = (p_pars[1]*M_PI/180);
                REAL psi         = phi;
                
                TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
                LEMC.SetElasticResponse(ER);
                LEMC.fYC.SetUp(phi, psi, cohesion, ER);
                
                TPMRSElastoPlastic <TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPMRSMemory> * material = new TPMRSElastoPlastic <TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPMRSMemory>(matid);
                material->SetDimension(dim);
                material->SetPlasticIntegrator(LEMC);
                
                material->SetSimulationData(sim_data);
                cmesh->InsertMaterialObject(material);
                return material;
            }
            break;
            case plasticity_parameters.ep_ds: {
                // Dimaggio Sandler data
                
                STATE G   = E / (2.0 * (1.0 + nu));
                STATE K   = E / (3.0 * (1.0 - 2 * nu));
                
                REAL A    = p_pars[0];
                REAL B    = p_pars[1];
                REAL C    = p_pars[2];
                REAL D    = p_pars[3];
                REAL R    = p_pars[4];
                REAL W    = p_pars[5];
                REAL X_0  = p_pars[6];
                REAL phi = 0, psi = 1.0, N = 0;
                
                REAL Pc = X_0/3.0;
                TPZTensor<REAL> sigma;
                sigma.Zero();
                
                sigma.XX() = Pc;
                sigma.YY() = Pc;
                sigma.ZZ() = Pc;
            
                
                TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
                LEDS.SetElasticResponse(ER);
                LEDS.fYC.SetUp(A, B, C, D, K, G, W, R, phi, N, psi);
                
                
                // Initial damage data
                REAL k_0;
                LEDS.InitialDamage(sigma, k_0);
                LEDS.fN.m_hardening = k_0;
                LEDS.fYC.SetInitialDamage(k_0);
                LEDS.fN.m_eps_t.Zero();
                LEDS.fN.m_eps_p.Zero();
                
                TPMRSElastoPlastic <TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>, TPMRSMemory> * material = new TPMRSElastoPlastic <TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>, TPMRSMemory>(matid);
                material->SetDimension(dim);
                material->SetPlasticIntegrator(LEDS);
                
                material->SetSimulationData(sim_data);
                cmesh->InsertMaterialObject(material);
                return material;
            }
                break;
            default:{
                TPZMaterial * material = NULL;
                std::cout << "Material not implemented. " << std::endl;
                DebugStop();
                return material;
            }
                break;
        }
    }
}

void RuningMonophasic(TPMRSSimulationData * sim_data){
    
#ifdef PZDEBUG
    sim_data->PrintGeometry();
#endif
    
    TPZManVector<TPZCompMesh * , 2 > mesh_vector(2);
    TPZCompMesh * cmesh_mixed = CMesh_Mixed(mesh_vector,sim_data);
    
    bool mustOptimizeBandwidth = false;
    TPMRSMonoPhasicAnalysis * analysis = new TPMRSMonoPhasicAnalysis;
    analysis->SetCompMesh(cmesh_mixed,mustOptimizeBandwidth);
    analysis->ConfigurateAnalysis(ELU, mesh_vector, sim_data);
    
    int n_time_steps = sim_data->ReportingTimes().size();
    std::string file_name("Monophasic.vtk");
    for (int it = 0; it < n_time_steps; it++) {
        analysis->ExecuteOneTimeStep();
        analysis->PostProcessTimeStep(file_name);
    }    
}


TPZCompMesh * CMesh_Flux(TPMRSSimulationData * sim_data)
{
    
    // Getting mesh dimension
    int dim = sim_data->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    int nstate = 1;
    TPZVec<STATE> sol;
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = sim_data->MaterialIds();
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
        cmesh->InsertMaterialObject(material);
        
        // Inserting boundary conditions
        int dirichlet = 0;
        TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
        int n_bc = material_ids[iregion].second.second.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second.second [ibc];
            TPZMaterial * bc = material->CreateBC(material, bc_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc);
        }
    }
    // Setting Hdiv approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(sim_data->DiffusionOrder());
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshFlux.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

TPZCompMesh * CMesh_PorePressure_disc(TPMRSSimulationData * sim_data)
{
    
    // Getting mesh dimension
    int dim = sim_data->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    int nstate = 1;
    TPZVec<STATE> sol;
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12> material_ids = sim_data->MaterialIds();
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
        cmesh->InsertMaterialObject(material);
    
    }
    
    // Setting Discontinuous approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(sim_data->DiffusionOrder());
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    
#ifndef USING_MKL2
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
#endif
    
#ifdef PZDEBUG
        std::ofstream out("CmeshPorePressure_Disc.txt");
        cmesh->Print(out);
#endif
    
    return cmesh;
}

TPZCompMesh * CMesh_Mixed(TPZManVector<TPZCompMesh * , 2 > & mesh_vector, TPMRSSimulationData * sim_data){
    
    mesh_vector[0] = CMesh_Flux(sim_data);
    mesh_vector[1] = CMesh_PorePressure_disc(sim_data);
    
    int dim = sim_data->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    
    std::map<int, std::string>::iterator it_bc_id_to_type;
    std::map< std::string,std::pair<int,std::vector<std::string> > >::iterator  it_condition_type_to_index_value_names;
    std::map<int , TPMRSInterpolator >::iterator it_bc_id_to_values;
    
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = sim_data->MaterialIds();
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        TPZMaterial * material = ConfigurateAndInsertVolumetricMaterialsRes(true, iregion, matid, sim_data, cmesh);
        
        // Inserting boundary conditions
        int n_bc = material_ids[iregion].second.second.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second.second [ibc];
            
            it_bc_id_to_type = sim_data->BCIdToConditionTypeReservoir().find(bc_id);
            it_bc_id_to_values = sim_data->BCIdToBCValuesReservoir().find(bc_id);
            it_condition_type_to_index_value_names = sim_data->ConditionTypeToBCIndexReservoir().find(it_bc_id_to_type->second);
            
            int bc_index = it_condition_type_to_index_value_names->second.first;
//            int n_bc_values = it_bc_id_to_values->second.n_functions();
            TPZFMatrix<STATE> val1(0,0,0.), val2(1,1,0.);
            REAL value = 0; // Values are currently interpolated using time functions
            val2(0,0) = value;
            
            // Memory not used but is required for coherence with geomechanics integration points
            TPZBndCondWithMem<TPMRSMonoPhasicMemory> * bc = new  TPZBndCondWithMem<TPMRSMonoPhasicMemory>(material, bc_id, bc_index, val1, val2);
            cmesh->InsertMaterialObject(bc);
        }
    }
    
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    cmesh->ApproxSpace().CreateWithMemory(true);
    cmesh->AutoBuild();
    
    TPZBuildMultiphysicsMesh::AddElements(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(mesh_vector, cmesh);
    
    int nel_res = cmesh->NElements();
    for (long el = 0; el < nel_res; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->PrepareIntPtIndices();
    }
    
#ifdef PZDEBUG
    std::ofstream out("CmeshReservoir.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

TPZCompMesh * CMesh_Primal(TPMRSSimulationData * sim_data){
    
    int dim = sim_data->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    
    std::map<int, std::string>::iterator it_bc_id_to_type;
    std::map< std::string,std::pair<int,std::vector<std::string> > >::iterator  it_condition_type_to_index_value_names;
    std::map<int , TPMRSInterpolator >::iterator it_bc_id_to_values;
    
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = sim_data->MaterialIds();
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        TPZMaterial * material = ConfigurateAndInsertVolumetricMaterialsRes(false, iregion, matid, sim_data, cmesh);
        
        // Inserting boundary conditions
        int n_bc = material_ids[iregion].second.second.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second.second [ibc];
            
            it_bc_id_to_type = sim_data->BCIdToConditionTypeReservoir().find(bc_id);
            it_bc_id_to_values = sim_data->BCIdToBCValuesReservoir().find(bc_id);
            it_condition_type_to_index_value_names = sim_data->ConditionTypeToBCIndexReservoir().find(it_bc_id_to_type->second);
            
            int bc_index = it_condition_type_to_index_value_names->second.first;
            TPZFMatrix<STATE> val1(0,0,0.), val2(1,1,0.);
            REAL value = 0; // Values are currently interpolated using time functions
            val2(0,0) = value;
            
            // Memory not used but is required for coherence with geomechanics integration points
            TPZBndCondWithMem<TPMRSMonoPhasicMemory> * bc = new  TPZBndCondWithMem<TPMRSMonoPhasicMemory>(material, bc_id, bc_index, val1, val2);
            cmesh->InsertMaterialObject(bc);
        }
    }
    
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(sim_data->DiffusionOrder());
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();

#ifdef PZDEBUG
    std::ofstream out("CmeshReservoirCG.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

TPZMaterial * ConfigurateAndInsertVolumetricMaterialsRes(bool IsMixedQ, int index, int matid, TPMRSSimulationData * sim_data, TPZCompMesh * cmesh){
    
    /// Scale factor
    REAL s  = sim_data->scale_factor_val();
    int dim = sim_data->Dimension();
    
    bool  is_crank_nicolson_Q = sim_data->Get_is_crank_nicolson_Q();
    
    std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters> chunk =    sim_data->MaterialProps()[index];
    
    // Reservoir parameters
    TPMRSPoroMechParameters poro_parameters(std::get<1>(chunk));
    std::vector<REAL> res_pars = poro_parameters.GetParameters();
    REAL c_f   = res_pars[3];
    REAL eta   = res_pars[4];
    REAL rho_0 = res_pars[5];


    if (IsMixedQ) {
        TPMRSMonoPhasic<TPMRSMemory> * material = new TPMRSMonoPhasic<TPMRSMemory>(matid,dim);
        material->SetSimulationData(sim_data);
        material->SetFluidProperties(rho_0, eta, c_f);
        material->SetScaleFactor(s);
        material->SetPorosityParameters(std::get<2>(chunk));
        material->SetPermeabilityParameters(std::get<3>(chunk));
        cmesh->InsertMaterialObject(material);
        if (is_crank_nicolson_Q) {
            material->SetCrank_Nicolson();
        }
        return material;
    }else{
        TPMRSMonoPhasicCG<TPMRSMemory> * material = new TPMRSMonoPhasicCG<TPMRSMemory>(matid,dim);
        material->SetSimulationData(sim_data);
        material->SetFluidProperties(rho_0, eta, c_f);
        material->SetScaleFactor(s);
        material->SetPorosityParameters(std::get<2>(chunk));
        material->SetPermeabilityParameters(std::get<3>(chunk));
        if (is_crank_nicolson_Q) {
            material->SetCrank_Nicolson();
        }
        cmesh->InsertMaterialObject(material);
        return material;
    }

}


TPZCompMesh * CMesh_Deformation(TPMRSSimulationData * sim_data){
    
    // Getting mesh dimension
    int dim = sim_data->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    int nstate = dim;
    TPZVec<STATE> sol;
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = sim_data->MaterialIds();
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
        cmesh->InsertMaterialObject(material);
        
        // Inserting boundary conditions
        int dirichlet = 0;
        TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
        int n_bc = material_ids[iregion].second.first.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second.first [ibc];
            TPZMaterial * bc = material->CreateBC(material, bc_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc);
        }
    }
    // Setting H1 approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(sim_data->ElasticityOrder());
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    //    std::ofstream out("CmeshDeformation.txt");
    //    cmesh->Print(out);
#endif
    
    return cmesh;
}


TPZCompMesh * CMesh_PorePressure(TPMRSSimulationData * sim_data)
{
    
    // Getting mesh dimension
    int dim = sim_data->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    int nstate = 1;
    TPZVec<STATE> sol;
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = sim_data->MaterialIds();
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
        cmesh->InsertMaterialObject(material);
        
        // Inserting boundary conditions
        int dirichlet = 0;
        TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
        int n_bc = material_ids[iregion].second.first.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second.first [ibc];
            TPZMaterial * bc = material->CreateBC(material, bc_id, dirichlet, val1, val2);
            
            cmesh->InsertMaterialObject(bc);
        }
    }
    // Setting H1 approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(sim_data->DiffusionOrder());
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    //    std::ofstream out("CmeshPorePressure.txt");
    //    cmesh->Print(out);
#endif
    
    return cmesh;
}


TPZCompMesh * CMesh_FullyCoupled(TPZManVector<TPZCompMesh * , 2 > & mesh_vector, TPMRSSimulationData * sim_data){
    
    mesh_vector[0] = CMesh_Deformation(sim_data);
    mesh_vector[1] = CMesh_PorePressure(sim_data);
    
    // Plane strain assumption
    int dim = sim_data->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    
    std::map<int, std::string>::iterator it_bc_id_to_type;
    std::map< std::string,std::pair<int,std::vector<std::string> > >::iterator  it_condition_type_to_index_value_names;
    std::map<int , TPMRSInterpolator >::iterator it_bc_id_to_values;
    
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = sim_data->MaterialIds();
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;

        TPZMaterial  * material = ConfigurateAndInsertVolumetricMaterialsFC(iregion, matid, sim_data, cmesh);
        
        // Inserting boundary conditions
        int n_bc = material_ids[iregion].second.first.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second.first [ibc];
            
            it_bc_id_to_type = sim_data->BCIdToConditionTypeGeomechanics().find(bc_id);
            it_bc_id_to_values = sim_data->BCIdToBCValuesGeomechanics().find(bc_id);
            it_condition_type_to_index_value_names = sim_data->ConditionTypeToBCIndexGeomechanics().find(it_bc_id_to_type->second);
            
            int bc_index = it_condition_type_to_index_value_names->second.first;
            int n_bc_values = it_bc_id_to_values->second.n_functions();
            TPZFMatrix<STATE> val1(0,0,0.), val2(n_bc_values,1,0.);
            for (int i = 0; i < n_bc_values; i++) {
                REAL value = 0.0; // Values are currently interpolated using time functions
                val2(i,0) = value;
            }
            
            TPZBndCondWithMem<TPMRSMemory> * bc = new  TPZBndCondWithMem<TPMRSMemory>(material, bc_id, bc_index, val1, val2);
            cmesh->InsertMaterialObject(bc);
        }
    }
    
    // Setting up multiphysics functions
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    // Transfer to multiphysic mesh
    TPZBuildMultiphysicsMesh::AddElements(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(mesh_vector, cmesh);
    
    long nel = cmesh->NElements();
    TPZVec<long> indices;
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel)
        {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
    
#ifdef PZDEBUG
    std::ofstream out("Cmesh_FullyCoupled.txt");
    cmesh->Print(out);
#endif
    cmesh->InitializeBlock();

    return cmesh;
    
}

TPZMaterial * ConfigurateAndInsertVolumetricMaterialsFC(int index, int matid, TPMRSSimulationData * sim_data, TPZCompMesh * cmesh){
    
    REAL s  = sim_data->scale_factor_val();
    int dim = sim_data->Dimension();
    
    std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters> chunk =    sim_data->MaterialProps()[index];
    
    // Elastic predictor
    TPMRSPoroMechParameters poro_parameters(std::get<1>(chunk));
    std::vector<REAL> pars = poro_parameters.GetParameters();
    REAL E  = pars[0];
    REAL nu = pars[1];
    
    TPZElasticResponse ER;
    ER.SetEngineeringData(E, nu);
    
    bool  is_crank_nicolson_Q = sim_data->Get_is_crank_nicolson_Q();
    // Reservoir parameters
    REAL c_f   = pars[3];
    REAL eta   = pars[4];
    REAL rho_0 = pars[5];
    
    // Plastic corrector
    TPMRSPlasticityParameters plasticity_parameters(std::get<4>(chunk));
    std::vector<REAL> p_pars = plasticity_parameters.GetParameters();
    
    if (p_pars.size() == 0) {
        // Elastic material
        TPZElasticCriterion Elastic;
        Elastic.SetElasticResponse(ER);
        
        TPMRSPoroElastoPlastic<TPZElasticCriterion, TPMRSMemory> * material = new TPMRSPoroElastoPlastic<TPZElasticCriterion, TPMRSMemory>(matid);
        material->SetDimension(dim);
        material->SetPlasticIntegrator(Elastic);
        material->SetSimulationData(sim_data);
        
        material->SetFluidProperties(rho_0, eta, c_f);
        material->SetScaleFactor(s);
        material->SetPorosityParameters(std::get<2>(chunk));
        material->SetPermeabilityParameters(std::get<3>(chunk));
        if (is_crank_nicolson_Q) {
            material->SetCrank_Nicolson();
        }
        
        cmesh->InsertMaterialObject(material);
        return material;
        
    }else{
        // Elastoplastic material
        switch (plasticity_parameters.GetModel()) {
            case plasticity_parameters.ep_mc: {
                // Mohr Coulomb data
                REAL cohesion    = p_pars[0];
                REAL phi         = (p_pars[1]*M_PI/180);
                REAL psi         = phi;
                
                TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
                LEMC.SetElasticResponse(ER);
                LEMC.fYC.SetUp(phi, psi, cohesion, ER);
                
                TPMRSPoroElastoPlastic <TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPMRSMemory> * material = new TPMRSPoroElastoPlastic <TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPMRSMemory>(matid);
                material->SetDimension(dim);
                material->SetPlasticIntegrator(LEMC);
                
                material->SetSimulationData(sim_data);
                
                material->SetFluidProperties(rho_0, eta, c_f);
                material->SetScaleFactor(s);
                material->SetPorosityParameters(std::get<2>(chunk));
                material->SetPermeabilityParameters(std::get<3>(chunk));
                if (is_crank_nicolson_Q) {
                    material->SetCrank_Nicolson();
                }
                
                cmesh->InsertMaterialObject(material);
                return material;
            }
                break;
            case plasticity_parameters.ep_ds: {
                // Dimaggio Sandler data
                
                STATE G   = E / (2.0 * (1.0 + nu));
                STATE K   = E / (3.0 * (1.0 - 2 * nu));
                
                REAL A    = p_pars[0];
                REAL B    = p_pars[1];
                REAL C    = p_pars[2];
                REAL D    = p_pars[3];
                REAL R    = p_pars[4];
                REAL W    = p_pars[5];
                REAL X_0  = p_pars[6];
                REAL phi = 0, psi = 1.0, N = 0;
                
                REAL Pc = X_0/3.0;
                TPZTensor<REAL> sigma;
                sigma.Zero();
                
                sigma.XX() = Pc;
                sigma.YY() = Pc;
                sigma.ZZ() = Pc;
                
                
                TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
                LEDS.SetElasticResponse(ER);
                LEDS.fYC.SetUp(A, B, C, D, K, G, W, R, phi, N, psi);
                
                
                // Initial damage data
                REAL k_0;
                LEDS.InitialDamage(sigma, k_0);
                LEDS.fN.m_hardening = k_0;
                LEDS.fYC.SetInitialDamage(k_0);
                LEDS.fN.m_eps_t.Zero();
                LEDS.fN.m_eps_p.Zero();
                
                TPMRSPoroElastoPlastic <TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>, TPMRSMemory> * material = new TPMRSPoroElastoPlastic <TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>, TPMRSMemory>(matid);
                material->SetDimension(dim);
                material->SetPlasticIntegrator(LEDS);
                
                material->SetSimulationData(sim_data);
                
                material->SetFluidProperties(rho_0, eta, c_f);
                material->SetScaleFactor(s);
                material->SetPorosityParameters(std::get<2>(chunk));
                material->SetPermeabilityParameters(std::get<3>(chunk));
                if (is_crank_nicolson_Q) {
                    material->SetCrank_Nicolson();
                }
                
                cmesh->InsertMaterialObject(material);
                return material;

            }
                break;
            default:{
                TPZMaterial * material = NULL;
                std::cout << "Material not implemented. " << std::endl;
                DebugStop();
                return material;
            }
                break;
        }
    }
}


void LEDSPorosityReductionPlot()
{
    
    // Getting Elastic Matrix, MC Mohr Coloumb PV
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    
    // Experimental data
    std::string dirname = PZSOURCEDIR;
    std::string file_name;
    file_name = dirname + "/Projects/PoropermCoupling/exp_data/fullstrain.txt";
    int64_t n_data = 2500; //2575
    TPZFMatrix<REAL> data = Read_Duplet(n_data, file_name);
//    data.Print(std::cout);
    
    // DS Dimaggio Sandler PV
    TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    
    /// Input data for shear enhanced compaction:

    REAL E  = 43457.2; // MPa * 1.025
    REAL nu = 0.357983; // MPa
    
    STATE G = E / (2. * (1. + nu));
    STATE K = E / (3. * (1. - 2 * nu));
    REAL CA      = 400;
    REAL CB      = 0.001;
    REAL CC      = 200;
    REAL CD      = 0.001;
    REAL CR      = 1.0;
    REAL CW      = 0.035;
    REAL phi = 0, psi = 1., N = 0;
    
    REAL Pc = -100.0;

    ER.SetEngineeringData(E, nu);
    
    // Mohr Coulomb data
    REAL mc_cohesion    = 25.0;
    REAL mc_phi         = 10.5*M_PI/180;
    REAL mc_psi         = mc_phi; // because MS do not understand
    
    LEMC.SetElasticResponse(ER);
    LEMC.fYC.SetUp(mc_phi, mc_psi, mc_cohesion, ER);
    
    LEDS.SetElasticResponse(ER);
    LEDS.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
    
    TPZTensor<REAL> epsilon_t,sigma,sigma_target;
    sigma.Zero();
    epsilon_t.Zero();
    
    
    sigma.XX() = Pc;
    sigma.YY() = Pc;
    sigma.ZZ() = Pc;
    
    // Initial damage data
    REAL k_0;
    LEDS.InitialDamage(sigma, k_0);
    LEDS.fN.m_hardening = k_0;
    
//    TPZFNMatrix<36,STATE> De_c(6,6,0.0),De(6,6,0.0),De_inv;
//    ER.SetUp(E, nu);
//    LEMC.SetElasticResponse(ER);
//    LEMC.fYC.SetUp(phi, psi, c, ER);
//    LEMC.ApplyStrainComputeSigma(epsilon_t, sigma, &De_c);
//    LEMC.ApplyStrainComputeSigma(epsilon_t, sigma, &De);
//    De_c.Inverse(De_inv, ECholesky);
    
    // For a given stress
//    STATE sigma_c = -137.9/3; // MPa
//    sigma_target.Zero();
    
    epsilon_t.Zero();
    
    TPZFNMatrix<2575,STATE> LEDS_epsilon_stress(n_data,2);
    for (int64_t id = 0; id < n_data; id++) {
        
//        sigma_target.XX() = sigma_c + data(id,1)*0.005;
//        sigma_target.YY() = sigma_c;
//        sigma_target.ZZ() = sigma_c;
//        Apply_Stress(LEDS, De, De_inv, sigma_target, epsilon_t);
        
        // For a given strain
        epsilon_t.XX() = data(id,0);
        epsilon_t.YY() = data(id,1);
        epsilon_t.ZZ() = data(id,1);
        
//        LEDS.ApplyStrainComputeSigma(epsilon_t, sigma_target);
        LEMC.ApplyStrainComputeSigma(epsilon_t, sigma_target);
        
        
//        LEDS_epsilon_stress(id,0) = epsilon_t.XX() - epsilon_t.YY();
//        LEDS_epsilon_stress(id,1) = sigma_target.XX() - sigma_target.YY();
        
//        LEDS_epsilon_stress(id,0) = epsilon_t.XX();
//        LEDS_epsilon_stress(id,1) = sigma_target.XX();
        
        LEDS_epsilon_stress(id,0) = epsilon_t.XX()+epsilon_t.YY()+epsilon_t.ZZ();
        LEDS_epsilon_stress(id,1) = (1/3.)*(sigma_target.XX()+sigma_target.YY()+sigma_target.ZZ());
    }

    LEDS_epsilon_stress.Print("data = ", std::cout,EMathematicaInput);
}


void Apply_Stress(TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> &LEDS, TPZFMatrix<REAL> &De, TPZFMatrix<REAL> &De_inv, TPZTensor<REAL> &sigma_target, TPZTensor<REAL> &epsilon)
{
    
    TPZPlasticState<STATE> plastic_state;
    plastic_state = LEDS.fN;
    
    STATE tol = 1.0e-2;
    int64_t n_iter = 50;
    STATE res_val;
    
    TPZFNMatrix<6,REAL> eps,eps_e_0(6,1,0.0),sigma_0(6,1,0.0);
    TPZFNMatrix<6,REAL> res(6,1,0.0);
    TPZTensor<REAL> sigma_x,dsigma,sigma_res,epsilon_e,epsilon_p,depsilon,depsilon_acum;

    depsilon_acum.Zero();
    sigma_res = sigma_target - sigma_x;
    epsilon   = plastic_state.m_eps_t;
    epsilon_p = plastic_state.m_eps_p;
    epsilon_e = epsilon - epsilon_p;
    
    epsilon_e.CopyTo(eps_e_0);
    De.Multiply(eps_e_0, sigma_0);
    sigma_x.CopyFrom(sigma_0);
    
    for (int64_t i = 0; i < n_iter; i++)
    {
        sigma_res.CopyTo(res);
        De_inv.Multiply(res, eps);
        depsilon.CopyFrom(eps);
//        depsilon_acum +=  depsilon;
//        epsilon = epsilon_e + epsilon_p + depsilon_acum;
        epsilon += depsilon;
        LEDS.ApplyStrainComputeSigma(epsilon,sigma_x);

        sigma_res = sigma_target - sigma_x;
        res_val = sigma_res.Norm();
        bool stop_criterion_Q = res_val < tol;
        if (stop_criterion_Q)
        { // Important Step
            std::cout << "Converged with i " << i << std::endl;
            break;
        }
        LEDS.fN = plastic_state;
        int aka = 0;
        aka = 1;
    }
}

TPZFMatrix<REAL> Read_Duplet(int n_data, std::string file)
{
    TPZFMatrix<REAL> data(n_data,2);
    std::ifstream in(file.c_str());
    int count = 0;
    while(in)
    {
        in >> data(count,0);
        in >> data(count,1);
        
        count++;
        if (count == n_data)
        {
            break;
        }
    }
    return data;
}
