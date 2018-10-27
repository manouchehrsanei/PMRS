
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

// Geometry
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"
#include "tpzhierarquicalgrid.h"

// Computational mesh
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzintel.h"
#include "pzlog.h"

// Materials
#include "pzl2projection.h"
#include "pzbndcond.h"
#include "TPZBndCondWithMem.h"

// Monophasic
#include "TPMRSMonoPhasic_impl.h"
#include "TPMRSMonoPhasicCG_impl.h"
#include "TPMRSMonoPhasicAnalysis.h"
#include "TPMRSMonoPhasicMemory.h"


// Geomechanics
#include "TPMRSElastoPlastic_impl.h"
#include "TPMRSGeomechanicAnalysis.h"
#include "TPMRSMemoryPoroElast.h"
#include "TPMRSElastoPlasticMemory.h"

// Segregated solver
#include "TPMRSMemory.h"
#include "TPMRSSegregatedAnalysis.h"


// Elasticity
#include "TPZElasticCriterion.h"

// Plasticity
#include "TPZPlasticStepPV.h"
#include "TPZSandlerExtended.h"
#include "TPZYCMohrCoulombPV.h"


#include "pzporoelastoplasticmem.h"
#include "TPZPlasticStepPV.h"
#include "TPZSandlerDimaggio.h"

// Analysis
#include "pzpostprocanalysis.h"
#include "pzanalysis.h"
#include "TPMRSAnalysis.h"


// Matrix
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzfstrmatrix.h"

// Simulation data structure
#include "TPMRSSimulationData.h"
#include "TPMRSUndrainedParameters.h"
#include "TPMRSPoroMechParameters.h"
#include "TPMRSPhiParameters.h"
#include "TPMRSKappaParameters.h"
#include "TPMRSPlasticityParameters.h"

// ElastoPlastic Materials
#include "TPMRSCouplPoroElast.h"
#include "TPMRSCouplPoroPlast.h"


#ifdef LOG4CXX
static LoggerPtr log_data(Logger::getLogger("pz.PMRS"));
#endif

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


/// PoroElastic Full Coupling
// L2 mesh for Deformation
TPZCompMesh * CMesh_Deformation(TPMRSSimulationData * sim_data);
// L2 mesh for Pore Pressure
TPZCompMesh * CMesh_PorePressure(TPMRSSimulationData * sim_data);
// Multiphysics Coupling
TPZCompMesh * CMesh_FullCoupling(TPZManVector<TPZCompMesh * , 2 > & mesh_vector, TPMRSSimulationData * sim_data);


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

// Method that makes the poroelastic full coupling
void RuningFullCoupling(TPMRSSimulationData * sim_data);


/// Restructuring implementation of Reservoir Geomechanics Simulator

// Method that makes use of TPMRSMonophasic for parabolic solutions
void RuningMonophasic(TPMRSSimulationData * sim_data);

// Method that makes use of TPMRSElastoPlastic
void RuningGeomechanics(TPMRSSimulationData * sim_data);

// Method that makes use of TPMRSMonophasic and TPMRSElastoPlastic with a common memory
void RuningSegregatedSolver(TPMRSSimulationData * sim_data);




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



int main(int argc, char *argv[])
{
    
    
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
            cout << "Usage: " << argv[0] << " Myinputfile.xml " << endl;
            cout <<    "Program stop: not xml file found \n" << endl;
            DebugStop();
        }
        
        if (argc == 2)
        {
            cout << "Control File used : " << argv[1] << "\n" << endl;
            simulation_file        = argv[1];
        }
    }
    
    // Simulation data to be configurated
    
    TPMRSSimulationData * sim_data = new TPMRSSimulationData;
    sim_data->ReadSimulationFile(simulation_file);
    
//    RuningGeomechanics(sim_data);
//    RuningMonophasic(sim_data);
    
#ifdef USING_BOOST
    boost::posix_time::ptime int_case_t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    RuningSegregatedSolver(sim_data);
    
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



void RuningFullCoupling(TPMRSSimulationData * sim_data)
{
  
#ifdef PZDEBUG
    //    sim_data->PrintGeometry();
#endif
    
    // Create multiphysisc mesh
    TPZManVector<TPZCompMesh * , 2 > mesh_vector(2);
    TPZCompMesh * cmesh_poro_perm_coupling = CMesh_FullCoupling(mesh_vector,sim_data);
    
    
    // The initial condition is set up to zero for Deformation and Pore Pressure
    // Create and run the Transient analysis
    
    bool mustOptimizeBandwidth = true;
    int number_threads = sim_data->n_threads();
    
#ifdef PZDEBUG
    number_threads = 0;
#endif
    
    TPMRSAnalysis * time_analysis = new TPMRSAnalysis;
    time_analysis->SetCompMesh(cmesh_poro_perm_coupling,mustOptimizeBandwidth);
    time_analysis->SetSimulationData(sim_data);
    time_analysis->SetMeshvec(mesh_vector);
    time_analysis->AdjustVectors();
    
#ifdef USING_Pardiso
    //    TPZSymetricSpStructMatrix struct_mat(cmesh_poro_perm_coupling); // Symm Pardiso MKL flag
    TPZSpStructMatrix struct_mat(cmesh_poro_perm_coupling); // NonSymm Pardiso MKL flag
#else
    
    TPZSkylineNSymStructMatrix struct_mat(cmesh_poro_perm_coupling);
    //    TPZSkylineStructMatrix struct_mat(cmesh_poro_perm_coupling);
    //    TPZFStructMatrix struct_mat(cmesh_poro_perm_coupling);
    
    //    TPZParFrontStructMatrix<TPZFrontSym<STATE> > struct_mat(cmesh_poro_perm_coupling);
    //    struct_mat.SetDecomposeType(ELDLt);
    
#endif
    
    
    TPZStepSolver<STATE> step;
    struct_mat.SetNumThreads(number_threads);
    step.SetDirect(ELU);
    time_analysis->SetSolver(step);
    time_analysis->SetStructuralMatrix(struct_mat);
    
    TPZVec<REAL> x(3);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    std::string file_ss_name("AxialStrainVsDiffStress.nb");
    std::string file_sp_name("VolumStrainVsPorosity.nb");
    std::string file_sk_name("VolumStrainVsPermeability.nb");
    std::string file_spex_name("VolumStrainVsPorePressure.nb");
    
    // Run Transient analysis
    time_analysis->Run_Evolution(x);
    //    time_analysis->PlotStrainStress(file_ss_name);
    //    time_analysis->PlotStrainPorosity(file_sp_name);
    //    time_analysis->PlotStrainPermeability(file_sk_name);
    //    time_analysis->PlotStrainPressure(file_spex_name);
    std::cout << " Execution finished" << std::endl;

}


void RuningSegregatedSolver(TPMRSSimulationData * sim_data){
    
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
 
    TPMRSSegregatedAnalysis * segregated_analysis = new TPMRSSegregatedAnalysis;
    if (sim_data->Get_is_dual_formulation_Q()) {
        segregated_analysis->ConfigurateAnalysis(ECholesky, ELU, sim_data, cmesh_geomechanic, cmesh_res, mesh_vector);
    }else{
        segregated_analysis->ConfigurateAnalysis(ECholesky, ELU, sim_data, cmesh_geomechanic, cmesh_res, mesh_vector);
    }

    segregated_analysis->ConfigurateBConditions(true);
    segregated_analysis->ExecuteStaticSolution();
    segregated_analysis->ConfigurateBConditions(false);
    segregated_analysis->ExecuteTimeEvolution();
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
        analysis->ExecuteOneTimeStep();
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
    std::map<int , std::vector<REAL> >::iterator it_bc_id_to_values;
    
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
            int n_bc_values = it_bc_id_to_values->second.size();
            TPZFMatrix<STATE> val1(0,0,0.), val2(dim,1,0.);
            for (int i = 0; i < n_bc_values; i++) {
                REAL value = it_bc_id_to_values->second[i];
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
    REAL E = e_pars[0];
    REAL nu = e_pars[1];
    
    // Updating bulk modulus for porosity model
    std::get<2>(sim_data->MaterialProps()[index]).SetBulkModulus(E, nu);
    
    TPZElasticResponse ER;
    ER.SetUp(E, nu);

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
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
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
    std::map<int , std::vector<REAL> >::iterator it_bc_id_to_values;
    
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
            
            it_bc_id_to_type = sim_data->BCIdToConditionTypeReservoirs().find(bc_id);
            it_bc_id_to_values = sim_data->BCIdToBCValuesReservoirs().find(bc_id);
            it_condition_type_to_index_value_names = sim_data->ConditionTypeToBCIndexReservoirs().find(it_bc_id_to_type->second);
            
            int bc_index = it_condition_type_to_index_value_names->second.first;
            int n_bc_values = it_bc_id_to_values->second.size();
            TPZFMatrix<STATE> val1(0,0,0.), val2(1,1,0.);
            REAL value = it_bc_id_to_values->second[n_bc_values-1];
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
    std::map<int , std::vector<REAL> >::iterator it_bc_id_to_values;
    
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
            
            it_bc_id_to_type = sim_data->BCIdToConditionTypeReservoirs().find(bc_id);
            it_bc_id_to_values = sim_data->BCIdToBCValuesReservoirs().find(bc_id);
            it_condition_type_to_index_value_names = sim_data->ConditionTypeToBCIndexReservoirs().find(it_bc_id_to_type->second);
            
            int bc_index = it_condition_type_to_index_value_names->second.first;
            int n_bc_values = it_bc_id_to_values->second.size();
            TPZFMatrix<STATE> val1(0,0,0.), val2(1,1,0.);
            REAL value = it_bc_id_to_values->second[n_bc_values-1];
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
    
    REAL s = 1.0e-6;
    REAL c = 0.0;
    int dim = sim_data->Dimension();
    
    std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters> chunk =    sim_data->MaterialProps()[index];
    
    // Reservoir parameters
    TPMRSPoroMechParameters poro_parameters(std::get<1>(chunk));
    std::vector<REAL> res_pars = poro_parameters.GetParameters();
    REAL eta = res_pars[4];
    REAL rho_0 = res_pars[5];

    if (IsMixedQ) {
        TPMRSMonoPhasic<TPMRSMemory> * material = new TPMRSMonoPhasic<TPMRSMemory>(matid,dim);
        material->SetSimulationData(sim_data);
        material->SetFluidProperties(rho_0, eta, c);
        material->SetScaleFactor(s);
        material->SetPorosityParameters(std::get<2>(chunk));
        material->SetPermeabilityParameters(std::get<3>(chunk));
        cmesh->InsertMaterialObject(material);
        return material;
    }else{
        TPMRSMonoPhasicCG<TPMRSMemory> * material = new TPMRSMonoPhasicCG<TPMRSMemory>(matid,dim);
        material->SetSimulationData(sim_data);
        material->SetFluidProperties(rho_0, eta, c);
        material->SetScaleFactor(s);
        material->SetPorosityParameters(std::get<2>(chunk));
        material->SetPermeabilityParameters(std::get<3>(chunk));
        cmesh->InsertMaterialObject(material);
        return material;
    }

}

TPZCompMesh * CMesh_FullCoupling(TPZManVector<TPZCompMesh * , 2 > & mesh_vector, TPMRSSimulationData * sim_data){
    
    mesh_vector[0] = CMesh_Deformation(sim_data);
    mesh_vector[1] = CMesh_PorePressure(sim_data);
    
    // Plane strain assumption
    int dim = sim_data->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    
    std::map<int, std::string>::iterator it_bc_id_to_type;
    std::map< std::string,std::pair<int,std::vector<std::string> > >::iterator  it_condition_type_to_index_value_names;
    std::map<int , std::vector<REAL> >::iterator it_bc_id_to_values;
    
    
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = sim_data->MaterialIds();
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;

        TPMRSCouplPoroElast * material = new TPMRSCouplPoroElast(matid,dim);

        DebugStop();
//        int kmodel = 0;
//        REAL Ey_r = sim_data->Get_young();
//        REAL nu_r = sim_data->Get_nu();
//        REAL porosity = sim_data->Get_porosity();
//        REAL k = sim_data->Get_k();
//        REAL alpha_r = sim_data->Get_alpha();
//        REAL Se = sim_data->Get_Se();
//        REAL eta = sim_data->Get_eta();
//        REAL rho_f = sim_data->Get_rho_f();
//        REAL rho_s = sim_data->Get_rho_s();
//        REAL MC_coh = sim_data->Get_mc_coh();
//        REAL MC_phi = sim_data->Get_mc_phi();
//        REAL MC_psi = sim_data->Get_mc_psi();
//        
//        material->SetSimulationData(sim_data);
//        material->SetPorolasticParametersEngineer(Ey_r, nu_r);
//        material->SetBiotParameters(alpha_r, Se);
//        material->SetParameters(k, porosity, eta);
//        material->SetKModel(kmodel);
//        material->SetDensityFluidRock(rho_f, rho_s);
//        material->SetMohrCoulombParameters(MC_coh, MC_phi, MC_psi);
//        cmesh->InsertMaterialObject(material);
        
        DebugStop();
        // Inserting boundary conditions
        int n_bc = material_ids[iregion].second.first.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second.first [ibc];
            
            it_bc_id_to_type = sim_data->BCIdToConditionTypeReservoirs().find(bc_id);
            
            it_bc_id_to_values = sim_data->BCIdToBCValuesReservoirs().find(bc_id);
            it_condition_type_to_index_value_names = sim_data->ConditionTypeToBCIndexReservoirs().find(it_bc_id_to_type->second);
            
            int bc_index = it_condition_type_to_index_value_names->second.first;
            int n_bc_values = it_bc_id_to_values->second.size();
            TPZFMatrix<STATE> val1(0,0,0.), val2(4,1,0.);
            for (int i = 0; i < n_bc_values; i++) {
                REAL value = it_bc_id_to_values->second[i];
                val2(i,0) = value;
            }
            
            TPZMaterial * bc = material->CreateBC(material, bc_id, bc_index, val1, val2);
            
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
//    std::ofstream out("PorePermCoupling.txt");
//    cmesh->Print(out);
#endif
    cmesh->InitializeBlock();

    return cmesh;
    
}


void LEDSPorosityReductionPlot()
{
    
    // Getting Elastic Matrix
    // MC Mohr Coloumb PV
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
    
    
    /**
     * Input data for shear enhanced compaction:
     *
     */
    
    REAL E =43457.2; // MPa * 1.025
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

    ER.SetUp(E, nu);
    
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
