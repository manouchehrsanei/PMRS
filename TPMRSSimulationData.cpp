//
//  TPMRSSimulationData.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#include "TPMRSSimulationData.h"

TPMRSSimulationData::TPMRSSimulationData()
{
    m_dt                                  = 0.0;
    m_n_steps                             = 0  ;
    m_reporting_times.resize(0);
    m_time                                = 0.0;
    m_n_iteraions                         =   0;
    m_epsilon_res                         = 1.0;
    m_epsilon_cor                         = 1.0;
    m_n_fss_iterations                    =   0;
    m_n_enf_fss_iterations                =   0;
    m_max_plastic_strain                  = 0.0;
    m_n_threads                           =   0;
    m_scale_factor                        = 0.0;
    m_is_dual_formulation_Q               = true;
    m_transfer_current_to_last_solution_Q = false;
    m_h_level                             =   0;
    m_elasticity_order                    =   0;
    m_diffusion_order                     =   0;
    m_dimesion                            =   0;
    m_geometry_file                       =  "";
    m_geometry                            = NULL;
    m_vtk_file                            =  "";
    m_vtk_resolution                      =   0;
    m_n_outputs_geo                       =   0;
    m_n_outputs_res                       =   0;
    m_s_names_res.Resize(0);
    m_s_names_geo.Resize(0);
    m_v_names_res.Resize(0);
    m_v_names_geo.Resize(0);
    m_t_names_geo.Resize(0);
    m_g.Resize(0);
    m_n_regions                          =  0;
    m_mat_ids.Resize(0);
    m_mat_props.Resize(0);
    m_is_initial_state_Q                 = false;
    m_is_current_state_Q                 = false;
    m_is_crank_nicolson_Q                = false;
    m_n_sub_step_level                   = 0;
    m_must_use_sub_stepping_Q            = false;
    
}

TPMRSSimulationData::~TPMRSSimulationData()
{
    
}

void TPMRSSimulationData::ReadSimulationFile(char *simulation_file)
{
    
    TiXmlDocument document(simulation_file);
    bool file_ok_Q = false;
    
    file_ok_Q = document.LoadFile();
    if (file_ok_Q)
    {
        std::cout << "This Xml is ok! -> " << simulation_file << std::endl;
    }
    else
    {
        std::cout << "Failed to load file \n"      << simulation_file << std::endl;
        std::cout << "Check the given path or your xml structure. \n" << std::endl;
        DebugStop();
    }
    
    
    /// TiXmlElement dummy object
    TiXmlElement * container;
    const char * char_container;
    TiXmlHandle doc_handler( & document );
    
    /// Begin:: Geometry description
    container = doc_handler.FirstChild("CaseData").FirstChild("Mesh").FirstChild("MeshFile").ToElement();
    const char * geometry_file = container->Attribute("mesh_file");
    m_geometry_file = geometry_file;
    
    this->ReadGeometry();
    int dimension = m_geometry->Dimension();
    m_dimesion = dimension;
    /// End:: Geometry description
    
    
    /// Begin:: Time controls
    container = doc_handler.FirstChild("CaseData").FirstChild("TimeControls").FirstChild("StepSize").ToElement();
    char_container = container->Attribute("dt");
    REAL dt = std::atof(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("TimeControls").FirstChild("StepNumber").ToElement();
    char_container = container->Attribute("n_time_steps");
    int n_stpes = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("TimeControls").FirstChild("CrankNicolsonQ").ToElement();
    char_container = container->Attribute("useQ");
    bool is_crank_nicolson_Q = std::atoi(char_container);
    
    SetTimeControls(n_stpes,dt,is_crank_nicolson_Q);
    /// End:: Time controls
    
    
    /// Begin:: Newton method controls
    container = doc_handler.FirstChild("CaseData").FirstChild("NewtonControls").FirstChild("Iterations").ToElement();
    char_container = container->Attribute("n_iterations");
    int n_iterations = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("NewtonControls").FirstChild("Residue").ToElement();
    char_container = container->Attribute("res_tolerance");
    REAL epsilon_res = std::atof(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("NewtonControls").FirstChild("Correction").ToElement();
    char_container = container->Attribute("cor_tolerance");
    REAL epsilon_cor = std::atof(char_container);
    
    SetNumericControls(n_iterations,epsilon_res,epsilon_cor);
    /// End:: Newton method controls
    
    
    /// Begin:: Fixed Stress Split Scheme
    container = doc_handler.FirstChild("CaseData").FirstChild("FixedStressSplit").FirstChild("FssIterations").ToElement();
    char_container = container->Attribute("n_max_fss_iterations");
    int n_fss_iterations = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("FixedStressSplit").FirstChild("EnfFssIterations").ToElement();
    char_container = container->Attribute("n_enforce_fss_iterations");
    int n_enf_fss_iterations = std::atoi(char_container);
    
    SetFixedStressSplitSchemes(n_fss_iterations,n_enf_fss_iterations);
    /// End:: Fixed Stress Split Scheme

    
    /// Begin:: Substeps controls
    container = doc_handler.FirstChild("CaseData").FirstChild("SubSteps").FirstChild("MaxPlasticNorm").ToElement();
    char_container = container->Attribute("max_plastic_norm_value");
    REAL max_plastic_norm_val = std::atof(char_container);
    m_max_plastic_strain = max_plastic_norm_val;
    /// End:: Substeps controls
    
    
    /// Begin:: Parallel controls
    container = doc_handler.FirstChild("CaseData").FirstChild("ParallelControls").FirstChild("Numthreads").ToElement();
    char_container = container->Attribute("n_threads");
    int n_threads = std::atoi(char_container);
    m_n_threads = n_threads;
    /// End:: Parallel controls
    
    
    /// Begin:: Scale factor controls
    container = doc_handler.FirstChild("CaseData").FirstChild("ScaleFactor").FirstChild("Valscalefactor").ToElement();
    char_container = container->Attribute("scalfac_value");
    REAL scale_factor_value = std::atof(char_container);
    m_scale_factor = scale_factor_value;
    /// End:: Scale factor controls
    
    
    /// Begin:: Finite elements
    container = doc_handler.FirstChild("CaseData").FirstChild("FEM").FirstChild("MixedFormulationQ").ToElement();
    char_container = container->Attribute("useQ");
    bool is_mixed_formulation_Q = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("FEM").FirstChild("HRefine").ToElement();
    char_container = container->Attribute("h_level");
    int h_level = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("FEM").FirstChild("Elasticity").ToElement();
    char_container = container->Attribute("p_order");
    int elasticity_order = std::atoi(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("FEM").FirstChild("Diffusion").ToElement();
    char_container = container->Attribute("p_order");
    int diffusion_order = std::atoi(char_container);
    
    m_is_dual_formulation_Q     = is_mixed_formulation_Q;
    m_h_level                   = h_level;
    m_elasticity_order          = elasticity_order;
    m_diffusion_order           = diffusion_order;
    /// End:: Finite elements
    
    
    /// Begin:: Outputs
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("OutputFolder").ToElement();
    char_container = container->Attribute("name");
    system(char_container);
    std::cout << "Creating directory using : "<< char_container << "\n" << std::endl;
    
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("LogFolder").ToElement();
    char_container = container->Attribute("name");
    system(char_container);
    std::cout << "Creating directory using : "<< char_container << "\n" << std::endl;
    
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("PostProcessing").ToElement();
    const char * vtk_file = container->Attribute("vtk_file");
    m_vtk_file = vtk_file;
    
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("PostProcessing").ToElement();
    char_container = container->Attribute("n_divisions");
    int vtk_resolution = std::atoi(char_container);
    m_vtk_resolution = vtk_resolution;

    
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("PostProcessing").ToElement();
    char_container = container->Attribute("n_outputs_geo");
    int n_outputs_geo = std::atoi(char_container);
    m_n_outputs_geo = n_outputs_geo;

    
    container = doc_handler.FirstChild("CaseData").FirstChild("OutputControls").FirstChild("PostProcessing").ToElement();
    char_container = container->Attribute("n_outputs_res");
    int n_outputs_res = std::atoi(char_container);
    m_n_outputs_res = n_outputs_res;
    
    

    int is_r = 0;
    int iv_r = 0;
    container = doc_handler.FirstChild( "CaseData" ).FirstChild( "OutputControls" ).FirstChild( "OutputControlsRes" ).FirstChild( "Var" ).ToElement();
    for( ; container; container=container->NextSiblingElement())
    {
        char_container = container->Attribute("v_name");
        if (char_container) {
            iv_r++;
            m_v_names_res.Resize(iv_r,"");
            m_v_names_res[iv_r-1] = char_container;
        }
        char_container = container->Attribute("s_name");
        if (char_container) {
            is_r++;
            m_s_names_res.Resize(is_r,"");
            m_s_names_res[is_r-1] = char_container;
        }
    }
    
    if (is_r+iv_r!=n_outputs_res) {
        std::cout << "OutputControlsRes has a different number of variables that you request with n_outputs_res." << std::endl;
        DebugStop();
    }

    int is_g = 0;
    int iv_g = 0;
    int it_g = 0;
    container = doc_handler.FirstChild( "CaseData" ).FirstChild( "OutputControls" ).FirstChild( "OutputControlsGeo" ).FirstChild( "Var" ).ToElement();
    for( ; container; container=container->NextSiblingElement())
    {
        char_container = container->Attribute("t_name");
        if (char_container) {
            it_g++;
            m_t_names_geo.Resize(it_g,"");
            m_t_names_geo[it_g-1] = char_container;
        }
        char_container = container->Attribute("v_name");
        if (char_container) {
            iv_g++;
            m_v_names_geo.Resize(iv_g,"");
            m_v_names_geo[iv_g-1] = char_container;
        }
        char_container = container->Attribute("s_name");
        if (char_container) {
            is_g++;
            m_s_names_geo.Resize(is_g,"");
            m_s_names_geo[is_g-1] = char_container;
        }
    }
    
    if (is_g+iv_g+it_g!=n_outputs_geo) {
        std::cout << "OutputControlsGeo has a different number of variables that you request with n_outputs_geo." << std::endl;
        DebugStop();
    }
    
     /// End:: Outputs
    
    
    /// Begin:: Physics
    container = doc_handler.FirstChild("CaseData").FirstChild("Physics").FirstChild("GravityConstant").ToElement();
    char_container = container->Attribute("gravity");
    REAL g_c = std::atof(char_container);
    
    container = doc_handler.FirstChild("CaseData").FirstChild("Physics").FirstChild("GravityDirection").ToElement();
    char_container = container->Attribute("x_direction");
    REAL x_dir = std::atof(char_container);
    char_container = container->Attribute("y_direction");
    REAL y_dir = std::atof(char_container);
    m_g.Resize(dimension,0.0);
    m_g[0] = g_c * x_dir;
    m_g[1] = g_c * y_dir;
    
    if (dimension == 3)
    {
        char_container = container->Attribute("z_direction");
        REAL z_dir = std::atof(char_container);
        m_g[2] = g_c * z_dir;
    }
    
    /// End:: Physics
    
    
    /// Begin:: Regions and materials parameters
    container = doc_handler.FirstChild("CaseData").FirstChild("ReservoirRegions").FirstChild("RegionNumber").ToElement();
    char_container = container->Attribute("n_regions");
    int n_regions = std::atoi(char_container);
    
    m_n_regions = n_regions;
    m_mat_ids.Resize(n_regions);
    m_mat_props.Resize(n_regions);

    int iregion = 0;
    TiXmlElement * sub_container;
    container = doc_handler.FirstChild( "CaseData" ).FirstChild( "RegionsDefinition" ).FirstChild("RegionData").ToElement();
    for( ; container; container=container->NextSiblingElement())
    {
        
        char_container = container->Attribute("mat_id");
        int mat_id = std::atoi(char_container);
        char_container = container->Attribute("n_boundaries_geo");
        int n_boundaries_geo = std::atoi(char_container);
        char_container = container->Attribute("n_boundaries_res");
        int n_boundaries_res = std::atoi(char_container);
        
        m_mat_ids[iregion].first = mat_id;
        m_mat_ids[iregion].second.first.Resize(n_boundaries_geo);
        m_mat_ids[iregion].second.second.Resize(n_boundaries_res);
        
        sub_container = container->FirstChild("GeoBoundaries")->FirstChild("Boundary")->ToElement();
        int iboundary_g = 0;
        for( ; sub_container; sub_container=sub_container->NextSiblingElement())
        {
            char_container = sub_container->Attribute("bc_id");
            int bc_id = std::atoi(char_container);
            m_mat_ids[iregion].second.first[iboundary_g] = bc_id;
            iboundary_g++;
        }
        
        sub_container = container->FirstChild("ResBoundaries")->FirstChild("Boundary")->ToElement();
        int iboundary_r = 0;
        for( ; sub_container; sub_container=sub_container->NextSiblingElement())
        {
            char_container = sub_container->Attribute("bc_id");
            int bc_id = std::atoi(char_container);
            m_mat_ids[iregion].second.second[iboundary_r] = bc_id;
            iboundary_r++;
        }
        
        if ((iboundary_g != n_boundaries_geo) || (iboundary_r != n_boundaries_res)) {
        std::cout << "RegionData has a different number of boundary conditions that you request with n_boundaries_geo or n_boundaries_res." << std::endl;
            DebugStop();
        }
        
        std::vector<REAL> pars;
        TPMRSUndrainedParameters  u_pars;
        TPMRSPoroMechParameters   poro_pars;
        TPMRSPhiParameters        phi_pars;
        TPMRSKappaParameters      kappa_pars;
        TPMRSPlasticityParameters plasticity_pars;
        
        std::tuple <TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters, TPMRSPlasticityParameters> chunk;
        
        sub_container = container->FirstChild("InitialPoroMechParameters")->ToElement();
        pars.resize(4);
        char_container = sub_container->Attribute("Eyoung_u");
        if (!char_container) {
            std::cout << "Please provide Eyoung_u." << std::endl;
            DebugStop();
        }
        pars[0] = std::atof(char_container);
        char_container = sub_container->Attribute("nu_u");
        if (!char_container) {
            std::cout << "Please provide nu_u." << std::endl;
            DebugStop();
        }
        pars[1] = std::atof(char_container);
        char_container = sub_container->Attribute("phi_0");
        if (!char_container) {
            std::cout << "Please provide phi_0." << std::endl;
            DebugStop();
        }
        pars[2] = std::atof(char_container);
        char_container = sub_container->Attribute("kappa_0");
        if (!char_container) {
            std::cout << "Please provide kappa_0." << std::endl;
            DebugStop();
        }
        pars[3] = std::atof(char_container);
        u_pars.SetParameters(pars);
        
        sub_container = container->FirstChild("PoroMechParameters")->ToElement();
        pars.resize(7);
        char_container = sub_container->Attribute("Eyoung");
        if (!char_container) {
            std::cout << "Please provide Eyoung." << std::endl;
            DebugStop();
        }
        pars[0] = std::atof(char_container);
        
        char_container = sub_container->Attribute("nu");
        if (!char_container) {
            std::cout << "Please provide nu." << std::endl;
            DebugStop();
        }
        pars[1] = std::atof(char_container);
        
        char_container = sub_container->Attribute("alpha");
        if (!char_container) {
            std::cout << "Please provide alpha." << std::endl;
            DebugStop();
        }
        pars[2] = std::atof(char_container);
        
        char_container = sub_container->Attribute("c_f");
        if (!char_container) {
            std::cout << "Please provide c_f." << std::endl;
            DebugStop();
        }
        pars[3] = std::atof(char_container);
        
        char_container = sub_container->Attribute("eta");
        if (!char_container) {
            std::cout << "Please provide eta." << std::endl;
            DebugStop();
        }
        pars[4] = std::atof(char_container);
        
        char_container = sub_container->Attribute("rho_f");
        if (!char_container) {
            std::cout << "Please provide rho_f." << std::endl;
            DebugStop();
        }
        pars[5] = std::atof(char_container);
        
        char_container = sub_container->Attribute("rho_s");
        if (!char_container) {
            std::cout << "Please provide rho_s." << std::endl;
            DebugStop();
        }
        pars[6] = std::atof(char_container);
        poro_pars.SetParameters(pars);
        
        /// Porosity model data
        sub_container = container->FirstChild("PhiParameters")->ToElement();
        pars.resize(0);
        char_container = sub_container->Attribute("phi_model");
        if (!char_container) {
            std::cout << "Please provide phi_model name." << std::endl;
            DebugStop();
        }
        std::string phi_model(char_container);
        char_container = sub_container->Attribute("n_parameters");
        if (!char_container) {
            std::cout << "Please provide n_parameters >= 0." << std::endl;
            DebugStop();
        }
        int n_phi_pars = std::atoi(char_container);
        pars.resize(n_phi_pars);
        
        /// Reading porosity parameters by case
        if(n_phi_pars!=0){
        }
        phi_pars.SetModel(phi_model);
        phi_pars.SetParameters(pars);
        
        /// Permeability model data
        sub_container = container->FirstChild("KappaParameters")->ToElement();
        pars.resize(0);
        char_container = sub_container->Attribute("kappa_model");
        if (!char_container) {
            std::cout << "Please provide kappa_model name." << std::endl;
            DebugStop();
        }
        std::string kappa_model(char_container);
        kappa_pars.SetModel(kappa_model);
        char_container = sub_container->Attribute("n_parameters");
        if (!char_container) {
            std::cout << "Please provide n_parameters >= 0." << std::endl;
            DebugStop();
        }
        int n_kappa_pars = std::atoi(char_container);
        pars.resize(n_kappa_pars);
        
        /// Reading permeability parameters by case
        if(n_kappa_pars!=0){
            switch (kappa_pars.GetModel()) {
                case kappa_pars.k_constant:
                    {
                        /// nothing to say
                    }
                    break;
                case kappa_pars.k_petunin:
                {
                    char_container = sub_container->Attribute("a");
                    if (!char_container) {
                        std::cout << "PetuninModel::Please provide a parameter" << std::endl;
                        DebugStop();
                    }
                    REAL a = std::atof(char_container);
                    pars[0] = a;
                }
                    break;
                case kappa_pars.k_davies:
                {
                    char_container = sub_container->Attribute("c");
                    if (!char_container) {
                        std::cout << "DaviesModel::Please provide c parameter" << std::endl;
                        DebugStop();
                    }
                    REAL c = std::atof(char_container);
                    pars[0] = c;
                }
                    break;
                case kappa_pars.k_costa:
                {
                    char_container = sub_container->Attribute("a");
                    char_container = sub_container->Attribute("c");
                    if (!char_container) {
                        std::cout << "CostaModel::Please provide a $ c parameters" << std::endl;
                        DebugStop();
                    }
                    REAL a = std::atof(char_container);
                    REAL c = std::atof(char_container);
                    pars[0] = a;
                    pars[1] = c;
                }
                    break;
                case kappa_pars.k_nelson:
                {
                    char_container = sub_container->Attribute("a");
                    char_container = sub_container->Attribute("c");
                    if (!char_container) {
                        std::cout << "NelsonModel::Please provide a $ c parameters" << std::endl;
                        DebugStop();
                    }
                    REAL a = std::atof(char_container);
                    REAL c = std::atof(char_container);
                    pars[0] = a;
                    pars[1] = c;
                }
                    
                    break;
                case kappa_pars.k_bayles:
                {
                    char_container = sub_container->Attribute("a");
                    char_container = sub_container->Attribute("c");
                    if (!char_container) {
                        std::cout << "BaylesModel::Please provide a $ c parameters" << std::endl;
                        DebugStop();
                    }
                    REAL a = std::atof(char_container);
                    REAL c = std::atof(char_container);
                    pars[0] = a;
                    pars[1] = c;
                }
                    
                    break;
                    
                default:
                {
                    std::cout << "KappaModel::Please provide a permeability model example {kappa_model, n_parameters, a}." << std::endl;
                    DebugStop();
                }
                    break;
            }
            
        }
        kappa_pars.SetParameters(pars);
        
        /// Plasticity model data
        sub_container = container->FirstChild("PlasticityParameters")->ToElement();
        pars.resize(0);
        char_container = sub_container->Attribute("plasticity_model");
        if (!char_container) {
            std::cout << "Please provide plasticity_model name." << std::endl;
            DebugStop();
        }
        std::string plasticity_model(char_container);
        char_container = sub_container->Attribute("n_parameters");
        if (!char_container) {
            std::cout << "Please provide n_parameters >= 0." << std::endl;
            DebugStop();
        }
        int n_plas_pars = std::atoi(char_container);
        pars.resize(n_plas_pars);
        plasticity_pars.SetModel(plasticity_model);
        
        /// Reading plasticity parameters by case
        if(n_plas_pars!=0){
            switch (plasticity_pars.GetModel()) {
                case plasticity_pars.ep_mc:
                {
                    if (n_plas_pars!=2) {
                        std::cout << "MCModel::Please set n_parameters = 2" << std::endl;
                    }
                    char_container = sub_container->Attribute("cohesion");
                    if (!char_container) {
                        std::cout << "MCModel::Please provide the cohesion" << std::endl;
                        DebugStop();
                    }
                    REAL cohesion = std::atof(char_container);
                    pars[0] = cohesion;
                    
                    char_container = sub_container->Attribute("friction");
                    if (!char_container) {
                        std::cout << "MCModel::Please provide the friction" << std::endl;
                        DebugStop();
                    }
                    REAL friction = std::atof(char_container);
                    pars[1] = friction;
                }
                    break;
                case plasticity_pars.ep_ds:
                {
                    if (n_plas_pars!=7) {
                        std::cout << "DSModel::Please set n_parameters = 7" << std::endl;
                    }
                    char_container = sub_container->Attribute("a");
                    if (!char_container) {
                        std::cout << "DSModel::Please provide a" << std::endl;
                        DebugStop();
                    }
                    REAL a = std::atof(char_container);
                    pars[0] = a;
                    
                    char_container = sub_container->Attribute("b");
                    if (!char_container) {
                        std::cout << "DSModel::Please provide b" << std::endl;
                        DebugStop();
                    }
                    REAL b = std::atof(char_container);
                    pars[1] = b;
                    
                    char_container = sub_container->Attribute("c");
                    if (!char_container) {
                        std::cout << "DSModel::Please provide c" << std::endl;
                        DebugStop();
                    }
                    REAL c = std::atof(char_container);
                    pars[2] = c;
                    
                    char_container = sub_container->Attribute("d");
                    if (!char_container) {
                        std::cout << "DSModel::Please provide d" << std::endl;
                        DebugStop();
                    }
                    REAL d = std::atof(char_container);
                    pars[3] = d;
                    
                    char_container = sub_container->Attribute("r");
                    if (!char_container) {
                        std::cout << "DSModel::Please provide r" << std::endl;
                        DebugStop();
                    }
                    REAL r = std::atof(char_container);
                    pars[4] = r;
                    
                    char_container = sub_container->Attribute("w");
                    if (!char_container) {
                        std::cout << "DSModel::Please provide w" << std::endl;
                        DebugStop();
                    }
                    REAL w = std::atof(char_container);
                    pars[5] = w;
                    
                    char_container = sub_container->Attribute("x0");
                    if (!char_container) {
                        std::cout << "DSModel::Please provide x0" << std::endl;
                        DebugStop();
                    }
                    REAL x0 = std::atof(char_container);
                    pars[6] = x0;
                }
                    break;
                default:
                {
                    std::cout << "PlasticityModel::Please provide a plasticity model example {plasticity_model, n_parameters, cohesion, friction}." << std::endl;
                    DebugStop();
                }
                    break;
            }
            plasticity_pars.SetParameters(pars);
            
        }
        
        /// Assigning values to tuple using make_tuple()
        chunk = make_tuple(u_pars, poro_pars, phi_pars, kappa_pars, plasticity_pars);
        m_mat_props[iregion] = chunk;
        
        iregion++;
    }
    /// End:: Regions and materials parameters
    

    /// Begin:: BC for Geomechanic Simulator
    this->LoadBoundaryConditionsGeomechanics();
    
    std::pair<int, std::string > bc_id_to_type_chunk_geo;
    std::pair<int , TPMRSInterpolator > bc_id_to_values_chunk_geo;
    std::map< std::string,std::pair<int,std::vector<std::string> > >::iterator chunk_geo;
    
    /// Intitial condition
    TiXmlElement * sub_container_geo;
    container = doc_handler.FirstChild("CaseData").FirstChild("BCInitialGeomechanics").FirstChild("InitialGeomechanic").ToElement();
    for( ; container; container=container->NextSiblingElement())
    {

        char_container = container->Attribute("bc_id");
        int bc_id = std::atoi(char_container);

        char_container = container->Attribute("type");
        std::string condition(char_container);
        chunk_geo = m_condition_type_to_index_value_names_geo.find(condition);

        bool bc_condition_not_available_Q = chunk_geo == m_condition_type_to_index_value_names_geo.end();
        if (bc_condition_not_available_Q)
        {
            std::cout << " Geometry dimension =  " << dimension << std::endl;
            std::cout << " The boundary " << condition << " are not available " << std::endl;
            std::cout << " Please review your boundary conditions for Geomechanic Module. " << std::endl;
            DebugStop();
        }

        /// Association bc type with numerical values
        bc_id_to_values_chunk_geo.first = bc_id;
        bc_id_to_values_chunk_geo.second.Clear();
        
        char_container = container->Attribute("n_data");
        if (!char_container)
        {
            std::cout << " The boundary " << condition << "  needs the number of data to be interporlated. " << std::endl;
            std::cout << " Please review your boundary conditions for Geomechanic Module. " << std::endl;
            DebugStop();
        }
        int n_time_data = std::atoi(char_container);
        std::vector<std::pair<REAL, std::vector<REAL>>> bc_points_geo(n_time_data);
        
        /// Association bc type with numerical values
        bc_id_to_values_chunk_geo.first = bc_id;
        bc_id_to_values_chunk_geo.second.Clear();
        
        sub_container_geo = container->FirstChild("Data")->ToElement();
        int itime_data = 0;
        for( ; sub_container_geo; sub_container_geo=sub_container_geo->NextSiblingElement())
        {
            char_container = sub_container_geo->Attribute("t");
            if (!char_container)
            {
                std::cout << " the boundary " << condition << "  needs the time value t ." << std::endl;
                std::cout << " Please review your boundary conditions for Geomechanic Module. " << std::endl;
                DebugStop();
            }
            REAL time_value = std::atoi(char_container);
            
            int n_data = chunk_geo->second.second.size();
            bc_points_geo[itime_data].first = time_value;
            bc_points_geo[itime_data].second.resize(n_data);
            for (int i = 0; i < n_data; i++)
            {
                char_container = sub_container_geo->Attribute(chunk_geo->second.second[i].c_str());
                if (!char_container)
                {
                    std::cout << " the boundary " << condition << "  needs the value " << chunk_geo->second.second[i] << std::endl;
                    std::cout << " Please review your boundary conditions for Geomechanic Module. " << std::endl;
                    DebugStop();
                }
                REAL bc_value = std::atof(char_container);
                bc_points_geo[itime_data].second[i] = bc_value;
            }
            itime_data++;
        }
        
        if (itime_data != n_time_data) {
            std::cout << " the boundary " << condition << " is not properly uploaded. " << std::endl;
            std::cout << " It was provided n_time_data = " << n_time_data << std::endl;
            std::cout << " It was uploaded this number of time data = " << itime_data << std::endl;
            std::cout << " Please review your boundary conditions for Geomechanic Module. " << std::endl;
            DebugStop();
        }
        TPMRSInterpolator bc_interpolator_geo(bc_points_geo);
        bc_id_to_values_chunk_geo.second = bc_interpolator_geo;
        m_bc_id_to_values_geo_initial.insert(bc_id_to_values_chunk_geo);

        /// Association bc identifier with bc type
        bc_id_to_type_chunk_geo.first = bc_id;
        bc_id_to_type_chunk_geo.second = condition;
        m_bc_id_to_type_geo_initial.insert(bc_id_to_type_chunk_geo);

    }
    
    /// Recurrent condition
    container = doc_handler.FirstChild("CaseData").FirstChild("BCGeomechanics").FirstChild("Geomechanic").ToElement();
    for( ; container; container=container->NextSiblingElement())
    {
        
        char_container = container->Attribute("bc_id");
        int bc_id = std::atoi(char_container);
        
        char_container = container->Attribute("type");
        std::string condition(char_container);
        chunk_geo = m_condition_type_to_index_value_names_geo.find(condition);
        
        bool bc_condition_not_available_Q = chunk_geo == m_condition_type_to_index_value_names_geo.end();
        if (bc_condition_not_available_Q)
        {
            std::cout << " Geometry dimension =  " << dimension << std::endl;
            std::cout << " The boundary " << condition << " are not available " << std::endl;
            std::cout << " Please review your boundary conditions for Geomechanic Module. " << std::endl;
            DebugStop();
        }
        
        char_container = container->Attribute("n_data");
        if (!char_container)
        {
            std::cout << " The boundary " << condition << "  needs the number of data to be interporlated. " << std::endl;
            std::cout << " Please review your boundary conditions for Geomechanic Module. " << std::endl;
            DebugStop();
        }
        int n_time_data = std::atoi(char_container);
        std::vector<std::pair<REAL, std::vector<REAL>>> bc_points_geo(n_time_data);
        
        /// Association bc type with numerical values
        bc_id_to_values_chunk_geo.first = bc_id;
        bc_id_to_values_chunk_geo.second.Clear();
        
        sub_container_geo = container->FirstChild("Data")->ToElement();
        int itime_data = 0;
        for( ; sub_container_geo; sub_container_geo=sub_container_geo->NextSiblingElement())
        {
                char_container = sub_container_geo->Attribute("t");
                if (!char_container)
                {
                    std::cout << " the boundary " << condition << "  needs the time value t ." << std::endl;
                    std::cout << " Please review your boundary conditions for Geomechanic Module. " << std::endl;
                    DebugStop();
                }
            REAL time_value = std::atoi(char_container);
                
            int n_data = chunk_geo->second.second.size();
            bc_points_geo[itime_data].first = time_value;
            bc_points_geo[itime_data].second.resize(n_data);
            for (int i = 0; i < n_data; i++)
            {
                char_container = sub_container_geo->Attribute(chunk_geo->second.second[i].c_str());
                if (!char_container)
                {
                    std::cout << " the boundary " << condition << "  needs the value " << chunk_geo->second.second[i] << std::endl;
                    std::cout << " Please review your boundary conditions for Geomechanic Module. " << std::endl;
                    DebugStop();
                }
                REAL bc_value = std::atof(char_container);
                bc_points_geo[itime_data].second[i] = bc_value;
            }
            itime_data++;
        }
        
        if (itime_data != n_time_data) {
            std::cout << " the boundary " << condition << " is not properly uploaded. " << std::endl;
            std::cout << " It was provided n_time_data = " << n_time_data << std::endl;
            std::cout << " It was uploaded this number of time data = " << itime_data << std::endl;
            std::cout << " Please review your boundary conditions for Geomechanic Module. " << std::endl;
            DebugStop();
        }
        TPMRSInterpolator bc_interpolator_geo(bc_points_geo);
        bc_id_to_values_chunk_geo.second = bc_interpolator_geo;
        m_bc_id_to_values_geo.insert(bc_id_to_values_chunk_geo);
        
        /// Association bc identifier with bc type
        bc_id_to_type_chunk_geo.first = bc_id;
        bc_id_to_type_chunk_geo.second = condition;
        m_bc_id_to_type_geo.insert(bc_id_to_type_chunk_geo);
        
    }
    /// End:: Regions and materials parameters of Geomechanic Simulator
    
    /// Begin:: Regions and materials parameters of Reservoir Simulator
    this->LoadBoundaryConditionsReservoirs();
    
    std::pair<int, std::string > bc_id_to_type_chunk_res;
    std::pair<int , TPMRSInterpolator > bc_id_to_values_chunk_res;
    std::map< std::string,std::pair<int,std::vector<std::string> > >::iterator chunk_res;
    
    /// Initial boundary conditions
    TiXmlElement * sub_container_res;
    container = doc_handler.FirstChild("CaseData").FirstChild("BCInitialReservoir").FirstChild("InitialReservoir").ToElement();
    for( ; container; container=container->NextSiblingElement())
    {
        
        char_container = container->Attribute("bc_id");
        int bc_id = std::atoi(char_container);
        
        char_container = container->Attribute("type");
        std::string condition(char_container);
        chunk_res = m_condition_type_to_index_value_names_res.find(condition);
        
        bool bc_condition_not_available_Q = chunk_geo == m_condition_type_to_index_value_names_res.end();
        if (bc_condition_not_available_Q)
        {
            std::cout << " Geometry dimension =  " << dimension << std::endl;
            std::cout << " The boundary " << condition << " are not available " << std::endl;
            std::cout << " Please review your boundary conditions for Reservoir Module. " << std::endl;
            DebugStop();
        }
        
        char_container = container->Attribute("n_data");
        if (!char_container)
        {
            std::cout << " The boundary " << condition << "  needs the number of data to be interporlated. " << std::endl;
            std::cout << " Please review your boundary conditions for Reservoir Module. " << std::endl;
            DebugStop();
        }
        int n_time_data = std::atoi(char_container);
        std::vector<std::pair<REAL, std::vector<REAL>>> bc_points_res(n_time_data);
        
        /// Association bc type with numerical values
        bc_id_to_values_chunk_res.first = bc_id;
        bc_id_to_values_chunk_res.second.Clear();
        
        sub_container_res = container->FirstChild("Data")->ToElement();
        int itime_data = 0;
        for( ; sub_container_res; sub_container_res=sub_container_res->NextSiblingElement())
        {
            char_container = sub_container_res->Attribute("t");
            if (!char_container)
            {
                std::cout << " the boundary " << condition << "  needs the time value t ." << std::endl;
                std::cout << " Please review your boundary conditions for Reservoir Module. " << std::endl;
                DebugStop();
            }
            REAL time_value = std::atoi(char_container);
            
            int n_data = chunk_res->second.second.size();
            bc_points_res[itime_data].first = time_value;
            bc_points_res[itime_data].second.resize(n_data);
            for (int i = 0; i < n_data; i++)
            {
                char_container = sub_container_res->Attribute(chunk_res->second.second[i].c_str());
                if (!char_container)
                {
                    std::cout << " the boundary " << condition << "  needs the value " << chunk_res->second.second[i] << std::endl;
                    std::cout << " Please review your boundary conditions for Reservoir Module. " << std::endl;
                    DebugStop();
                }
                REAL bc_value = std::atof(char_container);
                bc_points_res[itime_data].second[i] = bc_value;
            }
            itime_data++;
        }
        
        if (itime_data != n_time_data) {
            std::cout << " the boundary " << condition << " is not properly uploaded. " << std::endl;
            std::cout << " It was provided n_time_data = " << n_time_data << std::endl;
            std::cout << " It was uploaded this number of time data = " << itime_data << std::endl;
            std::cout << " Please review your boundary conditions for Reservoir Module. " << std::endl;
            DebugStop();
        }
        TPMRSInterpolator bc_interpolator_res(bc_points_res);
        bc_id_to_values_chunk_res.second = bc_interpolator_res;
        m_bc_id_to_values_res_initial.insert(bc_id_to_values_chunk_res);
        
        /// Association bc identifier with bc type
        bc_id_to_type_chunk_res.first = bc_id;
        bc_id_to_type_chunk_res.second = condition;
        m_bc_id_to_type_res_initial.insert(bc_id_to_type_chunk_res);
        
    }
    
    /// Recurrent boundary conditions
    container = doc_handler.FirstChild("CaseData").FirstChild("BCReservoirs").FirstChild("Reservoir").ToElement();
    for( ; container; container=container->NextSiblingElement())
    {
        
        char_container = container->Attribute("bc_id");
        int bc_id = std::atoi(char_container);
        
        char_container = container->Attribute("type");
        std::string condition(char_container);
        chunk_res = m_condition_type_to_index_value_names_res.find(condition);
        bool bc_condition_not_available_Q = chunk_res == m_condition_type_to_index_value_names_res.end();
        if (bc_condition_not_available_Q)
        {
            std::cout << " Geometry dimension =  " << dimension << std::endl;
            std::cout << " The boundary " << condition << " are not available " << std::endl;
            std::cout << " Please review your boundary conditions for Reservoir Module. " << std::endl;
            DebugStop();
        }
        
        char_container = container->Attribute("n_data");
        if (!char_container)
        {
            std::cout << " The boundary " << condition << "  needs the number of data to be interporlated. " << std::endl;
            std::cout << " Please review your boundary conditions for Reservoir Module. " << std::endl;
            DebugStop();
        }
        int n_time_data = std::atoi(char_container);
        std::vector<std::pair<REAL, std::vector<REAL>>> bc_points_res(n_time_data);
        
        /// Association bc type with numerical values
        bc_id_to_values_chunk_res.first = bc_id;
        bc_id_to_values_chunk_res.second.Clear();
        
        sub_container_res = container->FirstChild("Data")->ToElement();
        int itime_data = 0;
        for( ; sub_container_res; sub_container_res=sub_container_res->NextSiblingElement())
        {
            char_container = sub_container_res->Attribute("t");
            if (!char_container)
            {
                std::cout << " the boundary " << condition << "  needs the time value t ." << std::endl;
                std::cout << " Please review your boundary conditions for Reservoir Module. " << std::endl;
                DebugStop();
            }
            REAL time_value = std::atoi(char_container);
            
            int n_data = chunk_res->second.second.size();
            bc_points_res[itime_data].first = time_value;
            bc_points_res[itime_data].second.resize(n_data);
            for (int i = 0; i < n_data; i++)
            {
                char_container = sub_container_res->Attribute(chunk_res->second.second[i].c_str());
                if (!char_container)
                {
                    std::cout << " the boundary " << condition << "  needs the value " << chunk_res->second.second[i] << std::endl;
                    std::cout << " Please review your boundary conditions for Reservoir Module. " << std::endl;
                    DebugStop();
                }
                REAL bc_value = std::atof(char_container);
                bc_points_res[itime_data].second[i] = bc_value;
            }
            itime_data++;
        }
        
        if (itime_data != n_time_data) {
            std::cout << " the boundary " << condition << " is not properly uploaded. " << std::endl;
            std::cout << " It was provided n_time_data = " << n_time_data << std::endl;
            std::cout << " It was uploaded this number of time data = " << itime_data << std::endl;
            std::cout << " Please review your boundary conditions for Reservoir Module. " << std::endl;
            DebugStop();
        }
        TPMRSInterpolator bc_interpolator_res(bc_points_res);
        bc_id_to_values_chunk_res.second = bc_interpolator_res;
        m_bc_id_to_values_res.insert(bc_id_to_values_chunk_res);
        
        /// Association bc identifier with bc type
        bc_id_to_type_chunk_res.first = bc_id;
        bc_id_to_type_chunk_res.second = condition;
        m_bc_id_to_type_res.insert(bc_id_to_type_chunk_res);
        
    }
    /// End:: Regions and materials parameters of Reservoir Simulator
    
    /// Begin:: Apply uniform refinement
    this->UniformRefinement();
    /// End:: Apply uniform refinement
}

/// Brief Setup reporting times and time step size
void TPMRSSimulationData::SetTimeControls(int n_times, REAL dt, bool crank_nicolson_Q)
{
    
    m_n_steps    = n_times;
    m_dt         = dt;
    m_is_crank_nicolson_Q = crank_nicolson_Q;
    m_reporting_times.Resize(n_times, 0.0);
    for (int it = 0; it < n_times; it++)
    {
        m_reporting_times[it] = it*dt;
    }
    
}

/// Brief Setup reporting times and time step size
void TPMRSSimulationData::SetNumericControls(int n_iterations, REAL epsilon_res, REAL epsilon_cor)
{
    
    m_n_iteraions    =   n_iterations;
    m_epsilon_res    =   epsilon_res;
    m_epsilon_cor    =   epsilon_cor;
    
}

/// Brief Setup fixed stress split schemes
void TPMRSSimulationData::SetFixedStressSplitSchemes(int n_fss_iterations, int n_enf_fss_iterations)
{
    m_n_fss_iterations        =   n_fss_iterations;
    m_n_enf_fss_iterations    =   n_enf_fss_iterations;
}

/// Brief Print the all members
void TPMRSSimulationData::Print()
{
    DebugStop();
    std::cout << " TPMRSSimulationData class members : " << std::endl;
    std::cout << std::endl;
    std::cout << " m_dt = " << m_dt << std::endl;
    std::cout << " m_n_steps = " << m_n_steps << std::endl;
    std::cout << " m_reporting_times = " << m_reporting_times << std::endl;
    std::cout << " m_time = " << m_time << std::endl;
    std::cout << " m_n_iteraions = " << m_n_iteraions << std::endl;
    std::cout << " m_epsilon_res = " << m_epsilon_res << std::endl;
    std::cout << " m_epsilon_cor = " << m_epsilon_cor << std::endl;
    std::cout << " m_n_fss_iterations = " << m_n_fss_iterations << std::endl;
    std::cout << " m_n_enf_fss_iterations = " << m_n_enf_fss_iterations << std::endl;
    std::cout << " m_max_plastic_strain = " << m_max_plastic_strain << std::endl;
    std::cout << " m_n_threads = " << m_n_threads << std::endl;
    std::cout << " m_scale_factor = " << m_scale_factor << std::endl;
    std::cout << " m_is_dual_formulation_Q = " << m_is_dual_formulation_Q << std::endl;
    std::cout << " m_transfer_current_to_last_solution_Q = " << m_transfer_current_to_last_solution_Q << std::endl;
    std::cout << " m_h_level = " << m_h_level << std::endl;
    std::cout << " m_elasticity_order = " << m_elasticity_order << std::endl;
    std::cout << " m_diffusion_order = " << m_diffusion_order << std::endl;
    std::cout << " m_dimesion = " << m_dimesion << std::endl;
    std::cout << " m_geometry_file = " << m_geometry_file << std::endl;
    std::cout << " m_geometry = " << m_geometry << std::endl;
    std::cout << " m_vtk_file = " << m_vtk_file << std::endl;
    std::cout << " m_vtk_resolution = " << m_vtk_resolution << std::endl;
    std::cout << " m_n_outputs_geo = " << m_n_outputs_geo << std::endl;
    std::cout << " m_n_outputs_res = " << m_n_outputs_res << std::endl;
    std::cout << " m_g = " << m_g << std::endl;
    std::cout << " m_n_regions = " << m_n_regions << std::endl;
    std::cout << " m_mat_ids = " << std::endl;
    DebugStop();
//    int n_data = m_mat_ids.size();
//    for (int i = 0; i < n_data; i++)
//    {
//        std::cout << " region material id = " << m_mat_ids[i].first << std::endl;
//        int n_bc = m_mat_ids[i].second.size();
//        std::cout << " bc material ids = " << m_mat_ids[i].first << std::endl;
//        for (int j = 0; j <n_bc; j++)
//        {
//            std::cout << " " << m_mat_ids[i].second [j];
//        }
//        std::cout << std::endl;
//    }
    
    std::cout << " m_mat_props = " << std::endl;
    DebugStop();
//    n_data = m_mat_props.size();
//    for (int i = 0; i < n_data; i++) {
//        std::cout << " region number = " << i << std::endl;
//        int n_bc = m_mat_props[i].size();
//        for (int j = 0; j <n_bc; j++) {
//            std::cout << " " << m_mat_props[i][j];
//        }
//        std::cout << std::endl;
//    }
    std::cout << " m_is_initial_state_Q = " <<m_is_initial_state_Q << std::endl;
    std::cout << " m_is_current_state_Q = " << m_is_current_state_Q << std::endl;
    std::cout << " ---------------------- " << std::endl;
    std::cout << std::endl;
    
}


/// Brief read the geometry
void TPMRSSimulationData::ReadGeometry()
{
    TPZGmshReader Geometry;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.0");
    m_geometry = Geometry.GeometricGmshMesh(m_geometry_file);
    Geometry.PrintPartitionSummary(std::cout);
#ifdef PZDEBUG
    if (!m_geometry)
    {
        std::cout << "The geometrical mesh was not generated." << std::endl;
        DebugStop();
    }
#endif
    
}

/// Brief print the geometry
void TPMRSSimulationData::PrintGeometry()
{
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name  << "geometry" << ".txt";
    vtk_name   << "geometry"  << ".vtk";
    std::ofstream textfile(text_name.str().c_str());
    m_geometry->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(m_geometry, vtkfile, true);

#ifdef PZDEBUG
    TPZCheckGeom checker(m_geometry);
    checker.CheckUniqueId();
    if(checker.PerformCheck())
    {
        DebugStop();
    }
#endif
    
}


/// Brief applying the boundary conditions for reservoir simulator
void TPMRSSimulationData::LoadBoundaryConditionsReservoirs()
{
 
    std::pair<std::string,std::pair<int,std::vector<std::string> > > chunkReser;

        /// Dirichlet for for diffusion
        chunkReser.first = "Dp"; // name
        chunkReser.second.first = 0; // index
        chunkReser.second.second.push_back("p");
        m_condition_type_to_index_value_names_res.insert(chunkReser);
        chunkReser.second.second.resize(0);
    
        /// Neumann for diffusion
        chunkReser.first = "Nq"; // name
        chunkReser.second.first = 1; // index
        chunkReser.second.second.push_back("qn");
        m_condition_type_to_index_value_names_res.insert(chunkReser);
        chunkReser.second.second.resize(0);
    
    return;
}


/// Brief applying the boundary conditions for geomechanics simulator
void TPMRSSimulationData::LoadBoundaryConditionsGeomechanics()
{
    std::pair<std::string,std::pair<int,std::vector<std::string> > > chunkGeo;
    
    if (m_dimesion == 2)
    {
        /// 2D conditions
        
        /// Dirichlet for elasticity
        chunkGeo.first = "Du"; // name
        chunkGeo.second.first = 2; // index
        chunkGeo.second.second.push_back("ux");
        chunkGeo.second.second.push_back("uy");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Dirichlet for elasticity in x_direction
        chunkGeo.first = "Dux"; // name
        chunkGeo.second.first = 3; // index
        chunkGeo.second.second.push_back("ux");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Dirichlet for elasticity in y_direction
        chunkGeo.first = "Duy"; // name
        chunkGeo.second.first = 4; // index
        chunkGeo.second.second.push_back("uy");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Neumann for elasticity
        chunkGeo.first = "Nt"; // name
        chunkGeo.second.first = 5; // index
        chunkGeo.second.second.push_back("tx");
        chunkGeo.second.second.push_back("ty");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Neumann for elasticity
        chunkGeo.first = "Ntn"; // name
        chunkGeo.second.first = 6; // index
        chunkGeo.second.second.push_back("tn");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Neumann for elasticity full tensor
        chunkGeo.first = "NS"; // name
        chunkGeo.second.first = 7; // index
        chunkGeo.second.second.push_back("sxx");
        chunkGeo.second.second.push_back("sxy");
        chunkGeo.second.second.push_back("sxz");
        chunkGeo.second.second.push_back("syy");
        chunkGeo.second.second.push_back("syz");
        chunkGeo.second.second.push_back("szz");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Dirichlet for elasticity (normal displacement)
        chunkGeo.first = "Dun"; // name
        chunkGeo.second.first = 8; // index
        chunkGeo.second.second.push_back("un");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
    }else
    {
        /// 3D conditions
        
        /// Dirichlet for elasticity
        chunkGeo.first = "Du"; // name
        chunkGeo.second.first = 2; // index
        chunkGeo.second.second.push_back("ux");
        chunkGeo.second.second.push_back("uy");
        chunkGeo.second.second.push_back("uz");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Dirichlet for elasticity in x_direction
        chunkGeo.first = "Dux"; // name
        chunkGeo.second.first = 3; // index
        chunkGeo.second.second.push_back("ux");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Dirichlet for elasticity in y_direction
        chunkGeo.first = "Duy"; // name
        chunkGeo.second.first = 4; // index
        chunkGeo.second.second.push_back("uy");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Neumann for elasticity
        chunkGeo.first = "Nt"; // name
        chunkGeo.second.first = 5; // index
        chunkGeo.second.second.push_back("tx");
        chunkGeo.second.second.push_back("ty");
        chunkGeo.second.second.push_back("tz");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Neumann for elasticity
        chunkGeo.first = "Ntn"; // name
        chunkGeo.second.first = 6; // index
        chunkGeo.second.second.push_back("tn");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Dirichlet for elasticity in z_direction
        chunkGeo.first = "Duz"; // name
        chunkGeo.second.first = 7; // index
        chunkGeo.second.second.push_back("uz");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
    
        /// Dirichlet for elasticity in x & y direction
        chunkGeo.first = "Duxy"; // name
        chunkGeo.second.first = 8; // index
        chunkGeo.second.second.push_back("ux");
        chunkGeo.second.second.push_back("uy");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Dirichlet for elasticity in x & z direction
        chunkGeo.first = "Duxz"; // name
        chunkGeo.second.first = 9; // index
        chunkGeo.second.second.push_back("ux");
        chunkGeo.second.second.push_back("uz");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Dirichlet for elasticity in y & z direction
        chunkGeo.first = "Duyz"; // name
        chunkGeo.second.first = 10; // index
        chunkGeo.second.second.push_back("uy");
        chunkGeo.second.second.push_back("uz");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Neumann for elasticity full tensor
        chunkGeo.first = "NS"; // name
        chunkGeo.second.first = 11; // index
        chunkGeo.second.second.push_back("sxx");
        chunkGeo.second.second.push_back("sxy");
        chunkGeo.second.second.push_back("sxz");
        chunkGeo.second.second.push_back("syy");
        chunkGeo.second.second.push_back("syz");
        chunkGeo.second.second.push_back("szz");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
        /// Dirichlet for elasticity (normal displacement)
        chunkGeo.first = "Dun"; // name
        chunkGeo.second.first = 12; // index
        chunkGeo.second.second.push_back("un");
        m_condition_type_to_index_value_names_geo.insert(chunkGeo);
        chunkGeo.second.second.resize(0);
        
    }
    
    return;
}

void TPMRSSimulationData::ReadRegionsAndMaterials(){
    DebugStop();
}

void TPMRSSimulationData::ReadBCForGeomechanicSimulator(){
    DebugStop();
}

void TPMRSSimulationData::ReadBCForReservoirSimulator(){
    DebugStop();
}

void TPMRSSimulationData::UniformRefinement() {
    
    TPZManVector<TPZGeoEl*> sons;
    for(int i=0; i < m_h_level; i++)
    {
        int64_t nels = m_geometry->NElements();
        for(int64_t elem = 0; elem < nels; elem++)
        {
            TPZGeoEl * gel = m_geometry->ElementVec()[elem];
            if(!gel || gel->HasSubElement())
                continue;
            gel->Divide(sons);
        }
    }
    m_geometry->ResetConnectivities();
    m_geometry->BuildConnectivity();
}
