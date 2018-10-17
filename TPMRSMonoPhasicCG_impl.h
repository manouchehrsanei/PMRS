//
//  TPMRSMonoPhasicCG_impl.h
//  PMRS
//
//  Created by Omar Durán on 9/22/18.
//


#include "TPMRSMonoPhasicCG.h"

template <class TMEM>
TPMRSMonoPhasicCG<TMEM>::TPMRSMonoPhasicCG() : m_phi_model(), m_kappa_model(){
    m_simulation_data   = NULL;
    m_dimension         = 0;
    m_c                 = 0;
    m_eta               = 0;
    m_rho_0             = 0;
    m_scale_factor      = 1;
}

template <class TMEM>
TPMRSMonoPhasicCG<TMEM>::TPMRSMonoPhasicCG(int mat_id, int dimension) :  TPZMatWithMem<TMEM>(mat_id), m_phi_model(), m_kappa_model(){
    m_simulation_data   = NULL;
    m_dimension = dimension;
    m_c                 = 0;
    m_eta               = 0;
    m_rho_0             = 0;
    m_scale_factor      = 1;
}

template <class TMEM>
TPMRSMonoPhasicCG<TMEM>::TPMRSMonoPhasicCG(const TPMRSMonoPhasicCG & other){
    m_simulation_data   = other.m_simulation_data;
    m_dimension         = other.m_dimension;
    m_c                 = other.m_c;
    m_eta               = other.m_eta;
    m_rho_0             = other.m_rho_0;
    m_scale_factor      = other.m_scale_factor;
    m_phi_model         = other.m_phi_model;
    m_kappa_model       = other.m_kappa_model;
}

template <class TMEM>
TPMRSMonoPhasicCG<TMEM> & TPMRSMonoPhasicCG<TMEM>::operator=(const TPMRSMonoPhasicCG & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_simulation_data   = other.m_simulation_data;
    m_dimension         = other.m_dimension;
    m_c                 = other.m_c;
    m_eta               = other.m_eta;
    m_rho_0             = other.m_rho_0;
    m_scale_factor      = other.m_scale_factor;
    m_phi_model         = other.m_phi_model;
    m_kappa_model       = other.m_kappa_model;
    return *this;
}

template <class TMEM>
TPMRSMonoPhasicCG<TMEM>::~TPMRSMonoPhasicCG(){
    
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::Print(std::ostream &out){
    
    out << " Material name : " << this->Name() << "\n";
    out << " Pointer to TPZSimulationData : " << m_simulation_data << "\n";
    out << " Material dimension : " << m_dimension << "\n";
    out << " Fluid compressibility : " << m_c << "\n";
    out << " Fluid viscosity : " << m_eta << "\n";
    out << " Fluid density : " << m_rho_0 << "\n";
    out << " Scale factor  : " << m_scale_factor << "\n";
    m_phi_model.Print(out);
    m_kappa_model.Print(out);
    out << "\t Base class print:\n";
    TPZMaterial::Print(out);
    
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::FillDataRequirements(TPZMaterialData &data){
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data){
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::UndrainedContribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    // Getting weight functions
    TPZFMatrix<REAL>  & phi_p     =  data.phi;
    int n_phi_p = phi_p.Rows();
    STATE p_n                  = data.sol[0][0];
    
    // Get the pressure at the integrations points
    long gp_index = data.intGlobPtIndex;
    TMEM & memory = this->MemItem(gp_index);
    STATE p_0      = memory.p_0();
    
    for (int ip = 0; ip < n_phi_p; ip++)
    {
        ef(ip) +=  weight * ( p_n - p_0 ) * phi_p(ip,0) ;
        
        for (int jp = 0; jp < n_phi_p; jp++)
        {
            ek(ip,jp) +=  weight * phi_p(jp,0) * phi_p(ip,0);
        }
        
    }
    
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    // Undarined contribute
    if(m_simulation_data->IsInitialStateQ()){
        this->UndrainedContribute(data, weight, ek, ef);
        return;
    }
    
    // Getting weight functions
    TPZFMatrix<REAL>  & phi_p     =  data.phi;
    int n_phi_p = phi_p.Rows();
    
    TPZFNMatrix<40,REAL> grad_phi_p(m_dimension,n_phi_p);
    TPZFNMatrix<9,REAL> grad_p(m_dimension,1);
    TPZAxesTools<REAL>::Axes2XYZ(data.dphix, grad_phi_p, data.axes);
    TPZAxesTools<REAL>::Axes2XYZ(data.dsol[0], grad_p, data.axes);
    
    
    // Get the pressure at the integrations points
    long gp_index = data.intGlobPtIndex;
    TMEM & memory = this->MemItem(gp_index);
    
    // Time
    STATE dt = m_simulation_data->dt();
    STATE p_n                  = data.sol[0][0];
    
    STATE p_0      = memory.p_0();
    STATE p        = memory.p();
    
    TPZFNMatrix<9,REAL> K(3,3),Kinv(3,3);
    
    K.Zero();
    K(0,0) = memory.kappa();
    K(1,1) = memory.kappa();
    K(2,2) = memory.kappa();
    
    Kinv.Zero();
    Kinv(0,0) = 1.0/memory.kappa();
    Kinv(1,1) = 1.0/memory.kappa();
    Kinv(2,2) = 1.0/memory.kappa();
    
    STATE rho = m_rho_0 * (1 + m_c*(p-p_0));
    STATE rho_n = m_rho_0 * (1 + m_c*(p_n-p_0));
    STATE drho_ndp_n = m_c;
    STATE lambda = rho_n/m_eta;
    
    // Defining local variables
    TPZFNMatrix<3,STATE> Kl_grad_p_(3,1);
    for (int i = 0; i < Dimension(); i++) {
        STATE dot = 0.0;
        for (int j =0; j < Dimension(); j++) {
            dot    += K(i,j)*grad_p(j,0);
        }
        Kl_grad_p_(i,0)     = lambda * dot;
    }
    
    // Integration point contribution
    TPZFNMatrix<3,STATE> phi_q_i(3,1), phi_q_j(3,1);
    
    if(!m_simulation_data->IsCurrentStateQ()){
        return;
    }
    
    STATE phi_n,dphi_ndp,phi;
    this->porosity(gp_index,phi_n,dphi_ndp,phi);
    
    TPZFNMatrix<3,STATE> Kl_grad_phi_j_(3,1);
    for (int ip = 0; ip < n_phi_p; ip++)
    {
        
        STATE Kl_grad_p_dot_grad_phi = 0.0;
        for (int i = 0; i < Dimension(); i++) {
            Kl_grad_p_dot_grad_phi    += Kl_grad_p_(i,0)*grad_phi_p(i,ip);
            
        }
        
        ef(ip) +=  weight * ( Kl_grad_p_dot_grad_phi + m_scale_factor * (1.0/dt) * ( phi_n*rho_n - phi*rho ) * phi_p(ip,0) );
        
        for (int jp = 0; jp < n_phi_p; jp++)
        {
        
            for (int i = 0; i < Dimension(); i++) {
                STATE dot = 0.0;
                for (int j =0; j < Dimension(); j++) {
                    dot    += K(i,j)*grad_phi_p(j,jp);
                }
                Kl_grad_phi_j_(i,0)     = lambda * dot;
            }
            
            STATE Kl_grad_phi_j_dot_grad_phi_j = 0.0;
            for (int i = 0; i < Dimension(); i++) {
                Kl_grad_phi_j_dot_grad_phi_j    += Kl_grad_phi_j_(i,0)*grad_phi_p(i,ip);
                
            }
            
            ek(ip,jp) +=  weight * ( Kl_grad_phi_j_dot_grad_phi_j + m_scale_factor * (1.0/dt) * ( phi_n * drho_ndp_n + dphi_ndp * rho_n ) * phi_p(jp,0)  * phi_p(ip,0) );
        }
        
    }
    
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
    

    if (m_simulation_data->Get_must_accept_solution_Q()) {
        long gp_index = data.intGlobPtIndex;
        
        if (m_simulation_data->GetTransferCurrentToLastQ()) {
            this->MemItem(gp_index).Setp(this->MemItem(gp_index).p_n()) ;
            return;
        }
        
        STATE p                  = data.sol[0][0];
        
        if (m_simulation_data->IsInitialStateQ()) {
            this->MemItem(gp_index).Setp_0(p);
        }
        
        if (m_simulation_data->IsCurrentStateQ()) {
            this->MemItem(gp_index).Setp_n(p);
        }else{
            this->MemItem(gp_index).Setp(p);
        }
        
    }

    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(data, weight, ek_fake, ef);
    
    
    return;
    
    
    
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    // Undarined contribute
    if(m_simulation_data->IsInitialStateQ()){
        return;
    }
    
    TPZFNMatrix<100,STATE> phi_p       = data.phi;
    int n_phi_p       = phi_p.Rows();
    
    REAL p  = data.sol[0][0];
    REAL BigNumber = TPZDiscontinuousGalerkin::gBigNumber;
    REAL Value = bc.Val2()(0,0);
    
    switch (bc.Type())
    {
  
        case 0 :    // Dirichlet BC  PD
        {
            REAL p_D = m_scale_factor * Value;
            for (int ip = 0; ip < n_phi_p; ip++)
            {
                ef(ip) += m_scale_factor * weight * BigNumber * (p - p_D) * phi_p(ip,0);
                
                for (int jp = 0; jp < n_phi_p; jp++)
                {
                    
                    ek(ip,jp) += m_scale_factor * weight * BigNumber * phi_p(jp,0) * phi_p(ip,0);
                }
            }
            
        }
            
            break;
            
        case 1 :    // Neumann BC  Normal flux qn_N
        {
            STATE qn_N = Value;
            for (int ip = 0; ip < n_phi_p; ip++)
            {
                ef(ip) += -1.0 * m_scale_factor * weight * qn_N * phi_p(ip,0);
            }
        }
            
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
    
    return;
    
}

template <class TMEM>
int TPMRSMonoPhasicCG<TMEM>::VariableIndex(const std::string &name) {
    if (!strcmp("p", name.c_str())) return 0;
    if (!strcmp("q", name.c_str())) return 1;
    if (!strcmp("div_q", name.c_str())) return 2;
    if (!strcmp("kappa", name.c_str())) return 3;
    if (!strcmp("phi", name.c_str())) return 4;
    if (!strcmp("order", name.c_str())) return 5;
    if (!strcmp("id", name.c_str())) return 6;
    return TPZMatWithMem<TMEM>::VariableIndex(name);
}

template <class TMEM>
int TPMRSMonoPhasicCG<TMEM>::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return m_dimension; // Vector
        case 2:
            return 1; // Scalar
        case 3:
            return 1; // Scalar
        case 4:
            return 1; // Scalar
        case 5:
            return 1; // Scalar
        case 6:
            return 1; // Scalar
    }
    return TPZMatWithMem<TMEM>::NSolutionVariables(var);
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout){
    
    long gp_index = data.intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    Solout.Resize( this->NSolutionVariables(var));
    
    switch (var) {
        case 0:
        {
            Solout[0] = memory.p_n();
        }
            break;
        case 3:
        {
            Solout[0] = memory.kappa();
        }
            break;
        case 4:
        {
            Solout[0] = memory.phi();
        }
            break;
            
        default:
        {
            std::cout << "TPMRSMonoPhasicCG<TMEM>:: Variable not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
    
    
    
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::porosity(long gp_index, REAL &phi_n, REAL &dphi_ndp, REAL &phi){
    
    TMEM & memory = this->MemItem(gp_index);
    
    REAL alpha = memory.GetAlpha();
    REAL Se = memory.GetSe();
    REAL phi_0 = memory.phi_0();
    
    REAL p_0 = memory.p_0();
    REAL p = memory.p();
    REAL p_n = memory.p_n();
    REAL sigma_t_v_0 = (memory.GetSigma_0().I1()/3) ;//- alpha * p_0;
    REAL sigma_t_v = (memory.GetSigma().I1()/3) ;//- alpha * p;
    REAL sigma_t_v_n = (memory.GetSigma_n().I1()/3) ;//- alpha * p;
    
    m_phi_model.Porosity(phi, dphi_ndp, phi_0, p, p_0, sigma_t_v, sigma_t_v_0, alpha, Se);
    m_phi_model.Porosity(phi_n, dphi_ndp, phi_0, p_n, p_0, sigma_t_v_n, sigma_t_v_0, alpha, Se);
    
    this->MemItem(gp_index).Setphi(phi); // Current phi, please rename it ot phi_n
}
