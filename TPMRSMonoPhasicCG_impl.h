//
//  TPMRSMonoPhasicCG_impl.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/22/18.
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
    m_theta_scheme      = 1;
}

template <class TMEM>
TPMRSMonoPhasicCG<TMEM>::TPMRSMonoPhasicCG(int mat_id, int dimension) :  TPZMatWithMem<TMEM>(mat_id), m_phi_model(), m_kappa_model(){
    m_simulation_data   = NULL;
    m_dimension         = dimension;
    m_c                 = 0;
    m_eta               = 0;
    m_rho_0             = 0;
    m_scale_factor      = 1;
    m_theta_scheme      = 1;
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
    m_theta_scheme      = other.m_theta_scheme;
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
    m_theta_scheme      = other.m_theta_scheme;
    return *this;
}

template <class TMEM>
TPMRSMonoPhasicCG<TMEM>::~TPMRSMonoPhasicCG(){
    
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::Print(std::ostream &out){
    
    out << " Material name : " << this->Name() << "\n";
    out << " Pointer to TPMRSSimulationData : " << m_simulation_data << "\n";
    out << " Material dimension : " << m_dimension << "\n";
    out << " Fluid compressibility : " << m_c << "\n";
    out << " Fluid viscosity : " << m_eta << "\n";
    out << " Fluid density : " << m_rho_0 << "\n";
    out << " Scale factor  : " << m_scale_factor << "\n";
    out << " Theta scheme  : " << m_theta_scheme << "\n";
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
    
    /// Getting weight functions
    TPZFMatrix<REAL>  & phi_p  = data.phi;
    int n_phi_p                = phi_p.Rows();
    STATE p_n                  = data.sol[0][0];
    
    /// Get the pressure at the integrations points
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
    
    /// Undarined contribute
    if(m_simulation_data->IsInitialStateQ()){
        this->UndrainedContribute(data, weight, ek, ef);
        return;
    }
    
    /// Getting weight functions
    TPZFMatrix<REAL>  & phi_p     =  data.phi;
    int n_phi_p = phi_p.Rows();
    
    TPZFNMatrix<40,REAL> grad_phi_p(m_dimension,n_phi_p);
    TPZFNMatrix<9,REAL> grad_p(m_dimension,1);
    TPZAxesTools<REAL>::Axes2XYZ(data.dphix, grad_phi_p, data.axes);
    TPZAxesTools<REAL>::Axes2XYZ(data.dsol[0], grad_p, data.axes);
    
    
    /// Get the pressure at the integrations points
    long gp_index  = data.intGlobPtIndex;
    TMEM & memory  = this->MemItem(gp_index);
    
    /// Time
    STATE dt       = m_simulation_data->dt();
    STATE p_n      = data.sol[0][0];
    
    STATE p_0      = memory.p_0();
    STATE p        = memory.p();
    
    TPZManVector<STATE,3> & last_Kl_grad_p        = memory.f_vec();
    
    STATE phi_n,dphi_ndp,phi;
    REAL phi_0 = memory.phi_0();
    this->porosity(gp_index,phi_n,dphi_ndp,phi);
    
    REAL kappa_n;
    REAL dkappa_ndphi,dkappa_ndp;
    this->permeability(gp_index, kappa_n, dkappa_ndphi, phi_n, phi_0);
    kappa_n *= (1.0/m_scale_factor);
    dkappa_ndp = (1.0/m_scale_factor)*dkappa_ndphi * dphi_ndp;
    
    TPZFNMatrix<9,REAL> K(3,3),dKdp(3,3);
    
    K.Zero();
    K(0,0) = (1.0/m_scale_factor)*memory.kappa_n();
    K(1,1) = (1.0/m_scale_factor)*memory.kappa_n();
    K(2,2) = (1.0/m_scale_factor)*memory.kappa_n();
    
    dKdp.Zero();
    dKdp(0,0) = dkappa_ndp;
    dKdp(1,1) = dkappa_ndp;
    dKdp(2,2) = dkappa_ndp;
    
    STATE rho        = m_rho_0 * (1 + (m_c/m_scale_factor)*(p-p_0));
    STATE rho_n      = m_rho_0 * (1 + (m_c/m_scale_factor)*(p_n-p_0));
    STATE drho_ndp_n = m_c/m_scale_factor;
    STATE lambda     = rho_n/m_eta;
    
    /// Defining local variables
    TPZFNMatrix<3,STATE> Kl_grad_p_(3,1),dKdpl_grad_p_(3,1);
    for (int i = 0; i < Dimension(); i++)
    {
        STATE dot     = 0.0;
        STATE dKdpdot = 0.0;
        for (int j = 0; j < Dimension(); j++)
        {
            dot        += K(i,j)*grad_p(j,0);
            dKdpdot    += dKdp(i,j)*grad_p(j,0);
        }
        
        Kl_grad_p_(i,0)      = lambda * dot;
        dKdpl_grad_p_(i,0)   = lambda * dKdpdot;
    }
    
    /// Integration point contribution
    TPZFNMatrix<3,STATE> phi_q_i(3,1), phi_q_j(3,1);
    
    if(!m_simulation_data->IsCurrentStateQ()){
        return;
    }
    
    TPZFNMatrix<3,STATE> Kl_grad_phi_j_(3,1);
    for (int ip = 0; ip < n_phi_p; ip++)
    {
        
        STATE Kl_grad_p_dot_grad_phi        = 0.0;
        STATE last_Kl_grad_p_dot_grad_phi   = 0.0;
        STATE dKdpl_grad_p_dot_grad_phi     = 0.0;
        for (int i = 0; i < Dimension(); i++)
        {
            Kl_grad_p_dot_grad_phi      += Kl_grad_p_(i,0)*grad_phi_p(i,ip);
            last_Kl_grad_p_dot_grad_phi += last_Kl_grad_p[i]*grad_phi_p(i,ip);
            dKdpl_grad_p_dot_grad_phi   += dKdpl_grad_p_(i,0)*grad_phi_p(i,ip);
            
        }

        REAL Kl_grad_p_dot_grad_phi_avg = m_theta_scheme*Kl_grad_p_dot_grad_phi + (1.0-m_theta_scheme)*last_Kl_grad_p_dot_grad_phi;
        ef(ip) +=  weight * ( Kl_grad_p_dot_grad_phi_avg + (1.0/dt) * ( phi_n*rho_n - phi*rho ) * phi_p(ip,0) );
        
        for (int jp = 0; jp < n_phi_p; jp++)
        {
        
            for (int i = 0; i < Dimension(); i++)
            {
                STATE dot = 0.0;
                for (int j =0; j < Dimension(); j++)
                {
                    dot    += K(i,j)*grad_phi_p(j,jp);
                }
                
                Kl_grad_phi_j_(i,0)     = lambda * dot;
            }
            
            STATE Kl_grad_phi_j_dot_grad_phi_j = 0.0;
            for (int i = 0; i < Dimension(); i++)
            {
                Kl_grad_phi_j_dot_grad_phi_j    += Kl_grad_phi_j_(i,0)*grad_phi_p(i,ip);
                
            }
            
            ek(ip,jp) +=  weight * ( m_theta_scheme*Kl_grad_phi_j_dot_grad_phi_j + m_theta_scheme*dKdpl_grad_p_dot_grad_phi + (1.0/dt) * ( phi_n * drho_ndp_n + dphi_ndp * rho_n ) * phi_p(jp,0)  * phi_p(ip,0) );
        }
        
    }
    
    if (m_simulation_data->GetTransferCurrentToLastQ()) {
        for (int i = 0; i < Dimension(); i++)
        {
            last_Kl_grad_p[i] = Kl_grad_p_(i,0);
        }
        this->MemItem(gp_index).Setf_vec(last_Kl_grad_p); // Update process for f function.
        this->MemItem(gp_index).Setp(this->MemItem(gp_index).p_n());
        this->MemItem(gp_index).Setq(this->MemItem(gp_index).q_n());
    }
    
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
    

    if (m_simulation_data->Get_must_accept_solution_Q()) {
        long gp_index = data.intGlobPtIndex;
        
        // Pressure variable
        STATE p = data.sol[0][0];
        
        // flux variable
        TPZManVector<REAL,3> q(3);
        {
            REAL p_0      = this->MemItem(gp_index).p_0();
            REAL p_n      = p;
            REAL rho_n    = m_rho_0 * (1 + m_c/m_scale_factor*(p_n-p_0)); //  Provide the compressibility in MPa
            REAL lambda   = rho_n/m_eta;
            REAL k        = this->MemItem(gp_index).kappa_n();
            TPZFNMatrix<9,REAL> grad_p(m_dimension,1);
            TPZAxesTools<REAL>::Axes2XYZ(data.dsol[0], grad_p, data.axes);
            for (int i = 0; i < m_dimension; i++)
            {
                q[i] = (- (1.0/m_scale_factor) * k * lambda * grad_p[i]);
            }
        }
        
        if (m_simulation_data->IsInitialStateQ()) {
            this->MemItem(gp_index).Setp_0(p);
        }
        
        if (m_simulation_data->IsCurrentStateQ()) {
            this->MemItem(gp_index).Setp_n(p);
            this->MemItem(gp_index).Setq_n(q);
        }else{
            this->MemItem(gp_index).Setp(p);
            this->MemItem(gp_index).Setq(q);
        }
        
    }

    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(data, weight, ek_fake, ef);
    
    
    return;
    
    
    
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    /// Undarined contribute
    if(m_simulation_data->IsInitialStateQ()){
        return;
    }
    
    int gp_index = data.intGlobPtIndex;
    TPZBndCondWithMem<TPMRSMonoPhasicMemory>  & bc_mem = dynamic_cast<TPZBndCondWithMem<TPMRSMonoPhasicMemory> & >(bc);
    TPMRSMonoPhasicMemory & memory = bc_mem.MemItem(gp_index);

    
    TPZFNMatrix<100,STATE> phi_p       = data.phi;
    int n_phi_p       = phi_p.Rows();
    
    REAL p  = data.sol[0][0];
    REAL BigNumber = TPZDiscontinuousGalerkin::gBigNumber;
    REAL Value = bc.Val2()(0,0);
    
    switch (bc.Type())
    {
  
        case 0 :    /// Dirichlet BC  PD
        {
            REAL p_D = Value;
            for (int ip = 0; ip < n_phi_p; ip++)
            {
                ef(ip) += weight * m_scale_factor * BigNumber * (p - p_D) * phi_p(ip,0);
                
                for (int jp = 0; jp < n_phi_p; jp++)
                {
                    
                    ek(ip,jp) += weight * m_scale_factor *  BigNumber * phi_p(jp,0) * phi_p(ip,0);
                }
            }
            
        }
            
            break;
            
        case 1 :    /// Neumann BC  Normal flux qn_N
        {
            
            STATE last_qn_N = memory.f();
            STATE current_qn = Value;
            STATE qn_N = m_theta_scheme*current_qn+(1.0-m_theta_scheme)*last_qn_N;
            for (int ip = 0; ip < n_phi_p; ip++)
            {
                ef(ip) += -1.0 * weight * qn_N * phi_p(ip,0);
            }
            if (m_simulation_data->GetTransferCurrentToLastQ()) {
                memory.Setf(current_qn);
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
    if (!strcmp("p"     , name.c_str())) return 0;
    if (!strcmp("phi"   , name.c_str())) return 1;
    if (!strcmp("kappa" , name.c_str())) return 2;
    if (!strcmp("div_q" , name.c_str())) return 3;
    if (!strcmp("qx"    , name.c_str())) return 4;
    if (!strcmp("qy"    , name.c_str())) return 5;
    if (!strcmp("qz"    , name.c_str())) return 6;
    if (!strcmp("q"    , name.c_str()))  return 7;
    return TPZMatWithMem<TMEM>::VariableIndex(name);
}

template <class TMEM>
int TPMRSMonoPhasicCG<TMEM>::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return 1; /// Scalar
        case 1:
            return 1; /// Scalar
        case 2:
            return 1; /// Scalar
        case 3:
            return 1; /// Scalar
        case 4:
            return 1; /// Scalar
        case 5:
            return 1; /// Scalar
        case 6:
            return 1; /// Scalar
        case 7:
            return m_dimension; /// Vector

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
        case 1:
        {
            Solout[0] = memory.phi_n();
        }
            break;
        case 2:
        {
            Solout[0] = memory.kappa_n();
        }
            break;
        case 3:
        {
            Solout[0] = memory.div_q_n(); //  Variable without meaning.
        }
            break;
        case 4:
        {
            Solout[0] = memory.q_n()[0];
        }
            break;
        case 5:
        {
            Solout[0] = memory.q_n()[1];

        }
            break;
        case 6:
        {
            Solout[0] = memory.q_n()[2];
        }
            break;
        case 7:
        {
            for (int i = 0; i < m_dimension; i++)
            {
                Solout[i] = memory.q_n()[i];
            }
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
    
    REAL alpha = memory.Alpha();
    REAL Kdr   = memory.Kdr();
    REAL phi_0 = memory.phi_0();
    
    REAL p_0   = memory.p_0();
    REAL p     = memory.p();
    REAL p_n   = memory.p_n();
    
    REAL sigma_t_v_0 = (memory.GetSigma_0().I1()/3) - alpha * p_0;
    REAL sigma_t_v   = (memory.GetSigma().I1()/3)  - alpha * p;
    
    REAL geo_delta_phi   = (alpha/Kdr)*(sigma_t_v-sigma_t_v_0);
    REAL geo_delta_phi_n = memory.delta_phi() + (alpha/Kdr)*(sigma_t_v-sigma_t_v_0);

    m_phi_model.Porosity(phi, dphi_ndp, phi_0, p, p_0, alpha, Kdr, geo_delta_phi);
    m_phi_model.Porosity(phi_n, dphi_ndp, phi_0, p_n, p_0, alpha, Kdr, geo_delta_phi_n);
    this->MemItem(gp_index).Setphi_n(phi_n);
}

template <class TMEM>
void TPMRSMonoPhasicCG<TMEM>::permeability(long gp_index, REAL &kappa_n, REAL &dkappa_ndphi,REAL &phi,REAL &phi_0){
    
    TMEM & memory = this->MemItem(gp_index);
    REAL kappa_0  = memory.kappa_0();
    m_kappa_model.Permeability(kappa_n, dkappa_ndphi, kappa_0, phi, phi_0);
    this->MemItem(gp_index).Setkappa_n(kappa_n);
}
