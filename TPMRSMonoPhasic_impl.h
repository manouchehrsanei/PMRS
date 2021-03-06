//
//  TPMRSMonoPhasic_impl.hpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPMRSMonoPhasic.h"

template <class TMEM>
TPMRSMonoPhasic<TMEM>::TPMRSMonoPhasic() : m_phi_model(), m_kappa_model(){
    m_simulation_data   = NULL;
    m_dimension         = 0;
    m_c                 = 0;
    m_eta               = 0;
    m_rho_0             = 0;
    m_scale_factor      = 1;
    m_theta_scheme      = 1;
}

template <class TMEM>
TPMRSMonoPhasic<TMEM>::TPMRSMonoPhasic(int mat_id, int dimension) : TPZMatWithMem<TMEM>(mat_id), m_phi_model(), m_kappa_model(){
    m_simulation_data   = NULL;
    m_dimension = dimension;
    m_c                 = 0;
    m_eta               = 0;
    m_rho_0             = 0;
    m_scale_factor      = 1;
    m_theta_scheme      = 1;
}

template <class TMEM>
TPMRSMonoPhasic<TMEM>::TPMRSMonoPhasic(const TPMRSMonoPhasic & other) : TPZMatWithMem<TMEM>(other){
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
TPMRSMonoPhasic<TMEM> & TPMRSMonoPhasic<TMEM>::operator=(const TPMRSMonoPhasic & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    TPZMatWithMem<TMEM>::operator=(other);
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
TPMRSMonoPhasic<TMEM>::~TPMRSMonoPhasic(){
    
}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::Print(std::ostream &out){
    
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
void TPMRSMonoPhasic<TMEM>::FillDataRequirements(TPZVec<TPZMaterialData> &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::UndrainedContribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    unsigned int q_b = 0;
    unsigned int p_b = 1;
    
    TPZFNMatrix<100,STATE> phi_ps       = datavec[p_b].phi;
    int nphi_q       = datavec[q_b].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    int firstq       = 0;
    int firstp       = nphi_q + firstq;
    
    /// Get the pressure at the integrations points
    long gp_index = datavec[q_b].intGlobPtIndex;
    TMEM & memory = this->MemItem(gp_index);
    
    STATE p_n      = datavec[p_b].sol[0][0];
    STATE p_0      = memory.p_0();
    
    for (int iq = 0; iq < nphi_q; iq++)
    {
        ef(iq + firstq) =  0.0;
        ek(iq + firstq,iq + firstq)  = 1.0;        
    }
    
    for (int ip = 0; ip < nphi_p; ip++)
    {
        
        ef(ip + firstp) += weight * ( p_n - p_0 ) * phi_ps(ip,0);
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(ip + firstp, jp + firstp) += weight * phi_ps(jp,0) * phi_ps(ip,0);
        }
        
    }
    
    
}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    /// Undarined contribute
    if(m_simulation_data->IsInitialStateQ()){
        this->UndrainedContribute(datavec, weight, ek, ef);
        return;
    }
    
    unsigned int q_b = 0;
    unsigned int p_b = 1;
    
    
    TPZFNMatrix<100,STATE> phi_qs       = datavec[q_b].phi;
    TPZFNMatrix<100,STATE> phi_ps       = datavec[p_b].phi;
    TPZFNMatrix<10,STATE> grad_q_axes   = datavec[q_b].dsol[0];
    TPZFNMatrix<40, REAL> div_phi = datavec[q_b].divphi;
    
    
    /// Get the pressure at the integrations points
    long gp_index = datavec[q_b].intGlobPtIndex;
    TMEM & memory = this->MemItem(gp_index);
    
    /// Time
    STATE dt = m_simulation_data->dt();
    TPZManVector<STATE,3> q   = datavec[q_b].sol[0];
    STATE p_n                 = datavec[p_b].sol[0][0];
    
    STATE p_c      = memory.p_n();
    STATE p_0      = memory.p_0();
    STATE p        = memory.p();
    
    STATE phi_n,dphi_ndp,phi;
    REAL phi_0 = memory.phi_0();
    this->porosity(gp_index,phi_n,dphi_ndp,phi);
    
    /// TPZFNMatrix<9,REAL> K(3,3),Kinv(3,3),Kinv_c(3,3),dKinvdp(3,3)
    REAL kappa_n;
    REAL dkappa_ndphi,dkappa_ndp;
    this->permeability(gp_index, kappa_n, dkappa_ndphi, phi_n, phi_0);
    dkappa_ndp = dkappa_ndphi * dphi_ndp;
    
    
//    K.Zero();
//    K(0,0) = kappa_n;
//    K(1,1) = kappa_n;
//    K(2,2) = kappa_n;
//    
//    Kinv.Zero();
//    Kinv(0,0) = 1.0/kappa_n;
//    Kinv(1,1) = 1.0/kappa_n;
//    Kinv(2,2) = 1.0/kappa_n;
//    
//    Kinv_c.Zero();
//    Kinv_c(0,0) = 1.0/memory.kappa_0();
//    Kinv_c(1,1) = 1.0/memory.kappa_0();
//    Kinv_c(2,2) = 1.0/memory.kappa_0();
//    
//    dKinvdp.Zero();
//    dKinvdp(0,0) = -dkappa_ndp/(kappa_n*kappa_n);
//    dKinvdp(1,1) = -dkappa_ndp/(kappa_n*kappa_n);
//    dKinvdp(2,2) = -dkappa_ndp/(kappa_n*kappa_n);
    
    int nphi_q       = datavec[q_b].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    int firstq       = 0;
    int firstp       = nphi_q + firstq;
    
    STATE rho        = m_rho_0 * (1 + (m_c/m_scale_factor)*(p-p_0));
    STATE rho_n      = m_rho_0 * (1 + (m_c/m_scale_factor)*(p_n-p_0));
    STATE drho_ndp_n = m_c/m_scale_factor;
    STATE lambda = rho_n/m_eta;
    
    /// Defining local variables
    TPZFNMatrix<3,STATE> Kl_inv_q(3,1),dK_invdp_q(3,1),Kl_inv_phi_q_j(3,1);
    
    for (int i = 0; i < Dimension(); i++) {
        Kl_inv_q(i,0)     = (1.0/lambda) * (1.0/kappa_n) * q[i];
        dK_invdp_q(i,0)   = (1.0/lambda) * (-dkappa_ndp/(kappa_n*kappa_n)) * q[i];
    }
    
    /// Integration point contribution
    TPZFNMatrix<3,STATE> phi_q_i(3,1), phi_q_j(3,1);
    
    if(!m_simulation_data->IsCurrentStateQ()){
        return;
    }
    
    int s_i, s_j;
    int v_i, v_j;
    
    for (int iq = 0; iq < nphi_q; iq++)
    {
        v_i = datavec[q_b].fVecShapeIndex[iq].first;
        s_i = datavec[q_b].fVecShapeIndex[iq].second;
        
        STATE Kl_inv_dot_q = 0.0;
        STATE dK_invdp_dot_q = 0.0;
        for (int i = 0; i < Dimension(); i++)
        {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[q_b].fNormalVec(i,v_i);
            Kl_inv_dot_q    += Kl_inv_q(i,0)*phi_q_i(i,0);
            dK_invdp_dot_q    += dK_invdp_q(i,0)*phi_q_i(i,0);
        }
        
        ef(iq + firstq) +=  weight * ( m_scale_factor * Kl_inv_dot_q - (p_n) * div_phi(iq,0) );
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            
            v_j = datavec[q_b].fVecShapeIndex[jq].first;
            s_j = datavec[q_b].fVecShapeIndex[jq].second;
            
            STATE Kl_inv_phi_q_j_dot_phi_q_j = 0.0;
            for (int j = 0; j < Dimension(); j++)
            {
                phi_q_j(j,0) = phi_qs(s_j,0) * datavec[q_b].fNormalVec(j,v_j);
                Kl_inv_phi_q_j(j,0) = (1.0/lambda) * (1.0/kappa_n) * phi_q_j(j,0);
                Kl_inv_phi_q_j_dot_phi_q_j += Kl_inv_phi_q_j(j,0)*phi_q_i(j,0);
            }
            
            
            ek(iq + firstq,jq + firstq) += m_scale_factor * weight * Kl_inv_phi_q_j_dot_phi_q_j;
        }
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(iq + firstq, jp + firstp) += weight * (m_scale_factor * dK_invdp_dot_q - div_phi(iq,0)) * phi_ps(jp,0);
        }
        
    }
    
    STATE current_div_q = datavec[q_b].divsol[0][0];
    STATE last_div_q = memory.f();
    STATE div_q = m_theta_scheme*current_div_q+(1.0-m_theta_scheme)*last_div_q;
    
    for (int ip = 0; ip < nphi_p; ip++)
    {
        
        ef(ip + firstp) += -1.0 * weight * (div_q + (1.0/dt) * ( phi_n*rho_n - phi*rho )) * phi_ps(ip,0);
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            ek(ip + firstp, jq + firstq) += -1.0 * weight * m_theta_scheme * div_phi(jq,0) * phi_ps(ip,0);
        }
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(ip + firstp, jp + firstp) += -1.0 * weight * ( (1.0/dt) * ((phi_n * drho_ndp_n + dphi_ndp * rho_n) * phi_ps(jp,0) ) ) * phi_ps(ip,0);
        }
        
    }
    
    if (m_simulation_data->GetTransferCurrentToLastQ()) {
        
        /// porosity update
        {
            
            REAL alpha = this->MemItem(gp_index).Alpha();
            REAL Kdr   = this->MemItem(gp_index).Kdr();
            REAL phi_0 = this->MemItem(gp_index).phi_0();
            
            REAL p_0   = this->MemItem(gp_index).p_0();
            REAL p_n   = this->MemItem(gp_index).p_n();
            
            REAL S = (1.0-alpha)*(alpha-phi_0)/Kdr;
            
            REAL eps_v_t_0   = this->MemItem(gp_index).GetPlasticState_0().m_eps_t.I1();
            REAL eps_v_t_n   = this->MemItem(gp_index).GetPlasticState_n().m_eps_t.I1();
            
            REAL phi_n = phi_0 + (alpha) * (eps_v_t_n-eps_v_t_0) + S * (p_n - p_0);
            this->MemItem(gp_index).Setphi_n(phi_n);
        }
        
        this->MemItem(gp_index).Setf(current_div_q);
        this->MemItem(gp_index).Setp(this->MemItem(gp_index).p_n());
        return;
    }
    
}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    unsigned int q_b = 0;
    unsigned int p_b = 1;
    
    if (m_simulation_data->Get_must_accept_solution_Q()) {
        
        long gp_index = datavec[0].intGlobPtIndex;
        TPZManVector<STATE,3> q  = datavec[q_b].sol[0];
        STATE p                  = datavec[p_b].sol[0][0];
        TPZFMatrix<STATE> dqdx   = datavec[q_b].dsol[0];
        STATE div_q              = datavec[q_b].divsol[0][0];
        
        if (m_simulation_data->IsInitialStateQ()) {
            this->MemItem(gp_index).Setp_0(p);
        }
        
        if (m_simulation_data->IsCurrentStateQ()) {
            this->MemItem(gp_index).Setp_n(p);
            this->MemItem(gp_index).Setq_n(q);
            this->MemItem(gp_index).Setdiv_q_n(div_q);
        }else{
            this->MemItem(gp_index).Setp(p);
            this->MemItem(gp_index).Setq(q);
        }
        
    }

    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
    return;
    
}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    /// Undarined contribute
    if(m_simulation_data->IsInitialStateQ()){
        return;
    }
    
    int q_b = 0;
    
    TPZFNMatrix<100,STATE> phi_qs       = datavec[q_b].phi;
    
    int nphi_q       = phi_qs.Rows();
    int first_q      = 0;
    
    TPZManVector<REAL,3> q  = datavec[q_b].sol[0];
    REAL BigNumber = TPZDiscontinuousGalerkin::gBigNumber;
    REAL Value = bc.Val2()(0,0);
    
    switch (bc.Type())
    {
            
        case 0 :    /// Dirichlet BC  PD
        {
            STATE p_D = Value;
            for (int iq = 0; iq < nphi_q; iq++)
            {
                ef(iq + first_q) += weight * p_D * phi_qs(iq,0);
            }
        }

            break;
            
        case 1 :    /// Neumann BC  QN
        {
            
            for (int iq = 0; iq < nphi_q; iq++)
            {
                STATE qn_N = Value, qn = q[0];
                ef(iq + first_q) += weight * BigNumber * (qn - qn_N) * phi_qs(iq,0);
                
                for (int jq = 0; jq < nphi_q; jq++)
                {
                    
                    ek(iq + first_q,jq + first_q) += weight * BigNumber * phi_qs(jq,0) * phi_qs(iq,0);
                }
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
int TPMRSMonoPhasic<TMEM>::VariableIndex(const std::string &name) {
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
int TPMRSMonoPhasic<TMEM>::NSolutionVariables(int var) {
    switch(var)
    {
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
void TPMRSMonoPhasic<TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout){
    
    long gp_index = data.intGlobPtIndex;
    TMEM & memory = this->MemItem(gp_index); 
    Solout.Resize( this->NSolutionVariables(var));
    
    TPZManVector<STATE,3> qb  = memory.q_n();
    
    
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
            Solout[0] = memory.div_q_n();            
        }
            break;
        case 4:
        {
            Solout[0] = qb[0];

        }
            break;
        case 5:
        {
            Solout[0] = qb[1];
        }
            break;
        case 6:
        {
            Solout[0] = qb[2];
        }
            break;
        case 7:
        {
            for (int i = 0; i < m_dimension; i++)
            {
                Solout[i] = qb[i];
            }
            
        }
            break;
        default:
        {
            std::cout << "TPMRSMonoPhasic<TMEM>:: Variable not implemented." << std::endl;
            DebugStop();
        }
            break;
    }

}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::porosity(long gp_index, REAL &phi_n, REAL &dphi_ndp, REAL &phi){
    
    TMEM & memory = this->MemItem(gp_index);
    
    REAL alpha = memory.Alpha();
    REAL Kdr   = memory.Kdr();
    REAL phi_0 = memory.phi_0();
    
    REAL p_0   = memory.p_0();
    REAL p     = memory.p();
    REAL p_n   = memory.p_n();
    
    REAL epsilon_p_v_0 = (this->MemItem(gp_index).GetPlasticState_0().m_eps_p.I1()/3);
    REAL epsilon_p_v   = (this->MemItem(gp_index).GetPlasticState().m_eps_p.I1()/3);
    REAL phi_p_0 = alpha * epsilon_p_v_0;
    REAL phi_p   = alpha * epsilon_p_v;
    
    REAL sigma_t_v_0 = (memory.GetSigma_0().I1()/3) - alpha * p_0;
    REAL sigma_t_v   = (memory.GetSigma().I1()/3)  - alpha * p;
    REAL geo_delta_phi   = (alpha/Kdr)*(sigma_t_v-sigma_t_v_0) + (phi_p - phi_p_0);
    REAL geo_delta_phi_n = memory.delta_phi() + (alpha/Kdr)*(sigma_t_v-sigma_t_v_0) + (phi_p - phi_p_0);
    
    m_phi_model.Porosity(phi, dphi_ndp, phi_0, p, p_0, alpha, Kdr, geo_delta_phi);
    m_phi_model.Porosity(phi_n, dphi_ndp, phi_0, p_n, p_0, alpha, Kdr, geo_delta_phi_n);
    this->MemItem(gp_index).Setphi_n(phi_n);
    
}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::permeability(long gp_index, REAL &kappa_n, REAL &dkappa_ndphi,REAL &phi,REAL &phi_0){
    
    TMEM & memory = this->MemItem(gp_index);
    REAL kappa_0  = memory.kappa_0();
    m_kappa_model.Permeability(kappa_n, dkappa_ndphi, kappa_0, phi, phi_0);
    this->MemItem(gp_index).Setkappa_n(kappa_n);
}
