//
//  TPMRSMonoPhasic_impl.hpp
//  PMRS
//
//  Created by Omar Durán on 9/12/18.
//

#include "TPMRSMonoPhasic.h"

template <class TMEM>
TPMRSMonoPhasic<TMEM>::TPMRSMonoPhasic(){
    m_simulation_data   = NULL;
    m_dimension         = 0;
    m_c                 = 0;
    m_eta               = 0;
    m_rho_0             = 0;
    m_scale_factor      = 1;
}

template <class TMEM>
TPMRSMonoPhasic<TMEM>::TPMRSMonoPhasic(int mat_id, int dimension) : TPZMatWithMem<TMEM>(mat_id){
    m_dimension = dimension;
}

template <class TMEM>
TPMRSMonoPhasic<TMEM>::TPMRSMonoPhasic(const TPMRSMonoPhasic & other){
    m_simulation_data   = other.m_simulation_data;
    m_dimension         = other.m_dimension;
    m_c                 = other.m_c;
    m_eta               = other.m_eta;
    m_rho_0             = other.m_rho_0;
    m_scale_factor      = other.m_scale_factor;
}

template <class TMEM>
TPMRSMonoPhasic<TMEM> & TPMRSMonoPhasic<TMEM>::operator=(const TPMRSMonoPhasic & other){
    
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
    return *this;
}

template <class TMEM>
TPMRSMonoPhasic<TMEM>::~TPMRSMonoPhasic(){
    
}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::Print(std::ostream &out){
    
    out << " Material name : " << this->Name() << "\n";
    out << " Pointer to TPZSimulationData : " << m_simulation_data << "\n";
    out << " Material dimension : " << m_dimension << "\n";
    out << " Fluid compressibility : " << m_c << "\n";
    out << " Fluid viscosity : " << m_eta << "\n";
    out << " Fluid density : " << m_rho_0 << "\n";
    out << " Scale factor  : " << m_scale_factor << "\n";
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
void TPMRSMonoPhasic<TMEM>::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi)
{
    int ublock = 0;
    int dim = this->Dimension();
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;   // For H1  test functions Q
    TPZFMatrix<STATE> dphiuH1       = datavec[ublock].dphi; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiuH1axes   = datavec[ublock].dphix; // Derivative For H1  test functions
    TPZFNMatrix<9,STATE> gradu = datavec[ublock].dsol[0];
    TPZFNMatrix<9,STATE> graduMaster;
    gradu.Transpose();
    
    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[ublock].axes);
    
    int nphiuHdiv = datavec[ublock].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[ublock].detjac;
    
    TPZFMatrix<STATE> Qaxes = datavec[ublock].axes;
    TPZFMatrix<STATE> QaxesT;
    TPZFMatrix<STATE> Jacobian = datavec[ublock].jacobian;
    TPZFMatrix<STATE> JacobianInverse = datavec[ublock].jacinv;
    
    TPZFMatrix<STATE> GradOfX;
    TPZFMatrix<STATE> GradOfXInverse;
    TPZFMatrix<STATE> VectorOnMaster;
    TPZFMatrix<STATE> VectorOnXYZ(3,1,0.0);
    Qaxes.Transpose(&QaxesT);
    QaxesT.Multiply(Jacobian, GradOfX);
    JacobianInverse.Multiply(Qaxes, GradOfXInverse);
    
    int ivectorindex = 0;
    int ishapeindex = 0;
    
    if (HDivPiola == 1)
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            for (int k = 0; k < dim; k++) {
                VectorOnXYZ(k,0) = datavec[ublock].fNormalVec(k,ivectorindex);
            }
            
            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;
            
            /* Contravariant Piola mapping preserves the divergence */
            for (int k = 0; k < dim; k++) {
                DivergenceofPhi(iq,0) +=  dphiuH1(k,ishapeindex)*VectorOnMaster(k,0);
            }
        }
        
    }
    else
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            /* Computing the divergence for constant jacobian elements */
            for (int k = 0; k < dim; k++) {
                DivergenceofPhi(iq,0) +=  datavec[ublock].fNormalVec(k,ivectorindex)*GradphiuH1(k,ishapeindex);
            }
        }
    }
    return;
}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    unsigned int q_b = 0;
    unsigned int p_b = 1;
    
    TPZFNMatrix<100,STATE> phi_qs       = datavec[q_b].phi;
    TPZFNMatrix<100,STATE> phi_ps       = datavec[p_b].phi;
    TPZFNMatrix<10,STATE> grad_q_axes = datavec[q_b].dsol[0];
    TPZFNMatrix<40,STATE> div_on_master;
    
    STATE jac_det;
    this->ComputeDivergenceOnMaster(datavec, div_on_master);
    jac_det = datavec[q_b].detjac;
    
    // Get the pressure at the integrations points
    long gp_index = datavec[0].intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    
    // Time
    STATE dt = m_simulation_data->dt();
    TPZManVector<STATE,3> q    = datavec[q_b].sol[0];
    STATE p_n                  = datavec[p_b].sol[0][0];
    
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
    
    int nphi_q       = datavec[q_b].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    int firstq       = 0;
    int firstp       = nphi_q + firstq;
    
    STATE rho = m_rho_0 * (1 + m_c*(p-p_0));
    STATE rho_n = m_rho_0 * (1 + m_c*(p_n-p_0));
    STATE drho_ndp_n = m_c;
    STATE lambda = rho_n/m_eta;
    
    // Defining local variables
    TPZFNMatrix<3,STATE> Kl_inv_q(3,1),Kdldp_inv_q(3,1),Kl_inv_phi_q_j(3,1);
    
    for (int i = 0; i < Dimension(); i++) {
        STATE dot = 0.0;
        for (int j =0; j < Dimension(); j++) {
            dot    += Kinv(i,j)*q[j];
        }
        Kl_inv_q(i,0)     = (1.0/lambda) * dot;
//        Kdldp_inv_q(i,0)  = (-(m_c*m_rho_0)/(m_eta*lambda*lambda)) * dot;
    }
    
    // Integration point contribution
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
        STATE Kdldp_inv_dot_q = 0.0;
        for (int i = 0; i < Dimension(); i++) {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[q_b].fNormalVec(i,v_i);
            Kl_inv_dot_q    += Kl_inv_q(i,0)*phi_q_i(i,0);
            Kdldp_inv_dot_q    += Kdldp_inv_q(i,0)*phi_q_i(i,0);
        }
        
        ef(iq + firstq) +=  weight * ( m_scale_factor * Kl_inv_dot_q - (1.0/jac_det) * (p_n) * div_on_master(iq,0) );
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            
            v_j = datavec[q_b].fVecShapeIndex[jq].first;
            s_j = datavec[q_b].fVecShapeIndex[jq].second;
            
            STATE Kl_inv_phi_q_j_dot_phi_q_j = 0.0;
            for (int j = 0; j < Dimension(); j++) {
                phi_q_j(j,0) = phi_qs(s_j,0) * datavec[q_b].fNormalVec(j,v_j);
                STATE dot = 0.0;
                for (int k = 0; k < Dimension(); k++) {
                    dot += Kinv(j,k)*phi_q_j(k,0);
                }
                Kl_inv_phi_q_j(j,0) = (1.0/lambda) * dot;
                Kl_inv_phi_q_j_dot_phi_q_j += Kl_inv_phi_q_j(j,0)*phi_q_i(j,0);
            }
            
            
            ek(iq + firstq,jq + firstq) += m_scale_factor * weight * Kl_inv_phi_q_j_dot_phi_q_j;
        }
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
//            ek(iq + firstq, jp + firstp) += weight * (m_scale_factor * Kdldp_inv_dot_q - (1.0/jac_det) * div_on_master(iq,0)) * phi_ps(jp,0);
            ek(iq + firstq, jp + firstp) += weight * (- (1.0/jac_det) * div_on_master(iq,0)) * phi_ps(jp,0);
        }
        
    }
    
    STATE div_q = (grad_q_axes(0,0) + grad_q_axes(1,1) + grad_q_axes(2,2));
    STATE phi_n,dphi_ndp,phi;
    this->porosity(gp_index,phi_n,dphi_ndp,phi);
    
    
    for (int ip = 0; ip < nphi_p; ip++)
    {
        
        ef(ip + firstp) += -1.0 * m_scale_factor * weight * (div_q + (1.0/dt) * ( phi_n*rho_n - phi*rho )) * phi_ps(ip,0);
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            ek(ip + firstp, jq + firstq) += -1.0 * m_scale_factor * weight * (1.0/jac_det) * div_on_master(jq,0) * phi_ps(ip,0);
        }
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(ip + firstp, jp + firstp) += -1.0 * m_scale_factor * weight * ( (1.0/dt) * ((phi_n * drho_ndp_n + dphi_ndp * rho_n) * phi_ps(jp,0) ) ) * phi_ps(ip,0);
        }
        
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
        
        if (m_simulation_data->IsInitialStateQ()) {
            this->GetMemory().get()->operator[](gp_index).Setp_0(p);
        }
        
        if (m_simulation_data->IsCurrentStateQ()) {
            this->GetMemory().get()->operator[](gp_index).Setp_n(p);
        }else{
            this->GetMemory().get()->operator[](gp_index).Setp(p);
        }
        
    }else{
        TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
        this->Contribute(datavec, weight, ek_fake, ef);
    }
    
    
    return;
    
    
    
}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!m_simulation_data->IsCurrentStateQ()) {
        return;
    }
    
    int q_b = 0;
    
    TPZFNMatrix<100,STATE> phi_qs       = datavec[q_b].phi;
    
    int nphi_q       = phi_qs.Rows();
    int first_q      = 0;
    
    TPZManVector<REAL,3> q  = datavec[q_b].sol[0];
    REAL BigNumber = TPZDiscontinuousGalerkin::gBigNumber;
    REAL Value = bc.Val2()(0,0);
    
    switch (bc.Type()) {
            
        case 0 :    // Dirichlet BC  PD
        {
            STATE p_D = Value;
            for (int iq = 0; iq < nphi_q; iq++)
            {
                ef(iq + first_q) += m_scale_factor * weight * p_D * phi_qs(iq,0);
            }
        }

            break;
        case 1 :    // Neumann BC  QN
        {
            
            for (int iq = 0; iq < nphi_q; iq++)
            {
                STATE qn_N = Value, qn = q[0];
                ef(iq + first_q) += m_scale_factor * weight * BigNumber * (qn - qn_N) * phi_qs(iq,0);
                
                for (int jq = 0; jq < nphi_q; jq++)
                {
                    
                    ek(iq + first_q,jq + first_q) += m_scale_factor * weight * BigNumber * phi_qs(jq,0) * phi_qs(iq,0);
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
int TPMRSMonoPhasic<TMEM>::NSolutionVariables(int var) {
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
void TPMRSMonoPhasic<TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout){
    
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
            std::cout << "TPMRSMonoPhasic<TMEM>:: Variable not implemented." << std::endl;
            DebugStop();
        }
            break;
    }

        
    
}

template <class TMEM>
void TPMRSMonoPhasic<TMEM>::porosity(long gp_index, REAL &phi_n, REAL &dphi_ndp, REAL &phi){
    
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    
    REAL phi_0 = memory.phi_0();
//    phi = phi_0;
//    phi_n = phi_0;
//    dphi_ndp = 0.0;
//    this->GetMemory().get()->operator[](gp_index).Setphi(phi);
//    return;
    
    REAL nu = m_simulation_data->Get_nu();
    REAL E = m_simulation_data->Get_young();
    REAL Kdr = E/(3.0*(1.0-2.0*nu));
    
    REAL alpha = memory.GetAlpha();
    REAL Se = memory.GetSe();
    
    REAL p_0 = memory.p_0();
    REAL p = memory.p();
    REAL p_n = memory.p_n();

    REAL sigma_v_0 = memory.GetSigma_0().I1()/3;
    REAL sigma_v = memory.GetSigma().I1()/3;
    REAL sigma_v_n = memory.GetSigma_n().I1()/3;
    
    phi     = phi_0 + (Se + (alpha*alpha)/Kdr)*(p-p_0) + (alpha/Kdr)*(sigma_v-sigma_v_0);
    phi_n   = phi_0 + (Se + (alpha*alpha)/Kdr)*(p_n-p_0) + (alpha/Kdr)*(sigma_v_n-sigma_v_0);
    
    dphi_ndp = (Se + (alpha*alpha)/Kdr);
    this->GetMemory().get()->operator[](gp_index).Setphi(phi); // Current phi, please rename it ot phi_n
}
