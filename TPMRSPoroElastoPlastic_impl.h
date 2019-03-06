//
//  TPMRSPoroElastoPlastic.h
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#include "TPMRSPoroElastoPlastic.h"


template <class T, class TMEM>
TPMRSPoroElastoPlastic<T,TMEM>::TPMRSPoroElastoPlastic() : TPZMatWithMem<TMEM>() {
    
    m_simulation_data = NULL;
    m_dimension         = 0;
    m_c                 = 0;
    m_eta               = 0;
    m_rho_0             = 0;
    m_scale_factor      = 1;
    m_theta_scheme      = 1;
}

template <class T, class TMEM>
TPMRSPoroElastoPlastic<T,TMEM>::TPMRSPoroElastoPlastic(int matid) : TPZMatWithMem<TMEM>(matid) {
    
    m_simulation_data = NULL;
    m_dimension         = 0;
    m_c                 = 0;
    m_eta               = 0;
    m_rho_0             = 0;
    m_scale_factor      = 1;
    m_theta_scheme      = 1;
}

template <class T, class TMEM>
TPMRSPoroElastoPlastic<T,TMEM>::~TPMRSPoroElastoPlastic(){
    
}

template <class T, class TMEM>
TPMRSPoroElastoPlastic<T,TMEM>::TPMRSPoroElastoPlastic(const TPMRSPoroElastoPlastic& other): TPZMatWithMem<TMEM>(other){
    m_simulation_data   = other.m_simulation_data;
    m_dimension         = other.m_dimension;
    m_plastic_integrator    = other.m_plastic_integrator;
    m_c                 = other.m_c;
    m_eta               = other.m_eta;
    m_rho_0             = other.m_rho_0;
    m_scale_factor      = other.m_scale_factor;
    m_phi_model         = other.m_phi_model;
    m_kappa_model       = other.m_kappa_model;
    m_theta_scheme      = other.m_theta_scheme;
}

template <class T, class TMEM>
TPMRSPoroElastoPlastic<T,TMEM> & TPMRSPoroElastoPlastic<T,TMEM>::TPMRSPoroElastoPlastic<T,TMEM>::operator = (const TPMRSPoroElastoPlastic& other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_simulation_data   = other.m_simulation_data;
    m_dimension         = other.m_dimension;
    m_plastic_integrator    = other.m_plastic_integrator;
    m_c                 = other.m_c;
    m_eta               = other.m_eta;
    m_rho_0             = other.m_rho_0;
    m_scale_factor      = other.m_scale_factor;
    m_phi_model         = other.m_phi_model;
    m_kappa_model       = other.m_kappa_model;
    m_theta_scheme      = other.m_theta_scheme;
    return *this;
    
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Print(std::ostream & out){
    
    out << " Material name : " << this->Name() << "\n";
    out << " Pointer to TPMRSSimulationData : " << m_simulation_data << "\n";
    out << " Material dimension : " << m_dimension << "\n";
    m_plastic_integrator.Print(out);
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

template <class T, class TMEM>
int TPMRSPoroElastoPlastic<T,TMEM>::NStateVariables(){
    return 1;
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::FillDataRequirements(TPZVec<TPZMaterialData > &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    DebugStop();
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    DebugStop();
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    DebugStop();
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    DebugStop();
}

template <class T, class TMEM>
int TPMRSPoroElastoPlastic<T,TMEM>::VariableIndex(const std::string &name){
    if (!strcmp("u"   , name.c_str())) return 0;
    if (!strcmp("s"   , name.c_str())) return 1;
    if (!strcmp("e"   , name.c_str())) return 2;
    if (!strcmp("ep"  , name.c_str())) return 3;
    if (!strcmp("s_t" , name.c_str())) return 4;
    if (!strcmp("p"     , name.c_str())) return 5;
    if (!strcmp("phi"   , name.c_str())) return 6;
    if (!strcmp("kappa" , name.c_str())) return 7;
    if (!strcmp("q"    , name.c_str()))  return 8;
    return TPZMatWithMem<TMEM>::VariableIndex(name);
}

template <class T, class TMEM>
int TPMRSPoroElastoPlastic<T,TMEM>::NSolutionVariables(int var){
    switch(var) {
        case 0:
            return m_dimension; /// Vector
        case 1:
            return 9; /// Tensor
        case 2:
            return 9; /// Tensor
        case 3:
            return 9; /// Tensor
        case 4:
            return 9; /// Tensor
        case 5:
            return 1; /// Scalar
        case 6:
            return 1; /// Scalar
        case 7:
            return 1; /// Scalar
        case 8:
            return m_dimension; /// Vector
            
    }
    return TPZMatWithMem<TMEM>::NSolutionVariables(var);
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    long gp_index = data.intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    Solout.Resize( this->NSolutionVariables(var));
    TPZTensor<REAL> epsilon_t = memory.GetPlasticState_n().m_eps_t;
    TPZTensor<REAL> epsilon_p = memory.GetPlasticState_n().m_eps_p;
    
    switch (var) {
        case 0:
        {
            for(int i = 0; i < m_dimension; i++){
                Solout[i] = memory.Getu_n()[i];
            }
        }
            break;
        case 1:
        {
            Solout[0] = memory.GetSigma_n().XX();
            Solout[1] = Solout[3] = memory.GetSigma_n().XY();
            Solout[2] = Solout[6] = memory.GetSigma_n().XZ();
            Solout[4] = memory.GetSigma_n().YY();
            Solout[5] = Solout[7] = memory.GetSigma_n().YZ();
            Solout[8] = memory.GetSigma_n().ZZ();
        }
            break;
        case 2:
        {
            Solout[0] = epsilon_t.XX();
            Solout[1] = Solout[3] = epsilon_t.XY();
            Solout[2] = Solout[6] = epsilon_t.XZ();
            Solout[4] = epsilon_t.YY();
            Solout[5] = Solout[7] = epsilon_t.YZ();
            Solout[8] = epsilon_t.ZZ();
        }
            break;
        case 3:
        {
            Solout[0] = epsilon_p.XX();
            Solout[1] = Solout[3] = epsilon_p.XY();
            Solout[2] = Solout[6] = epsilon_p.XZ();
            Solout[4] = epsilon_p.YY();
            Solout[5] = Solout[7] = epsilon_p.YZ();
            Solout[8] = epsilon_p.ZZ();
        }
            break;
        case 4:
        {
            Solout[0] = memory.GetSigma_n().XX() - memory.Alpha()*memory.p_n();
            Solout[1] = Solout[3] = memory.GetSigma_n().XY();
            Solout[2] = Solout[6] = memory.GetSigma_n().XZ();
            Solout[4] = memory.GetSigma_n().YY() - memory.Alpha()*memory.p_n();
            Solout[5] = Solout[7] = memory.GetSigma_n().YZ();
            Solout[8] = memory.GetSigma_n().ZZ() - memory.Alpha()*memory.p_n();
        }
            break;
        case 5:
        {
            Solout[0] = memory.p_n();
        }
            break;
        case 6:
        {
            Solout[0] = memory.phi_n();
        }
            break;
        case 7:
        {
            Solout[0] = memory.kappa_n();
        }
            break;
        case 8:
        {
            for (int i = 0; i < m_dimension; i++)
            {
                Solout[i] = memory.q_n()[i];
            }
        }
            break;
        default:
        {
            std::cout << "TPMRSPoroElastoPlastic<T,TMEM>:: Variable not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
    
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Epsilon(TPZMaterialData &data, TPZTensor<REAL> & epsilon_t){
    DebugStop();
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Sigma(TPZMaterialData &data, TPZTensor<REAL> & epsilon_t, TPZTensor<REAL> & sigma, TPZFMatrix<REAL> * Dep){
    DebugStop();
}
