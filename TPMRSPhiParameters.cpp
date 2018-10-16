//
//  TPMRSPhiParameters.cpp
//  PMRS
//
//  Created by Omar Durán on 10/15/18.
//

#include "TPMRSPhiParameters.h"


TPMRSPhiParameters::TPMRSPhiParameters(){
    m_model = p_none;
    m_parameters.resize(0);
    m_K = 0;
    Initialize();
}

TPMRSPhiParameters::TPMRSPhiParameters(const TPMRSPhiParameters & other){
    m_model         = other.m_model;
    m_parameters    = other.m_parameters;
    m_name_to_p_model = other.m_name_to_p_model;
    m_K = other.m_K;
}


const TPMRSPhiParameters & TPMRSPhiParameters::operator=(const TPMRSPhiParameters & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_model         = other.m_model;
    m_parameters    = other.m_parameters;
    m_name_to_p_model = other.m_name_to_p_model;
    m_K               = other.m_K;
    return *this;
    
}

TPMRSPhiParameters::~TPMRSPhiParameters(){
    
}

const std::string TPMRSPhiParameters::Name() const{
    return "TPMRSPhiParameters";
}

void TPMRSPhiParameters::Print(std::ostream &out) const{
    out << Name() << std::endl;
    out << "m_model = " << m_model << std::endl;
    for (int i  = 0; i < m_parameters.size(); i++) {
        out << "parameter number " << i << " = " << m_parameters[i] << std::endl;
    }
    out << "m_K = " << m_K << std::endl;
}

void TPMRSPhiParameters::SetModel(std::string model){
        
    switch (m_name_to_p_model[model])
    {
        case p_constant : {
            m_model = p_constant;
        }
            break;
        case p_linear : {
            m_model = p_linear;
        }
            break;
        default : {
            DebugStop();
        }
    }
}

void TPMRSPhiParameters::Porosity(REAL &phi, REAL &dphi_dp, REAL &phi_0, REAL &p, REAL &p_0, REAL &sigma_v, REAL &sigma_v_0, REAL &alpha, REAL &Se){
    
    switch (m_model)
    {
        case p_constant : {
            phi = phi_0;
            dphi_dp = 0.0;
        }
            break;
        case p_linear : {
            phi = phi_0 + (alpha/m_K) * (sigma_v-sigma_v_0) + (Se + (alpha*alpha)/m_K) * (p - p_0);
            dphi_dp = (Se + (alpha*alpha)/m_K);
        }
            break;
        default : {
            DebugStop();
        }
    }
    
}

void TPMRSPhiParameters::Initialize()
{
    m_name_to_p_model["none"] = p_none;
    m_name_to_p_model["Constant"] = p_constant;
    m_name_to_p_model["Linear"] = p_linear;

}


void TPMRSPhiParameters::SetParameters(std::vector<REAL> parameters){
    m_parameters = parameters;
}

std::vector<REAL> & TPMRSPhiParameters::GetParameters(){
    return m_parameters;
}

TPMRSPhiParameters::EPhiModel TPMRSPhiParameters::GetModel(){
    return m_model;
}

