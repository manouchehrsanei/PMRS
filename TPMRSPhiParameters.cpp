//
//  TPMRSPhiParameters.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 10/15/18.
//

#include "TPMRSPhiParameters.h"


TPMRSPhiParameters::TPMRSPhiParameters(){
    m_model = p_none;
    m_parameters.resize(0);
    Initialize();
}

TPMRSPhiParameters::TPMRSPhiParameters(const TPMRSPhiParameters & other){
    m_model           = other.m_model;
    m_parameters      = other.m_parameters;
    m_name_to_p_model = other.m_name_to_p_model;
}


const TPMRSPhiParameters & TPMRSPhiParameters::operator=(const TPMRSPhiParameters & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_model           = other.m_model;
    m_parameters      = other.m_parameters;
    m_name_to_p_model = other.m_name_to_p_model;
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
    for (int i  = 0; i < m_parameters.size(); i++)
    {
        out << "parameter number " << i << " = " << m_parameters[i] << std::endl;
    }
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

void TPMRSPhiParameters::Porosity(REAL &phi, REAL &dphi_dp, REAL &phi_0, REAL &p, REAL &p_0, REAL &sigma_v, REAL &sigma_v_0, REAL &alpha, REAL &Kdr){
    
    REAL S = (1.0-alpha)*(alpha-phi_0)/Kdr;
    switch (m_model)
    {
        case p_constant : {
            phi     = phi_0;
            dphi_dp = 0.0;
        }
            break;
        case p_linear : {
            phi     = phi_0 + (alpha/Kdr) * (sigma_v-sigma_v_0) + (S + (alpha*alpha)/Kdr) * (p - p_0);
            dphi_dp = (S + (alpha*alpha)/Kdr);
        }
            break;
        default : {
            DebugStop();
        }
    }
    
}

void TPMRSPhiParameters::Initialize()
{
    m_name_to_p_model["none"]     = p_none;
    m_name_to_p_model["Constant"] = p_constant;
    m_name_to_p_model["Linear"]   = p_linear;

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

