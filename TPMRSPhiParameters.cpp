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
    Initialize();
}

TPMRSPhiParameters::TPMRSPhiParameters(const TPMRSPhiParameters & other){
    m_model         = other.m_model;
    m_parameters    = other.m_parameters;
    Initialize();
}


const TPMRSPhiParameters & TPMRSPhiParameters::operator=(const TPMRSPhiParameters & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_model         = other.m_model;
    m_parameters    = other.m_parameters;
    Initialize();
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
    
}

void TPMRSPhiParameters::ConfigurateModel(std::string model, std::vector<REAL> & parameters){
        
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
    m_parameters = parameters;
}

void TPMRSPhiParameters::Porosity(REAL &phi, REAL &dphi_dp, REAL &phi_0, REAL &p, REAL &p_0, REAL &eps_v, REAL &eps_v_0, REAL &alpha, REAL &Se){
    
    switch (m_model)
    {
        case p_constant : {
            phi = phi_0;
            dphi_dp = 0.0;
        }
            break;
        case p_linear : {
            phi = phi_0 + alpha * (eps_v-eps_v_0) + Se * (p - p_0);
            dphi_dp = Se;
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
    
//    std::cout << "TPMRSPhiParameters::Initialization." << std::endl;
}


void TPMRSPhiParameters::SetParameters(std::vector<REAL> parameters){
    m_parameters = parameters;
}

std::vector<REAL> & TPMRSPhiParameters::GetParameters(){
    return m_parameters;
}

EPhiModel TPMRSPhiParameters::GetModel(){
    return m_model;
}

