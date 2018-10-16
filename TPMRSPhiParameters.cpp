//
//  TPMRSPhiParameters.cpp
//  PMRS
//
//  Created by Omar Durán on 10/15/18.
//

#include "TPMRSPhiParameters.h"


TPMRSPhiParameters::TPMRSPhiParameters(){
    m_model = e_none;
    m_parameters.resize(0);
}

TPMRSPhiParameters::TPMRSPhiParameters(const TPMRSPhiParameters & other){
    m_model         = other.m_model;
    m_parameters    = other.m_parameters;
}


const TPMRSPhiParameters::TPMRSPhiParameters & operator=(const TPMRSPhiParameters & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_model         = other.m_model;
    m_parameters    = other.m_parameters;
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
    switch (std::stoi(model))
    {
        case std::stoi("Linear") : {
            m_model = e_linear;
        }
            break;
        default : {
            DebugStop();
        }
    }
}


