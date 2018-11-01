//
//  TPMRSUndrainedParameters.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 10/14/18.
//

#include "TPMRSUndrainedParameters.h"


TPMRSUndrainedParameters::TPMRSUndrainedParameters(){
    m_parameters.resize(0);
}

TPMRSUndrainedParameters::TPMRSUndrainedParameters(const TPMRSUndrainedParameters & other){
    
    m_parameters = other.m_parameters;
}

const TPMRSUndrainedParameters & TPMRSUndrainedParameters::operator=(const TPMRSUndrainedParameters & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_parameters = other.m_parameters;
    return *this;
}

TPMRSUndrainedParameters::~TPMRSUndrainedParameters(){
    
}

const std::string TPMRSUndrainedParameters::Name() const {
    return "TPMRSUndrainedParameters";
}

void TPMRSUndrainedParameters::Print(std::ostream &out) const {
    out << Name() << std::endl;
    for (int i  = 0; i < m_parameters.size(); i++)
    {
        out << "parameter number " << i << " = " << m_parameters[i] << std::endl;
    }
}


void TPMRSUndrainedParameters::SetParameters(std::vector<REAL> & parameters){
    if (parameters.size() != 4)
    {
        std::cout << "TPMRSUndrainedParameters:: It is required 4 parameters." << std::endl;
        DebugStop();
    }
    m_parameters = parameters;
}

std::vector<REAL> & TPMRSUndrainedParameters::GetParameters(){
    return m_parameters;
}
