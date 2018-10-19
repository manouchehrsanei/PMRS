//
//  TPMRSPoroMechParameters.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 10/14/18.
//

#include "TPMRSPoroMechParameters.h"

TPMRSPoroMechParameters::TPMRSPoroMechParameters(){
    m_parameters.resize(0);
}


TPMRSPoroMechParameters::TPMRSPoroMechParameters(const TPMRSPoroMechParameters & other){
    
    m_parameters = other.m_parameters;
}

const TPMRSPoroMechParameters & TPMRSPoroMechParameters::operator=(const TPMRSPoroMechParameters & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_parameters = other.m_parameters;
    
    return *this;
}


TPMRSPoroMechParameters::~TPMRSPoroMechParameters(){
    
}

const std::string  TPMRSPoroMechParameters::Name() const{
    return "TPMRSPoroMechParameters";
}

void TPMRSPoroMechParameters::Print(std::ostream &out) const{
    out << Name() << std::endl;
    for (int i  = 0; i < m_parameters.size(); i++) {
        out << "parameter number " << i << " = " << m_parameters[i] << std::endl;
    }
}


void TPMRSPoroMechParameters::SetParameters(std::vector<REAL> & parameters){
    if (parameters.size() != 7) {
        std::cout << "TPMRSUndrainedParameters:: It is required 7 parameters." << std::endl;
        DebugStop();
    }
    m_parameters = parameters;
}

std::vector<REAL> & TPMRSPoroMechParameters::GetParameters(){
    return m_parameters;
}
