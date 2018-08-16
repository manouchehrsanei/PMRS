//
//  TPZPMRSPoroPermMemory.cpp
//  PMRS
//
//  Created by Omar Durán on 8/16/18.
//

#include "TPZPMRSPoroPermMemory.h"


TPZPMRSPoroPermMemory::TPZPMRSPoroPermMemory() : TPZElastoPlasticMem(){
    m_pressure = 0;
    m_kappa = 0;
    m_porosity = 0;
}

/// Copy constructor
TPZPMRSPoroPermMemory::TPZPMRSPoroPermMemory(const TPZPMRSPoroPermMemory & other){
    
    if(&other != this){
        fSigma          = other.fSigma;
        fPlasticState   = other.fPlasticState;
        fPlasticSteps   = other.fPlasticSteps;
        fPhi            = other.fPhi;
        fDisplacement   = other.fDisplacement;
        m_pressure      = other.m_pressure;
        m_kappa         = other.m_kappa;
        m_porosity      = other.m_porosity;
    }

}

/// Assignement constructor
const TPZPMRSPoroPermMemory & TPZPMRSPoroPermMemory::operator=(const TPZPMRSPoroPermMemory & other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    fSigma          = other.fSigma;
    fPlasticState   = other.fPlasticState;
    fPlasticSteps   = other.fPlasticSteps;
    fPhi            = other.fPhi;
    fDisplacement   = other.fDisplacement;
    m_pressure      = other.m_pressure;
    m_kappa         = other.m_kappa;
    m_porosity      = other.m_porosity;
    return *this;
    
}

/// Desconstructor
TPZPMRSPoroPermMemory::~TPZPMRSPoroPermMemory(){
    
}

/// Class name
const std::string TPZPMRSPoroPermMemory::Name()const{
    return "TPZPMRSPoroPermMemory";
}

/// Write class attributes
void TPZPMRSPoroPermMemory::Write(TPZStream &buf, int withclassid) const{
    TPZElastoPlasticMem::Write(buf,withclassid);
    buf.Write(&m_pressure);
    buf.Write(&m_kappa);
    buf.Write(&m_porosity);
}

/// Read class attributes
void TPZPMRSPoroPermMemory::Read(TPZStream &buf, void *context){
    TPZElastoPlasticMem::Read(buf,context);
    buf.Read(&m_pressure);
    buf.Read(&m_kappa);
    buf.Read(&m_porosity);
}

/// Print class attributes
void TPZPMRSPoroPermMemory::Print(std::ostream &out) const{
    
    out << Name();
    out << "\n Sigma = " << fSigma;
    out << "\n Plastic State = " << fPlasticState;
    out << "\n Plastic steps = " << fPlasticSteps;
    out << "\n Displacement = " << fDisplacement;
    out << "\n Phi = " << fPhi;
    out << "\n Pore Pressure = " << m_pressure;
    out << "\n Kappa = " << m_kappa;
    out << "\n Porosity = " << m_porosity;
}


int TPZPMRSPoroPermMemory::ClassId() const{
    return Hash("TPZPMRSPoroPermMemory");
}
