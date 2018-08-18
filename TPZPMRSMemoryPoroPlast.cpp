//
//  TPZPMRSMemoryPoroPlast.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/6/16.
//

#include "TPZPMRSMemoryPoroPlast.h"


TPZPMRSMemoryPoroPlast::TPZPMRSMemoryPoroPlast() : TPZElastoPlasticMem(){
    m_pressure = 0;
    m_kappa = 0;
    m_porosity = 0;
}

/// Copy constructor
TPZPMRSMemoryPoroPlast::TPZPMRSMemoryPoroPlast(const TPZPMRSMemoryPoroPlast & other){
    
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
const TPZPMRSMemoryPoroPlast & TPZPMRSMemoryPoroPlast::operator=(const TPZPMRSMemoryPoroPlast & other){
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
TPZPMRSMemoryPoroPlast::~TPZPMRSMemoryPoroPlast(){
    
}

/// Class name
const std::string TPZPMRSMemoryPoroPlast::Name()const{
    return "TPZPMRSMemoryPoroPlast";
}

/// Write class attributes
void TPZPMRSMemoryPoroPlast::Write(TPZStream &buf, int withclassid) const{
    TPZElastoPlasticMem::Write(buf,withclassid);
    buf.Write(&m_pressure);
    buf.Write(&m_kappa);
    buf.Write(&m_porosity);
}

/// Read class attributes
void TPZPMRSMemoryPoroPlast::Read(TPZStream &buf, void *context){
    TPZElastoPlasticMem::Read(buf,context);
    buf.Read(&m_pressure);
    buf.Read(&m_kappa);
    buf.Read(&m_porosity);
}

/// Print class attributes
void TPZPMRSMemoryPoroPlast::Print(std::ostream &out) const{
    
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


int TPZPMRSMemoryPoroPlast::ClassId() const{
    return Hash("TPZPMRSMemoryPoroPlast");
}
