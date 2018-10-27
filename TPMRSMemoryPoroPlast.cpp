//
//  TPMRSMemoryPoroPlast.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/6/16.
//

#include "TPMRSMemoryPoroPlast.h"


TPMRSMemoryPoroPlast::TPMRSMemoryPoroPlast() : TPZElastoPlasticMem(){
    m_pressure = 0;
    m_kappa = 0;
    m_porosity = 0;
}

/// Copy constructor
TPMRSMemoryPoroPlast::TPMRSMemoryPoroPlast(const TPMRSMemoryPoroPlast & other){
    
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
const TPMRSMemoryPoroPlast & TPMRSMemoryPoroPlast::operator=(const TPMRSMemoryPoroPlast & other){
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
TPMRSMemoryPoroPlast::~TPMRSMemoryPoroPlast(){
    
}

/// Class name
const std::string TPMRSMemoryPoroPlast::Name()const{
    return "TPMRSMemoryPoroPlast";
}

/// Write class attributes
void TPMRSMemoryPoroPlast::Write(TPZStream &buf, int withclassid) const{
    TPZElastoPlasticMem::Write(buf,withclassid);
    buf.Write(&m_pressure);
    buf.Write(&m_kappa);
    buf.Write(&m_porosity);
}

/// Read class attributes
void TPMRSMemoryPoroPlast::Read(TPZStream &buf, void *context){
    TPZElastoPlasticMem::Read(buf,context);
    buf.Read(&m_pressure);
    buf.Read(&m_kappa);
    buf.Read(&m_porosity);
}

/// Print class attributes
void TPMRSMemoryPoroPlast::Print(std::ostream &out) const{
    
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


int TPMRSMemoryPoroPlast::ClassId() const{
    return Hash("TPMRSMemoryPoroPlast");
}
