//
//  TPMRSCoupPoElaMemory.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 9/6/16.
//
//

#include "TPMRSCoupPoElaMemory.h"

/// Brief Default constructor
TPMRSCoupPoElaMemory::TPMRSCoupPoElaMemory()
{
    m_grad_u_n.Resize(3, 3);
    m_grad_u_n.Zero();
    m_epsilon_e_n.Resize(3, 3);
    m_epsilon_e_n.Zero();
    m_epsilon_p_n.Resize(3, 3);
    m_epsilon_p_n.Zero();
}

TPMRSCoupPoElaMemory::TPMRSCoupPoElaMemory(const TPMRSCoupPoElaMemory & other){
    
    m_grad_u_n     = other.m_grad_u_n;
    m_epsilon_e_n  = other.m_epsilon_e_n;
    m_epsilon_p_n  = other.m_epsilon_p_n;

}

const TPMRSCoupPoElaMemory & TPMRSCoupPoElaMemory::operator=(const TPMRSCoupPoElaMemory & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_grad_u_n     = other.m_grad_u_n;
    m_epsilon_e_n  = other.m_epsilon_e_n;
    m_epsilon_p_n  = other.m_epsilon_p_n;
    
    return *this;
}


TPMRSCoupPoElaMemory::~TPMRSCoupPoElaMemory(){
    
}

const std::string TPMRSCoupPoElaMemory::Name() const{
    return "TPMRSCoupPoElaMemory";
}


void TPMRSCoupPoElaMemory::Write(TPZStream &buf, int withclassid) const
{
    m_grad_u_n.Write(buf, withclassid);
    m_epsilon_e_n.Write(buf, withclassid);
    m_epsilon_p_n.Write(buf, withclassid);
}

void TPMRSCoupPoElaMemory::Read(TPZStream &buf, void *context)
{
    m_grad_u_n.Read(buf, context);
    m_epsilon_e_n.Read(buf, context);
    m_epsilon_p_n.Read(buf, context);
    
}

void TPMRSCoupPoElaMemory::Print(std::ostream &out) const
{
    out << Name();
    out << "\n Last Gradient of deformation = " << m_grad_u_n;
    out << "\n Last Elastic strain  =         " << m_epsilon_e_n;
    out << "\n Last Plastic strain  =         " << m_epsilon_p_n;
    
}

int TPMRSCoupPoElaMemory::ClassId() const{
    return Hash("TPMRSCoupPoElaMemory");
}

