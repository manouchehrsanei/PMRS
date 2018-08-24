//
//  TPZDarcyMem.cpp
//  PZ
//
//  Created by Manouchehr on Agust 23, 2018.
//
//

#include "TPZDarcyMem.h"
#include <stdio.h>
#include "pzreal.h"
#include "pzfmatrix.h"

TPZDarcyMem::TPZDarcyMem()
{
    m_pressure      = 0;
    m_dPorePressure = 0;
    m_kappa         = 0;
}

/// Copy constructor
TPZDarcyMem::TPZDarcyMem(const TPZDarcyMem & other)
{
    if(&other != this)
    {
        m_pressure        = other.m_pressure;
        m_dPorePressure   = other.m_dPorePressure;
        m_kappa           = other.m_kappa;
    }
}

/// Assignement constructor
const TPZDarcyMem & TPZDarcyMem::operator=(const TPZDarcyMem & other)
{
    // check for self-assignment
    if(&other == this)
    {
        return *this;
    }
    m_pressure        = other.m_pressure;
    m_dPorePressure   = other.m_dPorePressure;
    m_kappa           = other.m_kappa;
    return *this;
}

/// Desconstructor
TPZDarcyMem::~TPZDarcyMem()
{
    
}

/// Class name
const std::string TPZDarcyMem::Name()const
{
    return "TPZDarcyMem";
}

/// Write class attributes
void TPZDarcyMem::Write(TPZStream &buf, int withclassid) const
{
    buf.Write(&m_pressure);
    buf.Write(&m_dPorePressure);
    buf.Write(&m_kappa);
}

/// Read class attributes
void TPZDarcyMem::Read(TPZStream &buf, void *context)
{
    buf.Read(&m_pressure);
    buf.Write(&m_dPorePressure);
    buf.Read(&m_kappa);
}

/// Print class attributes
void TPZDarcyMem::Print(std::ostream &out) const
{
    out << Name();
    out << "\n Pore Pressure             = " << m_pressure;
    out << "\n Gradient of Pore Pressure = " << m_dPorePressure;
    out << "\n Absolute Permeability     = " << m_kappa;
}


int TPZDarcyMem::ClassId() const{
    return Hash("TPZDarcyMem");
}
