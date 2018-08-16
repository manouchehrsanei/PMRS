//
//  TPZPMRSPoroPermMemory.h
//  PMRS
//
//  Created by Omar Durán on 8/16/18.
//

#ifndef TPZPMRSPoroPermMemory_h
#define TPZPMRSPoroPermMemory_h

#include <stdio.h>
#include "TPZElastoPlasticMem.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Non linear poro mechanich memory
/////////////////////////////////////////////////////////////////////////////////////////////////////////

class TPZPMRSPoroPermMemory : public TPZElastoPlasticMem
{
    
public:
    
    /// Default constructor
    TPZPMRSPoroPermMemory();
    
    /// Copy constructor
    TPZPMRSPoroPermMemory(const TPZPMRSPoroPermMemory & other);
    
    /// Assignement constructor
    const TPZPMRSPoroPermMemory & operator=(const TPZPMRSPoroPermMemory & other);
    
    /// Desconstructor
    virtual ~TPZPMRSPoroPermMemory();
    
    /// Class name
    const std::string Name()const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout)const;
    
    /// Print class attributes
    friend std::ostream& operator<<( std::ostream& out, const TPZPMRSPoroPermMemory & s ){
        s.Print(out);
        return out;
    }
    
    /// pore pressure
    STATE m_pressure;
    
    /// absolute permeability
    STATE m_kappa;
    
    /// lagrangian porosity
    STATE m_porosity;
    
    virtual int ClassId() const;
    
};

#endif /* TPZPMRSPoroPermMemory_h */
