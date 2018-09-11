//
//  TPZPMRSMemoryPoroPlast.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/6/16.
//

#ifndef TPZPMRSMemoryPoroPlast_h
#define TPZPMRSMemoryPoroPlast_h

#include <stdio.h>
#include "TPZElastoPlasticMem.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Non linear poro mechanich memory: It is written for TPZPMRSMemoryPoroPlast
/////////////////////////////////////////////////////////////////////////////////////////////////////////

class TPZPMRSMemoryPoroPlast : public TPZElastoPlasticMem
{
    
public:
    
    /// Default constructor
    TPZPMRSMemoryPoroPlast();
    
    /// Copy constructor
    TPZPMRSMemoryPoroPlast(const TPZPMRSMemoryPoroPlast & other);
    
    /// Assignement constructor
    const TPZPMRSMemoryPoroPlast & operator=(const TPZPMRSMemoryPoroPlast & other);
    
    /// Desconstructor
    virtual ~TPZPMRSMemoryPoroPlast();
    
    /// Class name
    const std::string Name()const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout)const;
    
    /// Print class attributes
    friend std::ostream& operator<<( std::ostream& out, const TPZPMRSMemoryPoroPlast & s ){
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
    
    /**
     * Total (elastic+plastic) stress
     */
    TPZTensor<REAL> fSigma;
    
    /**
     * Plastic state vars
     */
    TPZPlasticState<REAL> fPlasticState;
    
    int fPlasticSteps;
    
    REAL fPhi;
    
    TPZManVector<REAL,3> fDisplacement;
    
};

#endif /* TPZPMRSMemoryPoroPlast_h */
