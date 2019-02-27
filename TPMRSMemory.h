//
//  TPMRSMemory.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#ifndef TPMRSMemory_h
#define TPMRSMemory_h

#include <stdio.h>
#include "TPMRSMonoPhasicMemory.h"
#include "TPMRSElastoPlasticMemory.h"

class TPMRSMemory : public TPMRSMonoPhasicMemory, public TPMRSElastoPlasticMemory {
    
private:
    
    /// Biot-Willis coefficient
    REAL m_alpha;
    
    /// Drained bulk modulus
    REAL m_Kdr;
    
    /// lagrangian porosity correction realted to geomechanics and fss iterations
    STATE m_delta_phi;
    
public:
    
    /// Default constructor
    TPMRSMemory();
    
    /// Copy constructor
    TPMRSMemory(const TPMRSMemory & other);
    
    /// Assignement constructor
    const TPMRSMemory & operator=(const TPMRSMemory & other);
    
    /// Destructor
    virtual ~TPMRSMemory();
    
    /// Class Name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPMRSMemory & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;
    
    void SetAlpha(REAL alpha){
        m_alpha  = alpha;
    }
    
    REAL Alpha(){
        return m_alpha;
    }
    
    void SetKdr(REAL Kdr){
        m_Kdr  = Kdr;
    }
    
    REAL Kdr(){
        return m_Kdr;
    }
    
    void Setdelta_phi(REAL delta_phi){
        m_delta_phi  = delta_phi;
    }
    
    REAL delta_phi(){
        return m_delta_phi;
    }
    
};

#endif /* TPMRSMemory_h */
