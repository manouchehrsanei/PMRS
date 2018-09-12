//
//  TPMRSMemory.h
//  PMRS
//
//  Created by Omar Durán on 9/9/18.
//

#ifndef TPMRSMemory_h
#define TPMRSMemory_h

#include <stdio.h>
#include "TPMRSMonoPhasicMemory.h"
#include "TPMRSElastoPlasticMemory.h"

class TPMRSMemory : public TPMRSMonoPhasicMemory, TPMRSElastoPlasticMemory {
    
public:
    
    /// Default constructor
    TPMRSMemory();
    
    /// Copy constructor
    TPMRSMemory(const TPMRSMemory & other);
    
    /// Assignement constructor
    const TPMRSMemory & operator=(const TPMRSMemory & other);
    
    /// Desconstructor
    virtual ~TPMRSMemory();
    
    /// Class name
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
    
};

#endif /* TPMRSMemory_h */
