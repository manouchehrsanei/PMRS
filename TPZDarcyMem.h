//
//  TPZDarcyMem.cpp
//  PZ
//
//  Created by Manouchehr on Agust 23, 2018.
//
//

#ifndef TPZDarcyMem_h
#define TPZDarcyMem_h

#include <stdio.h>
#include "pzreal.h"
#include "pzfmatrix.h"

////////////////////////////////////////////////////////////
//  TPZDarcyMem
////////////////////////////////////////////////////////////

class TPZDarcyMem 
{
    
public:
    
    /// @brief Default constructor
    TPZDarcyMem();
    
    /// @brief Copy constructor
    TPZDarcyMem(const TPZDarcyMem & other);
    
    /// @brief Assignement constructor
    const TPZDarcyMem & operator=(const TPZDarcyMem & other);
    
    /// @brief Desconstructor
    virtual ~TPZDarcyMem();
    
    /// @brief Class name
    const std::string Name()const;
    
    /// @brief Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// @brief Read class attributes
    void Read(TPZStream &buf, void *context);
    
    /// @brief Print class attributes
    virtual void Print(std::ostream &out = std::cout)const;
    
    /// @brief Print class attributes
    friend std::ostream& operator<<( std::ostream& out, const TPZDarcyMem & s )
    {
        s.Print(out);
        return out;
    }
    
    /// Pore Pressure
    STATE m_pressure;
    
    /// Gradient of Pore Pressure
    TPZVec<STATE> m_dPorePressure;
    
    /// Absolute Permeability
    STATE m_kappa;
    
    
    virtual int ClassId() const;
    
};

#endif /* TPZDarcyMem_h */
