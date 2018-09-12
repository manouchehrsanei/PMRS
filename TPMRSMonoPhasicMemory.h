//
//  TPMRSMonoPhasicMemory.h
//  PMRS
//
//  Created by Omar Durán on 9/9/18.
//

#ifndef TPMRSMonoPhasicMemory_h
#define TPMRSMonoPhasicMemory_h

#include <stdio.h>
#include <iostream>
#include <string>
#include "TPZSavable.h"
#include "Hash/TPZHash.h"
#include "TPZStream.h"
#include "TPZPersistenceManager.h"

class TPMRSMonoPhasicMemory{
    
private:
    
    /// pore pressure at initial state
    STATE m_p_0;
    
    /// pore pressure
    STATE m_p;
    
    /// pore pressure
    STATE m_p_n;
    
    /// absolute permeability at initial state
    STATE m_kappa_0;
    
    /// absolute permeability
    STATE m_kappa;
    
    /// lagrangian porosity at intial state
    STATE m_phi_0;
    
    /// lagrangian porosity
    STATE m_phi;
    
public:
    
    /// Default constructor
    TPMRSMonoPhasicMemory();
    
    /// Copy constructor
    TPMRSMonoPhasicMemory(const TPMRSMonoPhasicMemory & other);
    
    /// Assignement constructor
    const TPMRSMonoPhasicMemory & operator=(const TPMRSMonoPhasicMemory & other);
    
    /// Desconstructor
    virtual ~TPMRSMonoPhasicMemory();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPMRSMonoPhasicMemory & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;
    
    
    /// set and get methods
    
    /// Set pore pressure at initial state
    void Setp_0(STATE p_0)
    {
        m_p_0 = p_0;
    }
    
    /// Get pore pressure at initial state
    STATE p_0()
    {
        return m_p_0;
    }
    
    
    /// Set pore pressure at current state
    void Setp(STATE p)
    {
        m_p = p;
    }
    
    /// Get pore pressure at current state
    STATE p()
    {
        return m_p;
    }
    
    
    /// Set pore pressure at last state
    void Setp_n(STATE p_n)
    {
        m_p_n = p_n;
    }
    
    /// Get pore pressure at last state
    STATE p_n()
    {
        return m_p_n;
    }

    
    /// Set absolute permeability at initial state
    void Setkappa_0(STATE kappa_0)
    {
        m_kappa_0 = kappa_0;
    }
    
    /// Get absolute permeability at initial state
    STATE kappa_0()
    {
        return m_kappa_0;
    }

    
    /// Set absolute permeability at current state
    void Setkappa(STATE kappa)
    {
        m_kappa = kappa;
    }
    
    /// Get absolute permeability at current state
    STATE kappa()
    {
        return m_kappa;
    }

    
    /// Set lagrangian porosity at intial state
    void Setphi_0(STATE phi_0)
    {
        m_phi_0 = phi_0;
    }
    
    /// Get lagrangian porosity at intial state
    STATE phi_0()
    {
        return m_phi_0;
    }
    
    
    /// Set lagrangian porosity at current state
    void Setphi(STATE phi)
    {
        m_phi = phi;
    }
    
    /// Get lagrangian porosity at current state
    STATE phi()
    {
        return m_phi;
    }
    
    
    
};


#endif /* TPMRSMonoPhasicMemory_h */
