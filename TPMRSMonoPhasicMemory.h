//
//  TPMRSMonoPhasicMemory.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
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
    
    /// pore pressure at last state (n)
    STATE m_p;
    
    /// pore pressure at current state (n+1)
    STATE m_p_n;
    
    /// flux at last state (n)
    TPZManVector<STATE,3> m_q;
    
    /// flux at current state (n+1)
    TPZManVector<STATE,3> m_q_n;
    
    /// flux divergence (n+1)
    STATE m_div_q_n;
    
    /// absolute permeability at initial state
    STATE m_kappa_0;
    
    /// absolute permeability at current state (n+1)
    STATE m_kappa_n;
    
    /// lagrangian porosity at intial state
    STATE m_phi_0;
    
    /// lagrangian porosity at current state (n+1)
    STATE m_phi_n;
    
    /// Value function f(q,p) at last state (for CG right hand side and for Mixed div_q) (required for Crank-Nicolson method)
    STATE m_f;
    
    /// Value vector function f(q,p) at last state (for CG right hand side) (required for Crank-Nicolson method)
    TPZManVector<STATE,3> m_f_vec;
    
public:
    
    /// Default constructor
    TPMRSMonoPhasicMemory();
    
    /// Copy constructor
    TPMRSMonoPhasicMemory(const TPMRSMonoPhasicMemory & other);
    
    /// Assignement constructor
    const TPMRSMonoPhasicMemory & operator=(const TPMRSMonoPhasicMemory & other);
    
    /// Destructor
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
    
    
    /// Set pore pressure at last state
    void Setp(STATE p)
    {
        m_p = p;
    }
    
    /// Get pore pressure at last state
    STATE p()
    {
        return m_p;
    }
    
    
    /// Set pore pressure at current state
    void Setp_n(STATE p_n)
    {
        m_p_n = p_n;
    }
    
    /// Get pore pressure at current state
    STATE p_n()
    {
        return m_p_n;
    }

   
    /// Set flux at last state
    void Setq(TPZManVector<STATE,3> & q)
    {
        m_q = q;
    }
    
    /// Get flux at last state
    TPZManVector<STATE,3> & q(){
        return m_q;
    }
    
    /// Set flux at current state
    void Setq_n(TPZManVector<STATE,3> & q_n)
    {
        m_q_n = q_n;
    }
    
    /// Get flux at current state
    TPZManVector<STATE,3> & q_n(){
        return m_q_n;
    }
    
    /// Set flux divergence at current state
    void Setdiv_q_n(STATE div_q_n)
    {
        m_div_q_n = div_q_n;
    }
    
    /// Get flux divergence at current state
    STATE div_q_n()
    {
        return m_div_q_n;
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
    void Setkappa_n(STATE kappa_n)
    {
        m_kappa_n = kappa_n;
    }
    
    /// Get absolute permeability at current state
    STATE kappa_n()
    {
        return m_kappa_n;
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
    
    /// Get lagrangian porosity at current state
    STATE phi_n()
    {
        return m_phi_n;
    }
    
    /// Set lagrangian porosity at current state
    void Setphi_n(STATE phi_n)
    {
        m_phi_n = phi_n;
    }
    
   
    /// Get Value function f(q,p) at last state
    STATE f()
    {
        return m_f;
    }
    
    /// Set Value function f(q,p) at last state
    void Setf(STATE f)
    {
        m_f = f;
    }
    
    /// Get Value vector function f(q,p) at last state
    TPZManVector<STATE,3> & f_vec()
    {
        return m_f_vec;
    }
    
    /// Set Value vector function f(q,p) at last state
    void Setf_vec(TPZManVector<STATE,3> & f_vec)
    {
        m_f_vec = f_vec;
    }
    
    
};


#endif /* TPMRSMonoPhasicMemory_h */
