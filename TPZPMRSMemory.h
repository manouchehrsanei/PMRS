//
//  TPZPMRSMemory.h
//  PZ
//
//  Created by Omar and Manouchehr on 9/6/16.
//
//

#ifndef TPZPMRSMemory_h
#define TPZPMRSMemory_h

#include <stdio.h>
#include "pzreal.h"
#include "pzfmatrix.h"
#include "TPZElastoPlasticMem.h"


class TPZPMRSMemory
{
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory :
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Basis functions
    TPZElastoPlasticMem m_Mem;
    
    
    /** @brief Gradient of deformation at at n (last) state */
    TPZFNMatrix<9,REAL> m_grad_u_n;
    
   /** @brief elastic strain at n (last) state */
    TPZFMatrix<REAL> m_epsilon_e_n;

   /** @brief plastic strain at n (last) state */
    TPZFMatrix<REAL> m_epsilon_p_n;
    
    /** @brief Gradient of pore pressure at n (last) state */
    TPZFMatrix<REAL> m_grad_p_n;
    
    /** @brief Porosity at n (last) state */
    TPZFMatrix<REAL> m_phi_n;
    
    /** @brief Permeability at n (last) state */
    TPZFMatrix<REAL> m_k_n;
    
    
    
public:
    
    /** @brief Default constructor */
    TPZPMRSMemory();
    
    /** @brief Default destructor */
    ~TPZPMRSMemory();
    
    TPZPMRSMemory(const TPZPMRSMemory &copy)
    {
        
        m_grad_u_n    = copy.m_grad_u_n;
        m_epsilon_e_n = copy.m_epsilon_e_n;
        m_epsilon_p_n = copy.m_epsilon_p_n;
        m_grad_p_n    = copy.m_grad_p_n;
        m_phi_n       = copy.m_phi_n;
        m_k_n         = copy.m_k_n;


    }
    
    TPZPMRSMemory &operator=(const TPZPMRSMemory &other)
    {
        
        m_grad_u_n    = other.m_grad_u_n;
        m_epsilon_e_n = other.m_epsilon_e_n;
        m_epsilon_p_n = other.m_epsilon_p_n;
        m_grad_p_n    = other.m_grad_p_n;
        m_phi_n       = other.m_phi_n;
        m_k_n         = other.m_k_n;
        
        return *this;
    }
    
    void UpdateSolutionMemory()
    {
        //update pressure and total flux (un = unp1)
        DebugStop();
    }
    
    
    
    /** @brief Set Gradient of deformation at n (last) state */
    void Set_grad_u_n(TPZFMatrix<STATE> & grad_u_n)
    {
        m_grad_u_n = grad_u_n;
    }
    
    /** @brief Get Gradient of deformation at n (last) state */
    TPZFMatrix<STATE> & grad_u_n()
    {
        return m_grad_u_n;
    }
    
    
    /** @brief Set elastic strain at n (last) state */
    void Set_epsilon_e_n(TPZFMatrix<REAL> & epsilon_e_n)
    {
        m_epsilon_e_n = epsilon_e_n;
    }
    
    /** @brief Get elastic strain at n (last) state */
    TPZFMatrix<REAL> epsilon_e_n()
    {
        return m_epsilon_e_n;
    }
    
    
    /** @brief Set plastic strain at n (last) state */
    void Set_epsilon_p_n(TPZFMatrix<REAL> & epsilon_p_n)
    {
        m_epsilon_p_n = epsilon_p_n;
    }
    
    /** @brief Get plastic strain at n (last) state */
    TPZFMatrix<REAL> epsilon_p_n(){
        return m_epsilon_p_n;
    }
    
    
    /** @brief Set plastic strain at n (last) state */
    void Set_grad_p_n(TPZFMatrix<REAL> & grad_p_n)
    {
        m_grad_p_n = grad_p_n;
    }
    
    /** @brief Get plastic strain at n (last) state */
    TPZFMatrix<REAL> grad_p_n(){
        return m_grad_p_n;
    }
    
    
    /** @brief Set plastic strain at n (last) state */
    void Set_phi_n(TPZFMatrix<REAL> & phi_n)
    {
        m_phi_n = phi_n;
    }
    
    /** @brief Get plastic strain at n (last) state */
    TPZFMatrix<REAL> phi_n(){
        return m_phi_n;
    }
    
    
    /** @brief Set plastic strain at n (last) state */
    void Set_k_n(TPZFMatrix<REAL> & k_n)
    {
        m_k_n = k_n;
    }
    
    /** @brief Get plastic strain at n (last) state */
    TPZFMatrix<REAL> k_n(){
        return m_k_n;
    }
    
    
    /**
     * Total (elastic+plastic) stress
     */
    TPZTensor<REAL> fSigma;
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory :
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /**
     Write method not implemented.
     
     @param buf TPZStream buffer
     @param withclassid obsolete
     */
    void Write(TPZStream &buf, int withclassid) const
    {
//                buf.Write(&m_Pressure_n);
//                buf.Write(&m_Pressure);
        DebugStop();
    }
    
    /**
     Read method not implemented.
     
     @param buf TPZStream buffer
     @param context obsolete
     */
    void Read(TPZStream &buf, void *context)
    {
        //        buf.Read(&m_Pressure_n);
        //        buf.Read(&m_Pressure);
        DebugStop();
    }
    
    void Print(std::ostream &out) const
    {
        //        out << m_Pressure_n;
        //        out << m_Pressure;
        DebugStop();
    }

    
    /**
     Class unique identifier

     @return integer class id
     */
    virtual int ClassId() const;
    
};

inline std::ostream &operator<<(std::ostream &out,const TPZPMRSMemory &mem)
{
    mem.Print(out);
    return out;
}


#endif /* TPZPMRSMemory_h */
