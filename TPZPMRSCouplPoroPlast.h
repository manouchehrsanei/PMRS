//
//  TPZPMRSCouplPoroPlast.hpp
//  PZ
//
//  Created by Manouchehr on Jun 27, 2018.
//
//

#ifndef TPZPMRSCouplPoroPlast_H
#define TPZPMRSCouplPoroPlast_H

#include <iostream>

#include "TPZMaterial.h"
#include "pzdiscgal.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "TPZMatElastoPlastic2D.h"
#include "TPZMatWithMem.h"

#include "TPZSimulationData.h"

#include "TPZPMRSMemory.h"




template <class T, class TPZPMRSMemory>
class TPZPMRSCouplPoroPlast : public TPZMatElastoPlastic2D<T,TPZPMRSMemory>
{
    
protected:
    
    /** @brief define the simulation data */
    TPZSimulationData * m_SimulationData;
    
    /** @brief Problem dimension */
    int m_Dim;
    
    /** @brief body force */
    TPZManVector<REAL,3>  m_b;
    
    /** @brief Poison coeficient */
    REAL m_nu;
    REAL m_nuu;
    
    /** @brief first Lame Parameter */
    REAL m_lambda;
    REAL m_lambdau;

    /** @brief Bulk modulus */
    REAL m_K;
    REAL m_Ku;
    
    /** @brief Second Lame Parameter */
    REAL m_mu;
    
    /** @brief constants Biot poroelasticity */
    REAL m_alpha;
    
    /** @brief Storage coefficient poroelasticity */
    REAL m_Se;
    
    /** @brief Intact rock porosity */
    REAL m_porosity_0;
    
    /** @brief Initial Permeability of the rock */
    REAL m_k_0;
    
    /** @brief Fluid viscosity */
    REAL m_eta;
    
    
    /** @brief coehsion of the rock */
    REAL m_c;
    
    /** @brief Friction angle */
    REAL m_phi_f;
    
    
    /** @brief permeability coupling model  */
    int m_k_model;
    
    
    /** @brief Rock density */
    REAL m_rho_s;
    
    /** @brief Fluid density */
    REAL m_rho_f;
    
    
    /** @brief new ************************************************ */
    bool m_SetRunPlasticity;
    
    
    /** @brief State: one ou one+1 */
    enum EState { ELastState = 0, ECurrentState = 1 };
    EState gState;

    
public:
    
    TPZPMRSCouplPoroPlast();
    
    TPZPMRSCouplPoroPlast(int matid, int dim);
    
    ~TPZPMRSCouplPoroPlast();
    
    /** @brief Copy constructor $ */
    TPZPMRSCouplPoroPlast(const TPZPMRSCouplPoroPlast& other);
    
    /** @brief Copy assignemnt operator $ */
    TPZPMRSCouplPoroPlast & operator = (const TPZPMRSCouplPoroPlast& other);
    
    
    void SetLastState(){ gState = ELastState; }
    void SetCurrentState(){ gState = ECurrentState; }
    
    
    int ClassId() const;
    
    void Print(std::ostream & out);
    
    std::string Name() { return "TPZPMRSCouplPoroPlast"; }
    
    int Dimension() const {return m_Dim;}
    
    virtual int NStateVariables();
    
    /** @brief dimension of the model: */
    void SetDimension(int dimension)
    {
        m_Dim = dimension;
    }
    
    
    /** @brief Parameters of rock and fluid: */
    void SetParameters(REAL perm, REAL m_porosity, REAL eta)
    {
        m_k_0 = perm;
        m_eta = eta;
        m_porosity_0 = m_porosity;
    }
    
    /** @brief Set the simulation data */
    void SetSimulationData(TPZSimulationData * SimulationData)
    {
        m_SimulationData = SimulationData;
    }
    
    /** @brief Get the simulation data */
    TPZSimulationData * SimulationData()
    {
        return m_SimulationData;
    }
    
    /** @brief Set the porolastic parameters data */
    void SetPorolasticParameters(REAL l, REAL mu, REAL l_u)
    {
        m_lambda = l;
        m_mu = mu;
        m_lambdau = l_u;
        m_K = m_lambda + (2.0/3.0)*m_mu;
        m_Ku = m_lambdau + (2.0/3.0)*m_mu;
    }
    
    /** @brief Set the porolastic engineer parameters data */
    void SetPorolasticParametersEngineer(REAL Ey, REAL nu)
    {
        
        m_lambda = (Ey*nu)/((1.0+nu)*(1.0-2.0*nu));
        m_mu = (Ey)/(2.0*(1.0+nu));
        m_lambdau = (Ey*nu)/((1.0+nu)*(1.0-2.0*nu));
        m_K = m_lambda + (2.0/3.0)*m_mu;
        m_Ku = m_lambdau + (2.0/3.0)*m_mu;
    }
    
    /** @brief Set the Biot parameters data */
    void SetBiotParameters(REAL alpha, REAL Se)
    {
        if(alpha==0){
            std::cout << "Biot constan should be at leats equal to the intact porosity, alpha = " << alpha  << std::endl;
            DebugStop();
        }
        m_alpha = alpha;
        m_Se = Se;
    }
    
    
    /** @brief set the peremability models: */
    void SetKModel(int model)
    {
        m_k_model = model;
    }
    
    /** @brief return the peremability models: */
    int KModel()
    {
        return m_k_model;
    }
    
    /** @brief permeability correction model */
    REAL k_permeability(REAL &phi, REAL &k);
    
    /** @brief porosity correction model */
    REAL porosity_corrected_2D(TPZVec<TPZMaterialData> &datavec);
    
    REAL porosity_corrected_3D(TPZVec<TPZMaterialData> &datavec);
    
    
    /** @brief computation of effective sigma */
    void Compute_Sigma_n(TPZFMatrix<REAL> Grad_u_n, TPZFMatrix<REAL> Grad_u, TPZFMatrix<REAL> &e_e, TPZFMatrix<REAL> &e_p, TPZFMatrix<REAL> &S);
    
    /** @brief Principal Stress */
    void Principal_Stress(TPZFMatrix<REAL> S, TPZFMatrix<REAL> & PrinStres);
    
    
    
    
    
    
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    
    /** @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation. */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    virtual void ContributePlastic_2D(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
    void Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    void ContributeBC_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    void ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    int VariableIndex(const std::string &name);
    int NSolutionVariables(int var);
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);

    void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right)
    {
        DebugStop();
    }
    
    void ContributeInterface(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                             REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                             REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                               REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        DebugStop();
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
    {
        DebugStop();
    }
    
    
    
};


template <class T, class TPZPMRSMemory>
int TPZPMRSCouplPoroPlast<T,TPZPMRSMemory>::ClassId() const{
    return Hash("TPZPMRSCouplPoroPlast") ^ TPZMatElastoPlastic2D<T,TPZPMRSMemory>::ClassId() << 1;
}


#endif /* TPZPMRSCouplPoroPlast_hpp */
