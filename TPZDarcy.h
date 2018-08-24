//
//  TPZDarcyMem.cpp
//  PZ
//
//  Created by Manouchehr on Agust 24, 2018.
//
//

#ifndef TPZDarcy_h
#define TPZDarcy_h

#include <stdio.h>
#include "TPZMaterial.h"
#include "TPZMatWithMem.h"
#include "pzbndcond.h"
#include "pzvec.h"
#include "TPZSimulationData.h"
#include <iostream>


class TPZDarcy: public TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >  {
    
protected:
    
    /** @brief define the simulation data */
    TPZSimulationData * m_SimulationData;
    
    /** @brief Problem dimension */
    int m_Dim;
    
    /** @brief Initial Permeability of the rock */
    REAL m_k_0;
    
    /** @brief Fluid viscosity */
    REAL m_eta;
    

    
public:
    
    TPZDarcy();
    
    TPZDarcy(int matid, int dim);
    
    ~TPZDarcy();
    
    /** @brief Copy constructor $ */
    TPZDarcy(const TPZDarcy& other);
    
    /** @brief Copy assignemnt operator $ */
    TPZDarcy & operator = (const TPZDarcy& other);
    
    
    void Print(std::ostream & out);
    
    std::string Name() { return "TPZDarcy"; }
    
    /** @brief compute permeability (Kappa) */
    virtual void Compute_Kappa(TPZMaterialData &data, REAL &kappa);
    
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    /** @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation. */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    
    // ********* Set and Get some functions ********************************************
    
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
    
    /** @brief dimension of the model: */
    void SetDimension(int dimension)
    {
        m_Dim = dimension;
    }
    
    int Dimension() const {return m_Dim;}
    
    
    /** @brief Parameters of rock and fluid: */
    void SetParameters(REAL perm, REAL m_porosity, REAL eta)
    {
        m_k_0 = perm;
        m_eta = eta;
    }
    
};


#endif /* TPZDarcy_hpp */
