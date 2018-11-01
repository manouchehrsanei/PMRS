//
//  TPMRSCouplPoroElast.hpp
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#ifndef TPMRSCouplPoroElast_h
#define TPMRSCouplPoroElast_h

#include <stdio.h>
#include "TPZMaterial.h"
#include "TPZMatWithMem.h"
#include "TPMRSCoupPoElaMemory.h"
#include "pzdiscgal.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzvec.h"
#include "TPMRSSimulationData.h"
#include <iostream>
#include <cmath>
#include "pzlog.h"



class TPMRSCouplPoroElast : public TPZMatWithMem<TPMRSCoupPoElaMemory,TPZDiscontinuousGalerkin>
{
    
private:
    
    /// Pointer of Simulation data
    TPMRSSimulationData * m_SimulationData;
    
    /// Brief Problem dimension
    int m_Dim;
    
    /// Brief first Lame Parameter
    REAL m_lambda;
    
    /// Brief Second Lame Parameter
    REAL m_mu;
    
    /// Brief constants Biot poroelasticity
    REAL m_alpha;
    
    /// Brief Fluid viscosity
    REAL m_eta;
    
    /// Brief Storage coefficient poroelasticity
    REAL m_Se;
    
    /// Brief Intact rock porosity
    REAL m_porosity_0;
    
    /// Brief Initial Permeability of the rock
    REAL m_k_0;
    
    /// Brief permeability coupling model
    int m_k_model;
    
    /// Brief body force
    TPZManVector<REAL,3>  m_b;
    
    /// Brief Fluid density
    REAL m_rho_f;
    
    /// Brief Rock density
    REAL m_rho_s;
    
    /// Brief Compressibility of fluid
    REAL m_comp_f;
    
    
public:
    
    /// Default constructor
    TPMRSCouplPoroElast();
    
    /// Constructor based on a material id
    TPMRSCouplPoroElast(int matid, int dim);
    
    /// Destructor
    ~TPMRSCouplPoroElast();
    
    /// Brief Copy constructor
    TPMRSCouplPoroElast(const TPMRSCouplPoroElast& other);
    
    /// Brief Copy assignemnt operator
    TPMRSCouplPoroElast & operator = (const TPMRSCouplPoroElast& other);
    
    /// Print out the data associated with the material
    void Print(std::ostream & out);
    
    std::string Name() { return "TPMRSCouplPoroElast"; }
    
    /// Returns the number of state variables
    int NStateVariables();
    
    /// Brief permeability correction model
    REAL k_permeability(REAL &phi, REAL &k);
    
    /// Brief porosity correction model 2D
    REAL porosity_corrected_2D(TPZVec<TPZMaterialData> &datavec);
    
    /// Brief porosity correction model 3D
    REAL porosity_corrected_3D(TPZVec<TPZMaterialData> &datavec);
    
    /// Brief computation of effective sigma
    void Compute_Sigma_n(TPZFMatrix<REAL> Grad_u_n, TPZFMatrix<REAL> Grad_u, TPZFMatrix<REAL> &e_e, TPZFMatrix<REAL> &e_p, TPZFMatrix<REAL> &S);
    
    /// Brief Principal Stress
    void Principal_Stress(TPZFMatrix<REAL> T, TPZFMatrix<REAL> & S);
    
    /// Set the required data at each integration point
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    /// Brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    void ContributeBC_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    void ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /// Returns the variable index associated with the name.
    int VariableIndex(const std::string &name);
    
    /// Returns the number of variables associated with the variable indexed by var.
    int NSolutionVariables(int var);
    
    /// Returns the solution associated with the var index based on a finite element approximation.
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
    
    
    
    /// ********* Set and Get some functions ********************************************
    
    /// Brief Set the simulation data
    void SetSimulationData(TPMRSSimulationData * SimulationData)
    {
        m_SimulationData = SimulationData;
    }
    
    /// Brief Get the simulation data
    TPMRSSimulationData * SimulationData()
    {
        return m_SimulationData;
    }
    
    /// Brief Set the dimension of the model
    void SetDimension(int dimension)
    {
        m_Dim = dimension;
    }
    
    /// Brief Get the dimension of the model
    int Dimension() const
    {
        return m_Dim;
    }
    
    /// Set First Lamé Parameter
    void Setlambda(REAL Ey, REAL nu)
    {
        m_lambda = (Ey*nu)/((1.0+nu)*(1.0-2.0*nu));
    }
    
    /// Get First Lamé Parameter
    REAL lambda()
    {
        return m_lambda;
    }
    
    /// Set Second Lamé Parameter
    void Setmu(REAL Ey, REAL nu)
    {
        m_mu = (Ey)/(2.0*(1.0+nu));
    }
    
    /// Get Second Lamé Parameter
    REAL mu()
    {
        return m_mu;
    }
    
    /// Set Biot coefficient
    void Setalpha(REAL alphaCof)
    {
        m_alpha = alphaCof;
    }
    
    /// Get Biot coefficient
    REAL alpha()
    {
        return m_alpha;
    }
    
    /// Set Fluid dynamic viscosity
    void Seteta(REAL eta_f)
    {
        m_eta = eta_f;
    }
    
    /// Get Fluid dynamic viscosity
    REAL eta()
    {
        return m_eta;
    }
    
    /// Set Storage coefficient
    void SetSe(REAL SeCof)
    {
        m_Se = SeCof;
    }
    
    /// Get Storage coefficient
    REAL Se()
    {
        return m_Se;
    }
    
    /// Set Initial Porosity
    void Setporosity0(REAL porosityZero)
    {
        m_porosity_0 = porosityZero;
    }
    
    /// Get Initial Porosity
    REAL porosity0()
    {
        return m_porosity_0;
    }
    
    /// Set Initial Permeability
    void Setk0(REAL k_zer0)
    {
        m_k_0 = k_zer0;
    }
    
    /// Get Initial Permeability
    REAL k0()
    {
        return m_k_0;
    }
    
    /// Brief set the peremability models
    void SetKModel(int kmodel)
    {
        m_k_model = kmodel;
    }
    
    /// Brief Get the peremability models
    int KModel()
    {
        return m_k_model;
    }
    
    /// Set Density of fluid
    void Setrhof(REAL rho_f)
    {
        m_rho_f = rho_f;
    }
    
    /// Get Density of fluid
    REAL rhof()
    {
        return m_rho_f;
    }
    
    /// Set Density of fluid
    void Setrhos(REAL rho_s)
    {
        m_rho_s = rho_s;
    }
    
    /// Get Density of rock
    REAL rhos()
    {
        return m_rho_s;
    }
    
    /// Set Compressibility of fluid
    void Setcomf(REAL com_f)
    {
        m_comp_f = com_f;
    }
    
    /// Get Compressibility of fluid
    REAL comf()
    {
        return m_comp_f;
    }
    
    
    /// Set Body force
    void Setbforce(TPZManVector<REAL,3> b_force)
    {
        m_b = b_force;
    }
    
    /// Get Body force
    TPZManVector<REAL,3> bforce()
    {
        return m_b;
    }
    
    
};


#endif /* TPMRSCouplPoroElast_hpp */
