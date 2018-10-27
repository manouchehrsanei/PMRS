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
#include "TPMRSMemoryPoroElast.h"
#include "pzdiscgal.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzvec.h"
#include "TPMRSSimulationData.h"
#include <iostream>
#include <cmath>
#include "pzlog.h"



class TPMRSCouplPoroElast : public TPZMatWithMem<TPMRSMemoryPoroElast,TPZDiscontinuousGalerkin>
{
    
protected:
    
    /// Brief define the simulation data
    TPMRSSimulationData * m_SimulationData;
    
    /// Brief Problem dimension
    int m_Dim;
    
    /// Brief body force
    TPZManVector<REAL,3>  m_b;
    
    /// Brief permeability coupling model
    int m_k_model;
    
    /// Brief Poison coeficient
    REAL m_nu;
    REAL m_nuu;
    
    /// Brief first Lame Parameter
    REAL m_lambda;
    REAL m_lambdau;

    /// Brief Bulk modulus
    REAL m_K;
    REAL m_Ku;
    
    /// Brief Second Lame Parameter
    REAL m_mu;
    
    /// Brief constants Biot poroelasticity
    REAL m_alpha;
    
    /// Brief Storage coefficient poroelasticity
    REAL m_Se;
    
    /// Brief Intact rock porosity
    REAL m_porosity_0;
    
    /// Brief Initial Permeability of the rock
    REAL m_k_0;
    
    /// Brief Fluid viscosity
    REAL m_eta;
    
    /// Brief Fluid density
    REAL m_rho_f;
    
    /// Brief Rock density
    REAL m_rho_s;
    
    /// Brief Cohesion of Mohr-Coloumb
    REAL mc_coh;
    
    /// Brief Friction of Mohr-Coloumb
    REAL mc_phi;
    
    /// Brief Dilation of Mohr-Coloumb
    REAL mc_psi;
    

    
public:
    
    /// Default constructor
    TPMRSCouplPoroElast();
    
    /// Constructor
    TPMRSCouplPoroElast(int matid, int dim);
    
    /// Destructor
    ~TPMRSCouplPoroElast();
    
    /// Brief Copy constructor
    TPMRSCouplPoroElast(const TPMRSCouplPoroElast& other);
    
    /// Brief Copy assignemnt operator
    TPMRSCouplPoroElast & operator = (const TPMRSCouplPoroElast& other);
    
    
    void Print(std::ostream & out);
    
    std::string Name() { return "TPMRSCouplPoroElast"; }
    
    virtual int NStateVariables();
    
    
    /// Brief permeability correction model
    REAL k_permeability(REAL &phi, REAL &k);
    
    /// Brief porosity correction model
    REAL porosity_corrected_2D(TPZVec<TPZMaterialData> &datavec);
    
    REAL porosity_corrected_3D(TPZVec<TPZMaterialData> &datavec);
    
    /// Brief computation of effective sigma
    void Compute_Sigma_n(TPZFMatrix<REAL> Grad_u_n, TPZFMatrix<REAL> Grad_u, TPZFMatrix<REAL> &e_e, TPZFMatrix<REAL> &e_p, TPZFMatrix<REAL> &S);
    
    /// Brief Principal Stress
    void Principal_Stress(TPZFMatrix<REAL> T, TPZFMatrix<REAL> & S);
    
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    
    /// Brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
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
    
    /// Brief dimension of the model
    void SetDimension(int dimension)
    {
        m_Dim = dimension;
    }
    
    int Dimension() const {return m_Dim;}
    
    
    /// Brief set the peremability models
    void SetKModel(int model)
    {
        m_k_model = model;
    }
    
    /// Brief return the peremability models
    int KModel()
    {
        return m_k_model;
    }
    
    
    /// Brief Parameters of rock and fluid
    void SetParameters(REAL perm, REAL m_porosity, REAL eta)
    {
        m_k_0 = perm;
        m_eta = eta;
        m_porosity_0 = m_porosity;
    }
    
    
    /// Brief Set the porolastic parameters data
    void SetPorolasticParameters(REAL l, REAL mu, REAL l_u)
    {
        m_lambda = l;
        m_mu = mu;
        m_lambdau = l_u;
        m_K = m_lambda + (2.0/3.0)*m_mu;
        m_Ku = m_lambdau + (2.0/3.0)*m_mu;
    }
    
    /// Brief Set the porolastic engineer parameters data
    void SetPorolasticParametersEngineer(REAL Ey, REAL nu)
    {
        
        m_lambda = (Ey*nu)/((1.0+nu)*(1.0-2.0*nu));
        m_mu = (Ey)/(2.0*(1.0+nu));
        m_lambdau = (Ey*nu)/((1.0+nu)*(1.0-2.0*nu));
        m_K = m_lambda + (2.0/3.0)*m_mu;
        m_Ku = m_lambdau + (2.0/3.0)*m_mu;
    }
    
    /// Brief Set the Biot parameters data
    void SetBiotParameters(REAL alpha, REAL Se)
    {
        if(alpha==0){
            std::cout << "Biot constan should be at leats equal to the intact porosity, alpha = " << alpha  << std::endl;
            DebugStop();
        }
        m_alpha = alpha;
        m_Se = Se;
    }
    
    
    /// Brief Density of fluid and rock
    void SetDensityFluidRock(REAL rhof, REAL rhos)
    {
        m_rho_f = rhof;
        m_rho_s = rhos;
    }
    
    void SetMohrCoulombParameters(REAL coh, REAL phi, REAL psi)
    {
        mc_coh = coh;
        mc_phi = phi;
        mc_psi = psi;
    }
    
    
};


#endif /* TPMRSCouplPoroElast_hpp */
