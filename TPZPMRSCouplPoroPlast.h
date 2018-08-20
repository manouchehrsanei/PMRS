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

#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZSandlerDimaggio.h"

#include "TPZSimulationData.h"
#include "TPZElastoPlasticMem.h"
#include "pzporoelastoplasticmem.h"


// TPZElastoPlasticMem
// TPZCouplElasPlastMem


template <class T, class TMEM = TPZPoroElastoPlasticMem>
class TPZPMRSCouplPoroPlast : public TPZMatElastoPlastic2D<T,TMEM>
{
    
private:

    // boolean indicating whether the Young modulus needs to be computed at each point
    int m_VariableYoung;
    
    
protected:
    
    /** @brief define the simulation data */
    TPZSimulationData * m_SimulationData;
    
    /** @brief Problem dimension */
    int m_Dim;
    
    /** @brief body force */
    TPZManVector<REAL,3>  m_b;
    
    /** @brief permeability coupling model  */
    int m_k_model;
    
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
    
    /** @brief Fluid density */
    REAL m_rho_f;
    
    /** @brief Rock density */
    REAL m_rho_s;
    
    /** @brief Cohesion of Mohr-Coloumb */
    REAL mc_coh;
    
    /** @brief Friction of Mohr-Coloumb */
    REAL mc_phi;
    
    /** @brief Dilation of Mohr-Coloumb */
    REAL mc_psi;
    
    
    /** @brief new ************************************************ */
    bool m_SetRunPlasticity;

    /** @brief Flag to indicate if should update should zero displacement and EpsTotal. With this you can the solution vector means U, and not DeltaU */
    bool m_UpdateToUseFullDiplacement;
    
    
    
public:

    
    // Default constructor
    
    TPZPMRSCouplPoroPlast();
    
    TPZPMRSCouplPoroPlast(int matid, int dim);
    
    virtual ~TPZPMRSCouplPoroPlast();
    
    /** @brief Copy constructor $ */
    TPZPMRSCouplPoroPlast(const TPZPMRSCouplPoroPlast& other);
    
    /** @brief Copy assignemnt operator $ */
    TPZPMRSCouplPoroPlast & operator = (const TPZPMRSCouplPoroPlast& other);
    
    virtual int ClassId() const;
    
    void Print(std::ostream & out);
    
    virtual std::string Name() { return "TPZPMRSCouplPoroPlast"; }
    
    virtual int NStateVariables();
    
    /** @brief some functions for plasticity */
    virtual void ComputeStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &Strain);
    virtual void ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &DeltaStrain);
    virtual void ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep);
    virtual void ComputeStressVector(TPZMaterialData & data, TPZFMatrix<REAL> &Stress);
    virtual void ApplyDeltaStrain(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,TPZFMatrix<REAL> & Stress);
    
    
    /** @brief Sets/Unsets the internal memory data to be updated in the next assemble/contribute call */
    void SetUpdateToUseFullU(bool update = true);
    
    /** @brief permeability correction model */
    REAL k_permeability(REAL &phi, REAL &k);
    
    /** @brief porosity correction model */
    REAL porosity_corrected_2D(TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Poroelastic porosity correction from strains and pressure */
    REAL porosity_corrected_2D(TPZTensor<STATE> & eps_elastic, TPZTensor<STATE> & eps_plastic, STATE & pressure);
    
    REAL porosity_corrected_3D(TPZVec<TPZMaterialData> &datavec);
    
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    /** @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation. */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    virtual void Contribute_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    virtual void ContributePlastic_2D(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
    virtual void Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    virtual void ContributePlastic_3D(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ContributeBC_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    virtual int VariableIndex(const std::string &name);
    virtual int NSolutionVariables(int var);
    
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);

    void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right)
    {
        DebugStop();
    }
    
    virtual void ContributeInterface(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                             REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                             REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                               REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        DebugStop();
    }
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
    {
        DebugStop();
    }
    
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
    
    
    /** @brief Parameters of rock and fluid: */
    void SetParameters(REAL perm, REAL m_porosity, REAL eta)
    {
        m_k_0 = perm;
        m_eta = eta;
        m_porosity_0 = m_porosity;
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
    
    
    /** @brief Density of fluid and rock: */
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
    
    
    // ****************************************************************** plasticity ***********************
    
    /** @brief if IsPlasticity is true, it will calculate the contribution of a plastic material
     */
    void SetRunPlasticity(bool IsPlasticity = true);
    

    
    
};


template <class T, class TMEM>
int TPZPMRSCouplPoroPlast<T,TMEM>::ClassId() const{
    return Hash("TPZPMRSCouplPoroPlast") ^ TPZMatElastoPlastic2D<T,TMEM>::ClassId() << 1;
}


#endif /* TPZPMRSCouplPoroPlast_hpp */
