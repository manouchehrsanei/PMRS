//
//  TPZPMRSCoupling.cpp
//  PZ
//
//  Created by Manouchehr on Jun 27, 2018.
//
//

#ifndef TPZMatElastoPlastic2D_H
#define TPZMatElastoPlastic2D_H

#include <iostream>

#include "TPZMaterial.h"
#include "pzdiscgal.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "TPZMatElastoPlastic2D.h"
#include "TPZMatWithMem.h"



#include "TPZSimulationData.h"
#include "TPZPMRSMemory.h"



template <class T, class TMEM = TPZElastoPlasticMem>
class TPZPMRSElastoPlastic_2D : public TPZMatElastoPlastic2D<T,TMEM>
{
	
protected:
	
    
    /** @brief define the simulation data */
    TPZSimulationData * m_SimulationData;
    
    /** @brief Material Id */
    int m_matId;
    
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
    
    /** @Drucker Prager property */
    REAL m_eta_dp;
    REAL m_xi_dp;
    
    
    /** @brief permeability coupling model  */
    int m_k_model;
    
    /** @brief Uses plain stress
     * @note \f$m_PlaneStress = 1\f$ => Plain stress state
     * @note \f$m_PlaneStress != 1\f$ => Plain Strain state
     */
    int m_PlaneStress;
    
    
    /** @brief Rock density */
    REAL m_rho_s;
    
    /** @brief Fluid density */
    REAL m_rho_f;
 
   
    bool m_SetRunPlasticity;
  

	/** @brief State: one ou one+1 */
	enum EState { ELastState = 0, ECurrentState = 1 };
	EState gState;
    
    
    
public:

	
	TPZPMRSElastoPlastic_2D();
	
	TPZPMRSElastoPlastic_2D(int matid, int dim);
	
	virtual ~TPZPMRSElastoPlastic_2D();
    
    
    /** @brief Copy constructor $ */
    TPZPMRSElastoPlastic_2D(const TPZPMRSElastoPlastic_2D& other);
    
    /** @brief Copy assignemnt operator $ */
    TPZPMRSElastoPlastic_2D & operator = (const TPZPMRSElastoPlastic_2D& other);
    
    
    
    /** @brief if IsPlasticity is true, it will calculate the contribution of a plastic material
     */
    void SetRunPlasticity(bool IsPlasticity = true);
    
    int MatId()
    {
        return m_matId;
    }
    
    void SetLastState(){ gState = ELastState; }
    void SetCurrentState(){ gState = ECurrentState; }
    
    
    
        
    int ClassId() const;

	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZPMRSElastoPlastic_2D"; }
	
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
    
    
    /** @Drucker Prager parameters */
    void SetDruckerPragerParameters(REAL phi_f, REAL c)
    {
        m_phi_f = phi_f;
        m_c = c;
        
        // Outer edges condition
        m_eta_dp = 6.0*(sin(m_phi_f))/(sqrt(3.0)*(3.0-sin(m_phi_f)));
        m_xi_dp = 6.0*(cos(m_phi_f))/(sqrt(3.0)*(3.0-sin(m_phi_f)));
        
    }
    
    
    
    /** @brief Set plane problem
     * planestress = 1 => Plain stress state
     * planestress != 1 => Plain Strain state
     */
    void SetPlaneProblem(int planestress)
    {
        m_PlaneStress = planestress;
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
    
    
    REAL porosity_corrected(TPZVec<TPZMaterialData> &datavec);
    REAL porosity_corrected_3D(TPZVec<TPZMaterialData> &datavec);


    REAL k_permeability(REAL &phi, REAL &k);

	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
	 * @param datavec [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
	virtual void ContributePlastic(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);

	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc);
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
	
	
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
	
	virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
	
    
    void Read(TPZStream& buf, void* context)
    {
            DebugStop();
    }

    void Write(TPZStream& buf, int withclassid) const
    {
            DebugStop();
    }

};


template <class T, class TMEM>
int TPZPMRSElastoPlastic_2D<T,TMEM>::ClassId() const{
    return Hash("TPZH1PlasticFrac2D") ^ TPZMatElastoPlastic2D<T,TMEM>::ClassId() << 1;
}

#endif  /* TPZPMRSElastoPlastic_2D_hpp */
