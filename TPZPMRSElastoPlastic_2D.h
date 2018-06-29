//
//  TPZPMRSCoupling.cpp
//  PZ
//
//  Created by Manouchehr on 6/28/18.
//
//

#ifndef PZELASTOPLASTIC2D_H
#define PZELASTOPLASTIC2D_H

#include "TPZMaterial.h"
#include "TPZMatWithMem.h"
#include "TPZElastoPlasticMem.h"
#include "pzporoelastoplasticmem.h"
#include "TPZMatElastoPlastic.h"
#include "TPZMaterial.h"

#include "TPZSimulationData.h"



/**
 * Implements an elastoplastic material and uses the memory feature to store the damage variables
 * This material works only together with the Plasticity Library.
 */

template <class T, class TMEM = TPZElastoPlasticMem>
class  TPZPMRSElastoPlastic_2D : public TPZMatElastoPlastic<T,TMEM> //, TPZMatWithMem<TMEM>
{
public:
	
	//enum SOLUTIONVARS{ENone = -1};
	/**
	 * Default constructor
	 */
	TPZPMRSElastoPlastic_2D();
	
	/** Creates a material object and inserts it in the vector of
	 *  material pointers of the mesh. Upon return vectorindex
	 *  contains the index of the material object within the
	 *  vector
	 */
	TPZPMRSElastoPlastic_2D(int id ,  int PlaneStrainOrPlaneStress);
	
	/** Creates a material object based on the referred object and
	 *  inserts it in the vector of material pointers of the mesh.
	 *  Upon return vectorindex contains the index of the material
	 *  object within the vector
	 */
	TPZPMRSElastoPlastic_2D(const TPZPMRSElastoPlastic_2D<T,TMEM> &mat);
	
	virtual ~TPZPMRSElastoPlastic_2D();
	
	/** returns the name of the material*/
	virtual std::string Name();
	
	/**returns the integrable dimension of the material*/
	virtual int Dimension() const { return 2; }
	
	/** returns the number of state variables associated with the material*/
	virtual int NStateVariables() { return 2; }
	
    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out);

	/** print out the data associated with the material*/
	virtual void Print(std::ostream &out, const int memory);
	
	
	/**returns the solution associated with the var index based on
	 * the finite element approximation*/
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
	
	/**
	 * It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef);
	
	/**
	 * It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
	
    /**
     * This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition
     */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data);
    
	/**
	 * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);

    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);

	
	/** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
	 * @param DeltaStrain [out]
	 * @param data [in]
	 */
	virtual void ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &DeltaStrain);
    
	
	/** Calls the plasticity template aggregate applyStrainComputeDep method
	 *  @param data [in]
	 *  @param DeltaStrain [in]
	 *  @param Stress [out]
	 *  @param Dep [out]
	 */
	virtual void ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,
									TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep);
	
	/** Calls the plasticity template aggregate applyStrainComputeDep method
	 *  @param data [in]
	 *  @param DeltaStrain [in]
	 *  @param Stress [out]
	 */
	virtual void ApplyDeltaStrain(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,
									TPZFMatrix<REAL> & Stress);
	
	
	/**To create another material of the same type*/
	virtual TPZMaterial * NewMaterial();
	
	
	/**
	 * Unique identifier for serialization purposes
	 */
	public:
virtual int ClassId() const;

	
	/**
	 * Save the element data to a stream
	 */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/**
	 * Read the element data from a stream
	 */
	virtual void Read(TPZStream &buf, void *context);
    
	
protected:
    int fPlaneStrain;
    
    // new file
    
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
    
    /** @brief Permeability of the rock */
    REAL m_k;
    
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

    
};

template <class T, class TMEM>
int TPZPMRSElastoPlastic_2D<T, TMEM>::ClassId() const{
    return Hash("TPZMatElastoPlastic2D") ^ TPZMatElastoPlastic<T,TMEM>::ClassId() << 1;
}

#endif
