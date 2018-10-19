//
//  TPMRSElastoPlastic.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#ifndef TPMRSElastoPlastic_h
#define TPMRSElastoPlastic_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TPMRSSimulationData.h"
#include "TPZBndCondWithMem.h"
#include "pzaxestools.h"
#include "TPZTensor.h"
#include "TPMRSElastoPlasticMemory.h"

template <class T, class TMEM>
class TPMRSElastoPlastic : public TPZMatWithMem<TMEM>
{
    
private:
    
    /// Pointer of Simulation data
    TPMRSSimulationData * m_simulation_data;
    
    /// Material dimension
    int m_dimension;
    
    /// Plastic integrator object composed by a yield function and a elastic predictor
    T m_plastic_integrator;
    
public:
    
    /// Default constructor
    TPMRSElastoPlastic();
    
    /// Destructor
    virtual ~TPMRSElastoPlastic();
    
    /// Constructor base on material id and dimension
    TPMRSElastoPlastic(int mate_id);
    
    /// Constructor based on TPMRSElastoPlastic
    TPMRSElastoPlastic(const TPMRSElastoPlastic & other);
    
    /// Class name
    std::string Name(){
        return "TPMRSElastoPlastic<T,TMEM>";
    }
    
    /// Class Id
    int ClassId() const {
        return Hash("TPMRSElastoPlastic") ^ TPZMatWithMem<TMEM>::ClassId() << 1 ^ T().ClassId() << 2;
    }
    
    /// Dimension of model
    int Dimension() const { return m_dimension; }
    
    /// Number of state variables
    int NStateVariables() { return m_dimension; }
    
    /// Print function with considering memory
    void Print(std::ostream &out, const int memory);
    
    /// Print function
    void Print(std::ostream &out);
    
    /// Returns the variable index associated with the name
    int VariableIndex(const std::string &name);
    
    /// Number of solution related to variables
    int NSolutionVariables(int var);
    
    /// Calculate Secondary variables based on displacement and their derivatives
    void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
    
    /// Brief of contribute method
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc);
    
    void ContributeBC_3D(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc);
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef);
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef, TPZBndCond &bc);
    
    /// Write function
    void Write(TPZStream &buf, int withclassid) const;
    
    /// Read function
    void Read(TPZStream &buf, void *context);
    
    /// Fill material data parameter with necessary requirements for the Contribute
    void FillDataRequirements(TPZMaterialData &data);
    
    /// Fill material data parameter with necessary requirements for the Contribute BC
    void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPMRSSimulationData * simulation_data);
    
    /// Set the dimension of model
    void SetDimension(int dimension){
        m_dimension = dimension;
    }
    
    /// Get the dimension of model
    int GetDimension(){
        return m_dimension;
    }
    
    /// Set the plastic integrator
    void SetPlasticIntegrator(T & plastic_integrator){
        m_plastic_integrator = plastic_integrator;
    }
    
    /// Get the plastic integrator
    T & GetPlasticIntegrator(){
        return m_plastic_integrator;
    }
    
    /// Brief of epsilon function
    void Epsilon(TPZMaterialData &data, TPZTensor<REAL> & epsilon_t);
    
    /// Brief of sigma function
    void Sigma(TPZTensor<REAL> & epsilon_t, TPZTensor<REAL> & sigma, TPZFMatrix<REAL> * Dep = NULL);
    
    /// Brief of contribute Biot stress
    void Contribute_Biot_Stress(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef);
    
};

#endif /* TPMRSElastoPlastic_h */
