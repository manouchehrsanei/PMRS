//
//  TPMRSElastoPlastic.h
//  PMRS
//
//  Created by Omar Durán on 9/11/18.
//

#ifndef TPMRSElastoPlastic_h
#define TPMRSElastoPlastic_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TPZSimulationData.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TPZTensor.h"


template <class T, class TMEM>
class TPMRSElastoPlastic : public TPZMatWithMem<TMEM>
{
    
private:
    
    /// Pointer of Simulation data
    TPZSimulationData * m_simulation_data;
    
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
    
    std::string Name(){
        return "TPMRSElastoPlastic<T,TMEM>";
    }
    
    int ClassId() const {
        return Hash("TPMRSElastoPlastic") ^ TPZMatWithMem<TMEM>::ClassId() << 1 ^ T().ClassId() << 2;
    }
    
    int Dimension() const { return m_dimension; }
    
    int NStateVariables() { return m_dimension; }
    
    void Print(std::ostream &out, const int memory);
    
    void Print(std::ostream &out);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc);
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef);
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef, TPZBndCond &bc);
    
    void Write(TPZStream &buf, int withclassid) const;
    
    void Read(TPZStream &buf, void *context);
    
    void FillDataRequirements(TPZMaterialData &data);
    
    void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPZSimulationData * simulation_data);
    
    void SetDimension(int dimension){
        m_dimension = dimension;
    }
    
    int GetDimension(){
        return m_dimension;
    }
    
    void SetPlasticIntegrator(T & plastic_integrator){
        m_plastic_integrator = plastic_integrator;
    }
    
    T & GetPlasticIntegrator(){
        return m_plastic_integrator;
    }
    
    void Epsilon(TPZMaterialData &data, TPZTensor<REAL> & epsilon_t);
    
    void Sigma(TPZTensor<REAL> & epsilon_t, TPZTensor<REAL> & sigma, TPZFMatrix<REAL> * Dep = NULL);
    
};

#endif /* TPMRSElastoPlastic_h */
