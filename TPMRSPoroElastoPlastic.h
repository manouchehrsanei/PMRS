//
//  TPMRSPoroElastoPlastic.h
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#ifndef TPMRSPoroElastoPlastic_h
#define TPMRSPoroElastoPlastic_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TPMRSSimulationData.h"
#include "TPZBndCondWithMem.h"
#include "pzaxestools.h"
#include "TPZTensor.h"
#include "TPMRSMemory.h"

template <class T, class TMEM>
class TPMRSPoroElastoPlastic : public TPZMatWithMem<TMEM>
{
    
private:
    
    /// Pointer of simulation data object
    TPMRSSimulationData * m_simulation_data;
    
    /// Problem dimension
    int m_dimension;
    
    /// Plastic integrator object composed by a yield function and a elastic predictor
    T m_plastic_integrator;
    
    /// Fluid compressibility
    STATE m_c;
    
    /// Fluid dynamic viscosity
    STATE m_eta;
    
    /// Fluid density at reference state
    STATE m_rho_0;
    
    /// Scale factor to use numerically balance the penalty number
    STATE m_scale_factor;
    
    /// Defines the porosity model
    TPMRSPhiParameters m_phi_model;
    
    /// Defines the permeability model
    TPMRSKappaParameters m_kappa_model;
    
    /// Defines the use of Backward Euler or Crank–Nicolson method (theta_scheme = 1.0 or theta_scheme = 0.5)
    REAL m_theta_scheme;
    
    /// Defines the block index for displacement equations
    int m_u_b = 0;
    
    /// Defines the block index for pressure equations
    int m_p_b = 1;
    
    
public:
    
    /// Default constructor
    TPMRSPoroElastoPlastic();
    
    /// Constructor based on a material id
    TPMRSPoroElastoPlastic(int matid);
    
    /// Destructor
    ~TPMRSPoroElastoPlastic();
    
    /// Copy constructor
    TPMRSPoroElastoPlastic(const TPMRSPoroElastoPlastic& other);
    
    /// Copy assignemnt operator
    TPMRSPoroElastoPlastic & operator = (const TPMRSPoroElastoPlastic& other);
    
    /// Print out the data associated with the material
    void Print(std::ostream & out);
    
    /// Class Name
    std::string Name() { return "TPMRSPoroElastoPlastic"; }
    
    /// Returns the number of state variables
    int NStateVariables();
    
    /// Set the required data at each integration point
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    /// It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    void ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /// Contribute methods required for be able to complie (ˆ\ˆ)
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    
    /// Returns the variable index associated with the name.
    int VariableIndex(const std::string &name);
    
    /// Returns the number of variables associated with the variable indexed by var.
    int NSolutionVariables(int var);
    
    /// Returns the solution associated with the var index based on a finite element approximation.
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    
    
    /// ********* Access methods ********************************************
    
    /// Set the simulation data
    void SetSimulationData(TPMRSSimulationData * SimulationData)
    {
        m_simulation_data = SimulationData;
    }
    
    /// Get the simulation data
    TPMRSSimulationData * SimulationData()
    {
        return m_simulation_data;
    }
    
    /// Set the dimension of the model
    void SetDimension(int dimension)
    {
        m_dimension = dimension;
    }
    
    /// Get the dimension of the model
    int Dimension() const
    {
        return m_dimension;
    }
    
    /// ********* Access methods for Geomechanics properties ********************************************
    
    /// Set the plastic integrator
    void SetPlasticIntegrator(T & plastic_integrator){
        m_plastic_integrator = plastic_integrator;
    }
    
    /// Get the plastic integrator
    T & GetPlasticIntegrator(){
        return m_plastic_integrator;
    }
    
    /// Epsilon function
    void Epsilon(TPZMaterialData &data, TPZTensor<REAL> & epsilon_t);
    
    /// Sigma function
    void Sigma(TPZMaterialData &data, TPZTensor<REAL> & epsilon_t, TPZTensor<REAL> & sigma, TPZFMatrix<REAL> * Dep = NULL);
    
    
    /// ********* Access methods for Fluid properties ********************************************
    
    /// Set Fluid compressibility
    void Setc(STATE c_f)
    {
        m_c = c_f;
    }
    
    /// Get Fluid compressibility
    STATE c()
    {
        return m_c;
    }
    
    /// Set Fluid dynamic viscosity
    void Seteta(STATE eta_f)
    {
        m_eta = eta_f;
    }
    
    /// Get Fluid dynamic viscosity
    STATE eta()
    {
        return m_eta;
    }
    
    /// Set Fluid density at reference state
    void Setrho_0(STATE rho_0_f)
    {
        m_rho_0 = rho_0_f;
    }
    
    /// Get Fluid density at reference state
    STATE rho_0()
    {
        return m_rho_0;
    }
    
    /// Set scale factor
    void SetScaleFactor(STATE s)
    {
        m_scale_factor = s;
    }
    
    /// Get scale factor
    STATE ScaleFactor()
    {
        return m_scale_factor;
    }
    
    /// Set fluid properties
    void SetFluidProperties(STATE rho, STATE eta, STATE c){
        Setrho_0(rho);
        Seteta(eta);
        Setc(c);
    }
    
    /// Set Crank-Nicolson method
    void SetCrank_Nicolson(){
        m_theta_scheme = 0.5;
    }
    
    
    void porosity(long gp_index, REAL &phi_n, REAL &dphi_ndp, REAL &phi);
    
    void SetPorosityParameters(TPMRSPhiParameters phi_model){
        m_phi_model = phi_model;
    }
    
    TPMRSPhiParameters GetPorosityParameters(){
        return m_phi_model;
    }
    
    void permeability(long gp_index, REAL &kappa, REAL &dkappa_ndphi,REAL &phi,REAL &phi_0);
    
    void SetPermeabilityParameters(TPMRSKappaParameters kappa_model){
        m_kappa_model = kappa_model;
    }
    
    TPMRSKappaParameters GetPermeabilityParameters(){
        return m_kappa_model;
    }
    
};


#endif /* TPMRSPoroElastoPlastic_h */
