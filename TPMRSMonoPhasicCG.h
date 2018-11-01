//
//  TPMRSMonoPhasicCG.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/22/18.
//

#ifndef TPMRSMonoPhasicCG_h
#define TPMRSMonoPhasicCG_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TPMRSSimulationData.h"
#include "TPMRSMonoPhasicMemory.h"
#include "TPMRSPhiParameters.h"
#include "TPMRSKappaParameters.h"

template <class TMEM>
class TPMRSMonoPhasicCG : public TPZMatWithMem<TMEM> {
    
    /// Pointer of Simulation data
    TPMRSSimulationData * m_simulation_data;
    
    /// Dimension
    int m_dimension;
    
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
    
public:
    
    /// Default constructor
    TPMRSMonoPhasicCG();
    
    /// Constructor based on a material id
    TPMRSMonoPhasicCG(int mat_id, int dimension);
    
    /// Constructor based on a TPMRSMonoPhasicCG object
    TPMRSMonoPhasicCG(const TPMRSMonoPhasicCG & other);
    
    /// Constructor based on a TPMRSMonoPhasicCG object
    TPMRSMonoPhasicCG &operator=(const TPMRSMonoPhasicCG & other);
    
    /// Default destructor
    ~TPMRSMonoPhasicCG();
    
    /// Set the required data at each integration point
    void FillDataRequirements(TPZMaterialData &data);
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data);
    
    /// Returns the name of the material
    std::string Name() {
        return "TPMRSMonoPhasicCG";
    }
    
    /// Returns the integrable dimension of the material
    int Dimension() const {return m_dimension;}
    
    /// Returns the number of state variables associated with the material
    int NStateVariables() {return 1;}
    
    virtual TPZMaterial *NewMaterial()
    {
        return new TPMRSMonoPhasicCG(*this);
    }
    
    /// Print out the data associated with the material
    void Print(std::ostream &out = std::cout);
    
    /// Returns the variable index associated with the name
    int VariableIndex(const std::string &name);
    
    /// Returns the number of variables associated with the variable indexed by var.
    int NSolutionVariables(int var);
    
    /// Returns the solution associated with the var index based on a finite element approximation (Used for TPZPostProcAnalysis)
    void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
    
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPMRSSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }
    
    /// Returns the solution associated with the var index based on a finite element approximation
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
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
    
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    // Contribute Methods being used
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void UndrainedContribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
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


#endif /* TPMRSMonoPhasicCGCG_h */
