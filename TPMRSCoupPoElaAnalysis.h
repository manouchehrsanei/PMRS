//
//  TPMRSCoupPoElaAnalysis.hpp
//  PZ
//
//  Created by Omar and Manouchehr on 9/11/18.
//
//

#ifndef TPMRSCoupPoElaAnalysis_hpp
#define TPMRSCoupPoElaAnalysis_hpp

#include <stdio.h>

#include <stdio.h>
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPMRSSimulationData.h"
#include "TPZElastoPlasticMem.h"
#include "pzporoelastoplasticmem.h"
#include "pzadmchunk.h"

class TPMRSCoupPoElaAnalysis : public TPZAnalysis
{
    
private:
    
    /// Brief define the simulation data
    TPMRSSimulationData * m_SimulationData;
    
    /// Brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2
    TPZManVector<TPZCompMesh * , 2> m_meshvec;
    
    /// Brief Part of residue at n+1 state
    TPZFMatrix<STATE> m_R_n;
    
    /// Brief Part of residue at n (past) state
    TPZFMatrix<STATE> m_R;
    
    /// Brief Solution at n+1 state
    TPZFMatrix<STATE> m_X_n;
    
    /// Brief memory at n+1 state
    TPZAdmChunkVector<TPZPoroElastoPlasticMem> m_memory_n;
    
    /// Brief Solution at n (past) state
    TPZFMatrix<STATE> m_X;
    
    /// Brief memory at n (past) state
    TPZAdmChunkVector<TPZPoroElastoPlasticMem> m_memory;
    
    /// Brief Strain-Stress solution data
    TPZStack< std::pair<REAL,REAL> > m_strain_stress_duplets;
    
    /// Brief Strain-Porosity solution data
    TPZStack< std::pair<REAL,REAL> > m_strain_porosity_duplets;
    
    /// Brief Strain-Permeability solution data
    TPZStack< std::pair<REAL,REAL> > m_strain_permeability_duplets;
    
    /// Brief Strain-Pressure solution data
    TPZStack< std::pair<REAL,REAL> > m_strain_pressure_duplets;
    
    /// Brief Residue error
    STATE m_error;
    
    /// Brief Correction variation
    STATE m_dx_norm;
    
    /// Brief number of newton corrections
    int m_k_iterations;
    
    
public:
    
    /// Brief default Constructor
    TPMRSCoupPoElaAnalysis();
    
    /// Brief default Destructor
    ~TPMRSCoupPoElaAnalysis();
    
    /// Brief Copy constructor
    TPMRSCoupPoElaAnalysis(const TPMRSCoupPoElaAnalysis &copy);
    
    /// Brief Copy assignemnt operator
    TPMRSCoupPoElaAnalysis &operator=(const TPMRSCoupPoElaAnalysis &other);
    
    
    /// Brief Set Solution at n+1 (current) state
    void SetX_n(TPZFMatrix<STATE> &x)
    {
        m_X_n = x;
    }
    
    /// Brief Get Solution at n+1 (current) state
    TPZFMatrix<STATE> & X_n()
    {
        return m_X_n;
    }
    
    /// Brief Set memory at n+1 state
    void SetMemory_n(TPZAdmChunkVector<TPZPoroElastoPlasticMem> &memory)
    {
        m_memory_n = memory;
    }
    
    /// Brief Get memory at n+1 state
    TPZAdmChunkVector<TPZPoroElastoPlasticMem> & GetMemory_n()
    {
        return m_memory_n;
    }
    
    /// Brief Set Solution at n (last) state
    void SetX(TPZFMatrix<STATE> &x)
    {
        m_X = x;
    }
    
    /// Brief Get Solution at n (last) state
    TPZFMatrix<STATE> & X()
    {
        return m_X;
    }
    
    /// Brief Set memory at n state
    void SetMemory(TPZAdmChunkVector<TPZPoroElastoPlasticMem> &memory)
    {
        m_memory = memory;
    }
    
    /// Brief Get memory at n state
    TPZAdmChunkVector<TPZPoroElastoPlasticMem> & GetMemory()
    {
        return m_memory;
    }
    
    /// Brief Set the simulation data
    void SetSimulationData(TPMRSSimulationData * SimulationData)
    {
        m_SimulationData = SimulationData;
        m_meshvec.Resize(2);
    }
    
    /// Brief Get the space generator
    TPMRSSimulationData * SimulationData()
    {
        return m_SimulationData;
    }

    /// Brief Set vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure
    void SetMeshvec(TPZManVector<TPZCompMesh * , 2> &Meshvec)
    {
        m_meshvec = Meshvec;
    }
    
    /// Brief Get Vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure
    TPZManVector<TPZCompMesh * , 2> & Meshvec()
    {
        return m_meshvec;
    }
    
    /// Brief Resize and fill residue and solution vectors
    void AdjustVectors();
    
    /// Brief Set k iterations
    void Set_k_ietrarions(int k)
    {
        m_k_iterations = k;
    }
    
    /// Brief Get k iterations
    int k_ietrarions()
    {
        return m_k_iterations;
    }
    
    /// Brief Execute an Euler method step
    void ExcecuteOneStep();
    
    /// Brief Execute a Quasi Newton iteration
    void QuasiNewtonIteration();
    
    /// Brief PostProcessStandard results
    void PostProcessStep();
    
    /// Brief update last (n) state solution for PMRS_PoroElastic
    void UpdateState();
    
    /// Brief update current (n+1) state solution for PMRS_PoroElastic
    void Update_at_n_State();

    /// Brief execute the evolutionary problem
    void Run_Evolution(TPZVec<REAL> & x);
    
    /// Brief Compute the strain and the stress at x euclidean point for each time
    void AppendStrain_Stress(TPZVec<REAL> & x);
    
    /// Brief Compute the strain and the Porosity at x euclidean point for each time
    void AppendStrain_Pososity(TPZVec<REAL> & x);
    
    /// Brief Compute the strain and the Permeability at x euclidean point for each time
    void AppendStrain_Permeability(TPZVec<REAL> & x);
    
    /// Brief Compute the strain and the Pressure at x euclidean point for each time
    void AppendStrain_Pressure(TPZVec<REAL> & x);
    
    /// Brief Compute the strain and the stress at x euclidean point for each time
    void PlotStrainStress(std::string file_name);
    
    /// Brief Compute the strain and the Porosity at x euclidean point for each time
    void PlotStrainPorosity(std::string file_name);
    
    /// Brief Compute the strain and the Permeability at x euclidean point for each time
    void PlotStrainPermeability(std::string file_name);
    
    /// Brief Compute the strain and the Pressure at x euclidean point for each time
    void PlotStrainPressure(std::string file_name);
    
};


#endif /* TPMRSCoupPoElaAnalysis_hpp */
