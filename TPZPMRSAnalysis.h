//
//  TPZPMRSAnalysis.hpp
//  PZ
//
//  Created by Omar and Manouchehr on 9/11/18.
//
//

#ifndef TPZPMRSAnalysis_hpp
#define TPZPMRSAnalysis_hpp

#include <stdio.h>

#include <stdio.h>
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPMRSSimulationData.h"
#include "TPZElastoPlasticMem.h"
#include "pzporoelastoplasticmem.h"
#include "pzadmchunk.h"

class TPZPMRSAnalysis : public TPZAnalysis
{
    
private:
    
    /** @brief whether it is PoroElastic */
    bool IsPoroElastic = false;
    
    /** @brief define the simulation data */
    TPMRSSimulationData * m_SimulationData;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 2> m_meshvec;
    
    /** @brief Part of residue at n+1 state  */
    TPZFMatrix<STATE> m_R_n;
    
    /** @brief Part of residue at n (past) state  */
    TPZFMatrix<STATE> m_R;
    
    /** @brief Solution at n+1 state */
    TPZFMatrix<STATE> m_X_n;
    
    /** @brief memory at n+1 state */
    TPZAdmChunkVector<TPZPoroElastoPlasticMem> m_memory_n;
    
    /** @brief Solution at n (past) state */
    TPZFMatrix<STATE> m_X;
    
    /** @brief memory at n (past) state */
    TPZAdmChunkVector<TPZPoroElastoPlasticMem> m_memory;
    
    /** @brief Strain-Stress solution data */
    TPZStack< std::pair<REAL,REAL> > m_strain_stress_duplets;
    
    /** @brief Strain-Porosity solution data */
    TPZStack< std::pair<REAL,REAL> > m_strain_porosity_duplets;
    
    /** @brief Strain-Permeability solution data */
    TPZStack< std::pair<REAL,REAL> > m_strain_permeability_duplets;
    
    /** @brief Strain-Pressure solution data */
    TPZStack< std::pair<REAL,REAL> > m_strain_pressure_duplets;
    
    /** @brief Residue error */
    STATE m_error;
    
    /** @brief Correction variation */
    STATE m_dx_norm;
    
    /** @brief number of newton corrections */
    int m_k_iterations;
    
    
public:
    
    /** @brief default Constructor  */
    TPZPMRSAnalysis();
    
    /** @brief default Destructor  */
    ~TPZPMRSAnalysis();
    
    /** @brief Copy constructor $ */
    TPZPMRSAnalysis(const TPZPMRSAnalysis &copy);
    
    /** @brief Copy assignemnt operator $ */
    TPZPMRSAnalysis &operator=(const TPZPMRSAnalysis &other);
    
    /**
     * @brief    Implements Access methods:
     */
    
    
    /** @brief Set Solution at n+1 (current) state */
    void SetX_n(TPZFMatrix<STATE> &x)
    {
        m_X_n = x;
    }
    
    /** @brief Get Solution at n+1 (current) state */
    TPZFMatrix<STATE> & X_n()
    {
        return m_X_n;
    }
    
    /** @brief Set memory at n+1 state */
    void SetMemory_n(TPZAdmChunkVector<TPZPoroElastoPlasticMem> &memory)
    {
        m_memory_n = memory;
    }
    
    /** @brief Get memory at n+1 state */
    TPZAdmChunkVector<TPZPoroElastoPlasticMem> & GetMemory_n()
    {
        return m_memory_n;
    }
    
    /** @brief Set Solution at n (last) state */
    void SetX(TPZFMatrix<STATE> &x)
    {
        m_X = x;
    }
    
    /** @brief Get Solution at n (last) state */
    TPZFMatrix<STATE> & X()
    {
        return m_X;
    }
    
    /** @brief Set memory at n state */
    void SetMemory(TPZAdmChunkVector<TPZPoroElastoPlasticMem> &memory)
    {
        m_memory = memory;
    }
    
    /** @brief Get memory at n state */
    TPZAdmChunkVector<TPZPoroElastoPlasticMem> & GetMemory()
    {
        return m_memory;
    }
    
    /** @brief Set the simulation data */
    void SetSimulationData(TPMRSSimulationData * SimulationData)
    {
        m_SimulationData = SimulationData;
        m_meshvec.Resize(2);
    }
    
    /** @brief Get the space generator */
    TPMRSSimulationData * SimulationData()
    {
        return m_SimulationData;
    }

    
    /** @brief Set vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure */
    void SetMeshvec(TPZManVector<TPZCompMesh * , 2> &Meshvec)
    {
        m_meshvec = Meshvec;
    }
    /** @brief Get Vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure */
    TPZManVector<TPZCompMesh * , 2> & Meshvec()
    {
        return m_meshvec;
    }
    
    /** @brief Resize and fill residue and solution vectors */
    void AdjustVectors();
    
    /** @brief Set k iterations */
    void Set_k_ietrarions(int k)
    {
        m_k_iterations = k;
    }
    
    /** @brief Get k iterations */
    int k_ietrarions()
    {
        return m_k_iterations;
    }
    
    /** @brief Execute an Euler method step */
    void ExcecuteOneStep();
    
    /** @brief Execute a Quasi Newton iteration  */
    void QuasiNewtonIteration();
    
    /** @brief PostProcessStandard results */
    void PostProcessStepStandard();
    
    /** @brief PostProcess results */
    void PostProcessStep();
    
    /** @brief update last (n) state solution for PMRS_PoroElastic*/
    void Standard_UpdateState();
    
    /** @brief update current (n+1) state solution for PMRS_PoroElastic */
    void Standard_Update_at_n_State();
    
    /** @brief update last (n) state solution */
    void UpdateState();
    
    /** @brief update current (n+1) state solution */
    void Update_at_n_State();

    /** @brief execute the evolutionary problem */
    void Run_Evolution(TPZVec<REAL> & x);
    
    /** @brief Compute the strain and the stress at x euclidean point for each time */
    void AppendStrain_Stress(TPZVec<REAL> & x);
    
    /** @brief Compute the strain and the Porosity at x euclidean point for each time */
    void AppendStrain_Pososity(TPZVec<REAL> & x);
    
    /** @brief Compute the strain and the Permeability at x euclidean point for each time */
    void AppendStrain_Permeability(TPZVec<REAL> & x);
    
    /** @brief Compute the strain and the Pressure at x euclidean point for each time */
    void AppendStrain_Pressure(TPZVec<REAL> & x);
    
    /** @brief Compute the strain and the stress at x euclidean point for each time */
    void PlotStrainStress(std::string file_name);
    
    /** @brief Compute the strain and the Porosity at x euclidean point for each time */
    void PlotStrainPorosity(std::string file_name);
    
    /** @brief Compute the strain and the Permeability at x euclidean point for each time */
    void PlotStrainPermeability(std::string file_name);
    
    /** @brief Compute the strain and the Pressure at x euclidean point for each time */
    void PlotStrainPressure(std::string file_name);
    
};


#endif /* TPZPMRSAnalysis_hpp */
