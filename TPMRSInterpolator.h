//
//  TPMRSInterpolator.h
//  PMRS
//
//  Created by Omar and Manouchehr on 2/13/19.
//

#ifndef TPMRSInterpolator_h
#define TPMRSInterpolator_h

#include <stdio.h>
#include "pzreal.h"
#include "pzerror.h"

#include <utility>
#include <vector>


class TPMRSInterpolator {
    
private:
    
    /// Table of time values vs number of functions values being interpolated
    std::vector<std::pair<REAL, std::vector<REAL>>> m_points;

    /// Number of time values
    int m_n_time_data;
    
    /// Number of functions being interpolated
    int m_n_f_data;
    
    /// The requirements for the provided table of values
    bool CheckDataCoherence(std::vector<std::pair<REAL, std::vector<REAL>>> & points);
    
public:
    
    /// Default constructor
    TPMRSInterpolator();
    
    /// Constructor
    TPMRSInterpolator(std::vector<std::pair<REAL, std::vector<REAL>>> & points);
    
    /// Copy constructor
    TPMRSInterpolator(const TPMRSInterpolator & other);
    
    /// Assignement constructor
    const TPMRSInterpolator & operator=(const TPMRSInterpolator & other);
    
    /// Destructor
    ~TPMRSInterpolator();
    
    /// Evaluates linear interpolator
    std::vector<REAL> f(REAL time);
    
    /// Clear object data
    void Clear();
    
    /// Sets the points provided for perform the interpolation
    void SetPoints(std::vector<std::pair<REAL, std::vector<REAL>>> & points);
    
    /// Gets the points provided for perform the interpolation
    std::vector<std::pair<REAL, std::vector<REAL>>> & Points();
    
    /// Return the number of functions being interpolated
    int n_functions();
    
};

#endif /* TPMRSInterpolator_h */
