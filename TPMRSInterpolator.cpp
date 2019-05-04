//
//  TPMRSInterpolator.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 2/13/19.
//

#include "TPMRSInterpolator.h"


TPMRSInterpolator::TPMRSInterpolator(){
    m_points.clear();
    m_n_time_data   = 0;
    m_n_f_data      = 0;
}

TPMRSInterpolator::TPMRSInterpolator(std::vector<std::pair<REAL, std::vector<REAL>>> &  points){
    if (CheckDataCoherence(points)) {
        m_points      = points;
        m_n_time_data = points.size();
        m_n_f_data    = points[0].second.size();
    }else{
        DebugStop();
    }

}


TPMRSInterpolator::TPMRSInterpolator(const TPMRSInterpolator & other){
    m_points        = other.m_points;
    m_n_time_data   = other.m_n_time_data;
    m_n_f_data      = other.m_n_f_data;
}

const TPMRSInterpolator & TPMRSInterpolator::operator=(const TPMRSInterpolator & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_points        = other.m_points;
    m_n_time_data   = other.m_n_time_data;
    m_n_f_data      = other.m_n_f_data;
    return *this;
}

TPMRSInterpolator::~TPMRSInterpolator(){
    
}

std::vector<REAL> TPMRSInterpolator::f(REAL time){
    
    std::vector<REAL> f_values(m_n_f_data);
    if (m_n_time_data == 1) { // Constant values.
        for (int i = 0; i < m_n_f_data; i++) {
            f_values[i] =  m_points[0].second[i];
        }
        return f_values;
    }
    
    /// Add corresponging controls
    /// lower and upper bounds
    int t_interval = -1;
    for (int t = 0; t < m_n_time_data - 1; t++) {
        REAL t_i = m_points[t].first;
        REAL t_e = m_points[t+1].first;
        if ((IsZero(time-t_i) || t_i < time) && (IsZero(time-t_e) || time < t_e)) {
            t_interval = t;
            break;
        }
    }
    if (t_interval==-1) {
        std::string err{"TPMRSInterpolator:: out of time range"};
//        throw std::range_error(err);
        DebugStop();
    }
    
    REAL t_i = m_points[t_interval].first;
    REAL t_e = m_points[t_interval+1].first;
    for (int i = 0; i < m_n_f_data; i++) {
        REAL f_i = m_points[t_interval].second[i];
        REAL f_e = m_points[t_interval+1].second[i];
        REAL m = (f_e-f_i)/(t_e-t_i);
        f_values[i] =  m*(time-t_i)+f_i;
    }
    return f_values;
}

bool TPMRSInterpolator::CheckDataCoherence(std::vector<std::pair<REAL, std::vector<REAL>>> & points){
    int n_data = points.size();
    if (n_data==0) {
        std::string err{"TPMRSInterpolator:: there is no points to work on."};
//        throw std::range_error(err);
    }
    
    for (int t = 0; t < n_data - 1; t++) {
        REAL t_i = points[t].first;
        REAL t_e = points[t+1].first;
        if(IsZero(t_i-t_e)){
            std::string err{"TPMRSInterpolator:: two time values are equal."};
//            throw std::range_error(err);
            return false;
        }
    }
    int n_f_data = points[0].second.size();
    for (int t = 0; t < n_data; t++) {
        if(n_f_data != points[t].second.size()){
            std::string err{"TPMRSInterpolator:: this is not a cartesian array."};
//            throw std::range_error(err);
            return false;
        }
    }
    return true;
}

void TPMRSInterpolator::SetPoints(std::vector<std::pair<REAL, std::vector<REAL>>> & points){
    if (CheckDataCoherence(points)) {
        m_points      = points;
        m_n_time_data = points.size();
        m_n_f_data    = points[0].second.size();
    }else{
        DebugStop();
    }
}

std::vector<std::pair<REAL, std::vector<REAL>>> & TPMRSInterpolator::Points(){
    return m_points;
}

int TPMRSInterpolator::n_functions(){
    return m_n_f_data;
}

void TPMRSInterpolator::Clear(){
    
    for (int t = 0; t < m_n_time_data; t++) {
        m_points[t].second.clear();
    }
    m_points.clear();
    m_n_time_data   = 0;
    m_n_f_data      = 0;
}
