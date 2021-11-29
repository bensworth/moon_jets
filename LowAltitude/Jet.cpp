// Copyright 2017
//
// The Regents of the University of Colorado, a body corporate
//
// Created by Ben Southworth (Department of Applied Mathematics) and Sascha Kempf
// (Labratory for Atmospheric and Space Physics)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <cfloat>
#include <new>
#include <ctime>
#include <random>
#include <algorithm>
#include "genFunctions.h"
#include "Jet.h"
#include "Solver.h"
#ifdef _OPENMP
    #include <omp.h>
#endif
using namespace std;
using namespace genFunctions;


// Const array size variables; must be defined in cpp file, not header file
// These are hard coded for low altitude simulations.
const int Jet::m_nr = 196;
const int Jet::m_nphi = 45;
const int Jet::m_nvr = 35;
const int Jet::m_nvphi = 20;
const float Jet::m_min_altitude = 0.1;
const float Jet::m_max_altitude = 5;
const float Jet::m_max_velocity = 0.7;
const float Jet::m_max_rphi = 15;
const float Jet::m_max_vphi = 90;

/* Constructor */
Jet::Jet() : m_dataID(0), m_dr((m_max_altitude - m_min_altitude) / m_nr),
    m_dvr(m_max_velocity / m_nvr), m_dphi(m_max_rphi / m_nphi),
    m_dvphi(m_max_vphi / m_nvphi), m_maxAngleRad(DEG2RAD*m_max_rphi)
{

}

/* Deconstructor */
Jet::~Jet()
{

}

/* Function to set the number of variables being integrated. */
void Jet::SetNumVariables(const int & n)
{
    CONST_numVariables = n;
}

/* Function to set the longitude and latitude location of the Jet on Enceladus, and direction */
/* vector of jet's center.                                                                    */
void Jet::SetLocation(const float & latitude, const float & longitude, const float & jetZen,
    const float & jetAz)
{
    m_latitude   = latitude;
    m_longitude  = longitude - 180.; // Must have longitude in [-180,180]
    m_jetZenith  = jetZen;
    m_jetAzimuth = jetAz;
}

/* Function to take in vectors of Enceladus position and velocity, save them as member */
/* then compute the initial particle location and velocity due to Enceladus' velocity and */
/* Enceladus' rotation. Also stores member unit vector initPos with the directional */
/* vector of the particle launch location. */
void Jet::SetInitCond(const vector<double> & pos, const vector<double> & vel)
{
    double lon, lat, vRot, phi, theta;
    m_moonPos = pos;
    m_moonVel = vel;

    // Convert latitude/longitude to radians for this code.
    lon = DEG2RAD*m_longitude;
    lat = DEG2RAD*m_latitude;
    
    // Put longitude and latitude into Cartesian unit direction vector. Note, these
    // basis vectors rotate from planet centered frame to moon centered frame. 
    m_ez = {0., 0., 1.};
    m_ey = Cross(m_ez, m_moonPos);
    Normalize(m_ey);
    m_ex = Cross(m_ey, m_ez);
    Normalize(m_ex);

    // Since transformation is orthogonal, taking transpose we can transfer from
    // moon frame to planet frame - Unit vector of jet's location in inertial frame.
    m_jetPos = {cos(lon)*cos(lat)*m_ex[0] + sin(lon)*cos(lat)*m_ey[0],
                cos(lon)*cos(lat)*m_ex[1] + sin(lon)*cos(lat)*m_ey[1],
                sin(lat)};
    Normalize(m_jetPos);

    // Set initial dust position.
    m_parPos = { m_moonPos[0] + GLOBAL_radiusMoon*m_jetPos[0],
                 m_moonPos[1] + GLOBAL_radiusMoon*m_jetPos[1],
                 m_moonPos[2] + GLOBAL_radiusMoon*m_jetPos[2] };

    // Account for Enceladus rotational speed. 
    vector<double> rotationDir = Cross(m_ez, m_jetPos);
    Normalize(rotationDir);
    vRot = 2.*PI*cos(lat)*GLOBAL_radiusMoon / GLOBAL_periodMoon;

    // Set initial dust velocity.
    m_parNonEjVel = { m_moonVel[0] + vRot*rotationDir[0],
                      m_moonVel[1] + vRot*rotationDir[1],
                      m_moonVel[2] + vRot*rotationDir[2] };

    // Ejection direction of the jet (including tilt)
    m_vz = m_jetPos;
    m_vx = {0., 0., 1.};
    m_vy = Cross(m_vz, m_vx);
        Normalize(m_vy);
    m_vx = Cross(m_vy, m_vz);
        Normalize(m_vx);
    phi   = -m_jetAzimuth*DEG2RAD;        // because it's western longitude
    theta = PI/2. - m_jetZenith*DEG2RAD;  // zenith -> theta
    m_jetDir = { cos(phi)*cos(theta)*m_vx[0] + sin(phi)*cos(theta)*m_vy[0] + sin(theta)*m_vz[0],
                 cos(phi)*cos(theta)*m_vx[1] + sin(phi)*cos(theta)*m_vy[1] + sin(theta)*m_vz[1],
                 cos(phi)*cos(theta)*m_vx[2] + sin(phi)*cos(theta)*m_vy[2] + sin(theta)*m_vz[2] };

    // build orthonormal system for the jet in inertial frame, where m_jetDir 
    // is parallel to +z and +x is in apex direction.
    m_vz = m_jetDir;
    if(Dot(m_ex, m_vz) < 1.) m_vx = m_ex;
    else m_vx = m_ez;
    m_vy = Cross(m_vz, m_vx);
        Normalize(m_vy);
    m_vx = Cross(m_vy, m_vz);
        Normalize(m_vx);
}

/* Function to set the largest angle of ejection, and create a vector of launch directions */
/* uniformly spaced from 0 degrees to the supplied maximum angle, and from 0 to 2pi around. */
void Jet::SetMaxAngle(const float & theta)
{
    // m_maxAngle    = theta;
    // m_maxAngleRad = theta*DEG2RAD;
    std::cout << "Warning : max angle currently hard_coded.\n";
}

/* Function to set data about the simulation including jet Id, code version, plasma model, */
/* jet reference source, and number of cores.                                              */
void Jet::SetSimulationData(const int & jetId, const int & codeVersion, const int & bField, 
    const int & charging, const int & jetRef, const int & numCore)
{
    m_jetId       = jetId;
    m_codeVersion = codeVersion;
    m_bFieldModel = bField;
    m_charging    = charging; 
    m_jetRef      = jetRef;
    m_numCore     = numCore;
}

/* Function to set density grid data for file output purposes. */
void Jet::SetGridData(const int & gridMin_x, const int & gridMax_x,
    const int & gridMin_y, const int & gridMax_y, const int & gridMin_z,
    const int & gridMax_z, const float & gridSize_dx)
{
    // Grid size
    m_gridSize_dx = gridSize_dx;
    // x-dimension grid properties
    m_numGrid_x = (gridMax_x - gridMin_x) / gridSize_dx;
    m_minGrid_x = gridMin_x;
    // y-dimension grid properties
    m_numGrid_y = (gridMax_y - gridMin_y) / gridSize_dx;
    m_minGrid_y = gridMin_y;
    // z-dimension grid properties 
    m_numGrid_z = (gridMax_z - gridMin_z) / gridSize_dx;
    m_minGrid_z = gridMin_z;
}

/* Function to write density data in binary. */
void Jet::DensityWrite(const unordered_map<long int,pair<float,float> > & Density,
    const float & inclination, const int & numAzimuth, const int & partRadId,
    const int & monteCarlo)
{
    stringstream density_out;
    if(GLOBAL_bodyID == 1) {
        if(monteCarlo == 1) {
            density_out << "/scratch/summit/beso3770/Enc_JetResults/Jet" << m_jetId
                << "/D" << m_jetId << "_r" << partRadId << "_a" << inclination << "_BF"
                << m_charging << "-" << m_bFieldModel << "_rc"
                << m_critRad << "_vgas" << m_v_gas << ".dens";  
        }
        else {
            density_out << "/scratch/summit/beso3770/Enc_JetResults/Jet" << m_jetId
                << "/D" << m_jetId << "_r" << partRadId << "_a" << inclination << "_BF"
                << m_charging << "-" << m_bFieldModel << ".dens";   
        }
    }
    else if(GLOBAL_bodyID == 2) {
        if(monteCarlo == 1) {
            density_out << "/scratch/summit/beso3770/Eur_JetResults/Jet" << m_jetId
                << "/D" << m_jetId << "_r" << partRadId << "_a" << inclination << "_BF"
                << m_charging << "-" << m_bFieldModel << "_rc"
                << m_critRad << "_vgas" << m_v_gas << ".dens";  
        }
        else {
            density_out << "/scratch/summit/beso3770/Eur_JetResults/Jet" << m_jetId
                << "/D" << m_jetId << "_r" << partRadId << "_a" << inclination << "_BF"
                << m_charging << "-" << m_bFieldModel << ".dens";   
        }
    }
    ofstream dens_out(density_out.str(), ios::binary);

    // Header Data. 
    int tmp, numDens = Density.size();
    // Version
    Jet::binary_write(dens_out,tmp=1);      
    Jet::binary_write(dens_out,m_codeVersion);
    // Body ID
    Jet::binary_write(dens_out,tmp=18);
    Jet::binary_write(dens_out,GLOBAL_bodyID);
    // Monte Carlo
    Jet::binary_write(dens_out,tmp=19);     
    Jet::binary_write(dens_out,monteCarlo);
    if(monteCarlo == 1) {
        Jet::binary_write(dens_out,tmp=20);     
        Jet::binary_write(dens_out,m_v_gas);
        Jet::binary_write(dens_out,tmp=21);     
        Jet::binary_write(dens_out,m_critRad);
    }
    // Data ID
    Jet::binary_write(dens_out,tmp=0);      
    Jet::binary_write(dens_out,m_dataID);
    // B-field model
    Jet::binary_write(dens_out,tmp=2);      
    Jet::binary_write(dens_out,m_bFieldModel);
    // Charging
    Jet::binary_write(dens_out,tmp=23);     
    Jet::binary_write(dens_out,m_charging);
    // Jet parameter reference
    Jet::binary_write(dens_out,tmp=3);      
    Jet::binary_write(dens_out,m_jetRef);
    // Jet ID
    Jet::binary_write(dens_out,tmp=4);      
    Jet::binary_write(dens_out,m_jetId); 
    // Size ID
    Jet::binary_write(dens_out,tmp=5);      
    Jet::binary_write(dens_out,partRadId);
    // Number of cores
    Jet::binary_write(dens_out,tmp=6);      
    Jet::binary_write(dens_out,m_numCore);
    // Grid data
    Jet::binary_write(dens_out,tmp=8);      
    Jet::binary_write(dens_out,m_numGrid_x);
    Jet::binary_write(dens_out,m_numGrid_y);
    Jet::binary_write(dens_out,m_numGrid_z);
    // Smallest x,y,z values covered by grid
    Jet::binary_write(dens_out,tmp=22);     
    Jet::binary_write(dens_out,m_minGrid_x);
    Jet::binary_write(dens_out,m_minGrid_y);
    Jet::binary_write(dens_out,m_minGrid_z);
    // Dx
    Jet::binary_write(dens_out,tmp=10);     
    Jet::binary_write(dens_out,m_gridSize_dx);
    // Length of data 
    Jet::binary_write(dens_out,tmp=12);     
    Jet::binary_write(dens_out,numDens);
    // Angle inclination
    Jet::binary_write(dens_out,tmp=16);     
    Jet::binary_write(dens_out,inclination);
    // Number of azimuth points for inclination
    Jet::binary_write(dens_out,tmp=17);     
    Jet::binary_write(dens_out,numAzimuth);
    // Density profile indices
    Jet::binary_write(dens_out,tmp=25); 
    for(auto const & data : Density) {
        long int ind = data.first;
        Jet::binary_write(dens_out,ind);    
    }
    // Density at respective indices
    Jet::binary_write(dens_out, tmp=15);    
    for(auto const & data : Density) {
        float val = (data.second).first;
        Jet::binary_write(dens_out,val);
    }   
    // Aggregate particle charge at respective indices
    if(m_charging != 0) {
        Jet::binary_write(dens_out, tmp=24);    
        for(auto const & data : Density) {
            float val = (data.second).second;
            Jet::binary_write(dens_out,val);
        }   
    }
    dens_out.close();
}

/* Function that takes the azimuth and inclination to perpendicular angles of ejection  */
/* in a cone with center direction given by m_jetLong, m_jetLat (degrees). Returns unit */
/* direction vector in the inertial Saturnian frame.                                    */
vector<double> Jet::GetEjecVelocity(const float & azimuth, const float & inc)
{   
    double theta = PI/2. - inc*DEG2RAD; // Because we must have longitude in [-180,180]
    double phi   = -azimuth*DEG2RAD;

    vector<double> dirVec =
      { cos(phi)*cos(theta)*m_vx[0] + sin(phi)*cos(theta)*m_vy[0] + sin(theta)*m_vz[0],
        cos(phi)*cos(theta)*m_vx[1] + sin(phi)*cos(theta)*m_vy[1] + sin(theta)*m_vz[1],
        cos(phi)*cos(theta)*m_vx[2] + sin(phi)*cos(theta)*m_vy[2] + sin(theta)*m_vz[2] };

    return dirVec;
}

/* Function to set normalization for cos^2 angular distrivbution over [inc0,inc1]. Note, */
/* the distribution is over [0,m_maxAngle], and this normalies over a subinterval.       */
void Jet::SetAngDistConst(const float & inc0, const float & inc1)
{
    m_angDistConst = 1. / ( m_maxAngleRad*( sin(PI*inc1*DEG2RAD/m_maxAngleRad) - 
                        sin(PI*inc0*DEG2RAD/m_maxAngleRad) ) + PI*DEG2RAD*(inc1 - inc0) );
}

/* Function to calcuate CDF of cos^2 angular distribution for a given angle ang. Speed */
/* distribution is defined over [0,m_maxAngle], and normalized over [inc0,inc1]. Input */
/* and output are in radians.                                                          */
double Jet::AngleDistribution(const double & inc0, const double & inc1, const double & ang)
{
    return m_angDistConst * ( m_maxAngleRad*( sin(PI*ang/m_maxAngleRad) - 
                sin(PI*inc0/m_maxAngleRad) ) + PI*(ang - inc0) );
}

/* Function to sample a cos^2 angular distribution over [inc0,inc1] by taking the inverse */
/* given some uniform sample in (0,1). Input and output are given in degrees.             */
float Jet::SampleIncAngle(const float & inc0, const float & inc1, const float & sample)
{
    double a0  = inc0*DEG2RAD,
           b0  = inc1*DEG2RAD;
    double a   = a0,
           f_a = -sample,
           b   = b0,
           f_b = 1-sample,
           err = 1,
           tol = 1e-6,
           c, 
           f_c;

    // Simple bisection method to invert angular CDF and sample an angle, given a uniform
    // probability of shift in (0,1) 
    while(err > tol) {
        c = (a+b)/2.;
        f_c = Jet::AngleDistribution(a0,b0,c) - sample;

        if(abs(f_c) < tol) {
            break;
        }
        else if(f_c < 0) {
            a   = c;
            f_a = f_c;
        }
        else if(f_c > 0) {
            b   = c;
            f_b = f_c;
        } 
        err = abs(f_b-f_a);
    }
    // Convert to degrees and return.
    c *= RAD2DEG;
    return c;
}

/* Function that takes in an unordered map of local residence times and particle charges, */
/* to add into a second unordered map of aggregate residence times and charges.           */ 
void Jet::UpdateDensity(unordered_map<long int,pair<float,float> > & local,
    unordered_map<long int,pair<float,float> > & aggregate)
{
    for(auto const & data : local) {
        aggregate[data.first].first  += (data.second).first;
        aggregate[data.first].second += (data.second).second;
    }
}

/* Simulate specific vector of speeds and azimuth angles, for a fixed inclination to the jet */
/* directional vector. Save binary output. Uses Open MP Parallelization.                     */
void Jet::HoverSimOMP(Solver & systemSolver, const int &numAzimuth, const int &partRad_ind,
    const float &partRad, const int &initVel_ind, const float &initVel)
{
    // Initialize variables for computing and storing density profile and collision locations. 
    systemSolver.SetSize(partRad);
    float dphi = 360. / numAzimuth;

    // Construct dist_grid in Solver object
    systemSolver.CreateDistributionGrid(m_min_altitude, m_max_altitude,
        m_nr, m_max_velocity, m_nvr, m_max_rphi, m_nphi, m_max_vphi, m_nvphi);

    // Allocate master data array
    // NOTE : ordering here is different than inner data; inner arrays are
    // ordered for efficiency, this array is ordered for later conveniuence.
    float residenceTime[m_nr][m_nphi][m_nvr][m_nvphi] = {0};

    // Declare OMP parallelism outside of loop so can construct objects
    // once and reuse,
    #pragma omp parallel
    {
    Solver tempSolver = systemSolver;
    // Declare inner OMP variables
    double *y = new double[CONST_numVariables];
    float threadResidenceTime[m_nvphi][m_nphi][m_nvr][m_nr] = {0};

    #pragma omp for schedule(static)
    // OpenMP Parallel Loop over inclination angles
    // - here we use inclination for the OpenMP because since we are
    // storing data in effectively a 2d space-velocity grid, all azimuthal
    // angles for a fixed inclination should traverse similar data bins.
    for (int j=0; j<m_nphi; j++) {

        // Define time step as function of initial velocity and gridsize
        // to ensure particle is counted in most cells it traverses
        double dt = 0.5 * m_dr / initVel;

        // Inclination of ejection direction from orthogonal; take as
        // center of interval, e.g., 0.5 for m_dphi = 1 and j=0.
        double inclination = (j + 0.5) * m_dphi;


        // TODO : ideally, use 2-3 inclination angles / bin for more data.
        // double inc_weight = 
        double inc_weight = 1.0;
        double tempWeight = inc_weight * dt / (numAzimuth);


        // Loop over azimuthal angles
        for (int i=0; i<numAzimuth; i++) {
            
            // Update azimuth and inclination angle and get velocity direction vector.
            float azimuth  = i*dphi;
            vector<double> ejecDir(3); 
            vector<float>  colLoc(7);
            ejecDir = Jet::GetEjecVelocity(azimuth,inclination);

            // Update initial conditions. 
            y[0]  = m_parPos[0];
            y[1]  = m_parPos[1];
            y[2]  = m_parPos[2];
            y[3]  = m_parNonEjVel[0] + initVel*ejecDir[0];
            y[4]  = m_parNonEjVel[1] + initVel*ejecDir[1];
            y[5]  = m_parNonEjVel[2] + initVel*ejecDir[2];
            y[6]  = m_moonPos[0];
            y[7]  = m_moonPos[1];
            y[8]  = m_moonPos[2];
            y[9]  = m_moonVel[0];
            y[10] = m_moonVel[1];
            y[11] = m_moonVel[2];
            if(CONST_numVariables == 13) {
                y[12] = 0.;
            }
            // Simulate particle, add density profile to aggregate jet density, and   
            // collision coordinate to aggregate collision density. 
            tempSolver.HoverSim(dt,y,tempWeight,threadResidenceTime);
        }
    }
    delete [] y;

    // Update master residenceTime with data from inner loops
    #pragma omp critical
    {
        // TODO : check ordering indices correct w/ definition
        for (int i = 0; i < m_nr; i++) {
            for (int j = 0; j < m_nphi; j++) {
                for (int k = 0; k < m_nvr; k++) {
                    for (int l = 0; l < m_nvphi; l++) {
                        residenceTime[i][j][k][l] += threadResidenceTime[l][j][k][i];
                    }
                }
            }
        }
    }

    } // end of OMP

    // DEBUG
    double vtotal_vol = 0;
    double rtotal_vol = 0;
    // DEBUG

    // Construct volume normalization for velocity cone
    // (2pi/3) * (|v_{i+1}|^3 - |v_i|^3) * (cos(vphi_i)-\cos(vphi_{i+1}))
    m_velVolume.resize(m_nvr);
    double s0=0, s1, temp, phi0, phi1;
    double dvphi_rad = DEG2RAD*m_dvphi;
    for (int i=0; i<m_nvr; i++) {
        m_velVolume[i].resize(m_nvphi);
        s1 = s0 + m_dvr;
        temp = s1*s1*s1 - s0*s0*s0;
        s0 = s1;
        phi0 = PI;
        for (int j=0; j<m_nvphi; j++) {
            phi1 = phi0 - dvphi_rad;
            m_velVolume[i][j] = (2*PI/3.0) * temp * (std::cos(phi1) - std::cos(phi0));

            vtotal_vol += m_velVolume[i][j];    // DEBUG: sum volume

            m_velVolume[i][j] *= 1e9;   // Convert to m^3 here for numerical stabillity
            phi0 = phi1;
        }
    }

    // Construct volume normalization for location cone,
    // (2pi/3) * (|p_{i+1}|^3 - |p_i|^3) * (cos(phi_i)-\cos(phi_{i+1}))
    m_locVolume.resize(m_nr);
    double r0=m_min_altitude, r1;
    double dphi_rad = DEG2RAD*m_dphi;
    for (int i=0; i<m_nr; i++) {
        m_locVolume[i].resize(m_nphi);
        r1 = r0 + m_dr;
        temp = r1*r1*r1 - r0*r0*r0;
        r0 = r1;
        phi0 = PI;
        for (int j=0; j<m_nphi; j++) {
            phi1 = phi0 - dphi_rad;
            m_locVolume[i][j] = (2*PI/3.0) * temp * (std::cos(phi1) - std::cos(phi0));
            
            rtotal_vol += m_locVolume[i][j];    // DEBUG: sum volume 

            m_locVolume[i][j] *= 1e9;   // Convert to m^3 here for numerical stabillity
            phi0 = phi1;
        }
    }

    // DEBUG
    // Check that volumes sum to volume of entire cone for vel/location
    temp = m_max_altitude*m_max_altitude*m_max_altitude - m_min_altitude*m_min_altitude*m_min_altitude;
    double vvol = (2*PI/3.0) * temp * (std::cos(m_max_rphi*DEG2RAD) - std::cos(PI));
    temp = m_max_velocity*m_max_velocity*m_max_velocity;
    double rvol = (2*PI/3.0) * temp * (std::cos(m_max_vphi*DEG2RAD) - std::cos(PI));
    if (std::abs(vvol - vtotal_vol) > 1e-5) {
        std::cout << "\n\nWARNING, velocity volumes do not match!! Exact = " << vvol << ", sum = " << vtotal_vol << "\n\n";
    }
    if (std::abs(rvol - rtotal_vol) > 1e-5) {
        std::cout << "\n\nWARNING, spatial volumes do not match!! Exact = " << rvol << ", sum = " << rtotal_vol << "\n\n";
    }
    // DEBUG
    
    // Normalize residenceTime based on 6d-cell volume
    double temp_vol;
    for (int i = 0; i < m_nr; i++) {
        for (int j = 0; j < m_nphi; j++) {
            for (int k = 0; k < m_nvr; k++) {
                for (int l = 0; l < m_nvphi; l++) {
                    temp_vol = m_locVolume[i][j] * m_velVolume[k][l];
                    residenceTime[i][j][k][l] /= temp;
                }
            }
        }
    }

    // Save density profile
    Jet::HDF5DistWrite(residenceTime, numAzimuth, partRad_ind,
        partRad, initVel_ind, initVel);
}






