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


// TODO : do not think I use this anywhere, does not compile on mac
// #include "sys/types.h"
// #include "sys/sysinfo.h"
// struct sysinfo memInfo;


/* Constructor */
Jet::Jet()
{
	m_dataID = 0;
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
/* vector of jet's center. 																	  */
void Jet::SetLocation(const float & latitude, const float & longitude, const float & jetZen,
	const float & jetAz)
{
	m_latitude   = latitude;
	m_longitude  = longitude - 180.; // Must have longitude in [-180,180]
	m_jetZenith  = jetZen;
	m_jetAzimuth = jetAz;
}

/* Function to set the gas velocity (km/s) and critical radius (um) for the speed */
/* distribution. Only used for Monte Carlo sampling. 							  */
void Jet::SetSpeedDist(const float & v_gas, const float & critRad)
{
	m_v_gas   = v_gas;
	m_critRad = critRad;
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
	m_maxAngle    = theta;
	m_maxAngleRad = theta*DEG2RAD;
}

/* Function to set data about the simulation including jet Id, code version, plasma model, */
/* jet reference source, and number of cores.											   */
void Jet::SetSimulationData(const int & jetId, const int & codeVersion, const int & bField, 
	const int & charging, const int & jetRef, const int & numCore)
{
	m_jetId		  = jetId;
	m_codeVersion = codeVersion;
	m_bFieldModel = bField;
	m_charging    = charging; 
	m_jetRef	  = jetRef;
	m_numCore	  = numCore;
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

/* Function to write collision data in binary. */
void Jet::CollisionWrite(const vector<vector<float> > collisions, const float & inclination,
	const int & particlesLaunched, const int & particlesCollided, const int & partRadId,
	const int & monteCarlo)
{
	stringstream collision_out;
	
	if(GLOBAL_bodyID == 1) {
		if(monteCarlo == 1) { 
			collision_out << "/scratch/summit/beso3770/Enc_JetResults/Jet" << m_jetId
				<< "/C" << m_jetId << "_r" << partRadId << "_a" << inclination << "_BF"
				<< m_charging << "-" << m_bFieldModel << "_rc"
				<< m_critRad << "_vgas" << m_v_gas << ".coll";	
		}
		else {
			collision_out << "/scratch/summit/beso3770/Enc_JetResults/Jet" << m_jetId
				<< "/C" << m_jetId << "_r" << partRadId << "_a" << inclination << "_BF"
				<< m_charging << "-" << m_bFieldModel << ".coll";	
		}
	}
	else if(GLOBAL_bodyID == 2) {
		if(monteCarlo == 1) { 
			collision_out << "/scratch/summit/beso3770/Eur_JetResults/Jet" << m_jetId
				<< "/C" << m_jetId << "_r" << partRadId << "_a" << inclination << "_BF"
				<< m_charging << "-" << m_bFieldModel << "_rc"
				<< m_critRad << "_vgas" << m_v_gas << ".coll";	
		}
		else {
			collision_out << "/scratch/summit/beso3770/Eur_JetResults/Jet" << m_jetId
				<< "/C" << m_jetId << "_r" << partRadId << "_a" << inclination << "_BF"
				<< m_charging << "-" << m_bFieldModel << ".coll";	
		}
	}
	ofstream coll_out(collision_out.str(), ios::binary);

	// Header Data. 
	int tmp, rows = collisions.size(), cols = collisions[0].size(); 
	// Version
	Jet::binary_write(coll_out,tmp=1);		
	Jet::binary_write(coll_out,m_codeVersion);
	// Body ID
	Jet::binary_write(coll_out,tmp=18);
	Jet::binary_write(coll_out,GLOBAL_bodyID);
	// Monte Carlo
	Jet::binary_write(coll_out,tmp=19);		
	Jet::binary_write(coll_out,monteCarlo);
	if(monteCarlo == 1) {
		Jet::binary_write(coll_out,tmp=20);		
		Jet::binary_write(coll_out,m_v_gas);
		Jet::binary_write(coll_out,tmp=21);		
		Jet::binary_write(coll_out,m_critRad);
	}
	// Data ID
	Jet::binary_write(coll_out,tmp=0);		
	Jet::binary_write(coll_out,m_dataID);
	// B-field model
	Jet::binary_write(coll_out,tmp=2);		
	Jet::binary_write(coll_out,m_bFieldModel);
	// Charging
	Jet::binary_write(coll_out,tmp=23);		
	Jet::binary_write(coll_out,m_charging);
	// Jet parameter reference
	Jet::binary_write(coll_out,tmp=3);		
	Jet::binary_write(coll_out,m_jetRef);
	// Jet ID
	Jet::binary_write(coll_out,tmp=4);		
	Jet::binary_write(coll_out,m_jetId); 
	// Size ID
	Jet::binary_write(coll_out,tmp=5);		
	Jet::binary_write(coll_out,partRadId);
	// Number of cores
	Jet::binary_write(coll_out,tmp=6);		
	Jet::binary_write(coll_out,m_numCore);
	// Number particles launched
	Jet::binary_write(coll_out,tmp=8);		
	Jet::binary_write(coll_out,particlesLaunched);
	// Number particles collided
	Jet::binary_write(coll_out,tmp=9);		
	Jet::binary_write(coll_out,particlesCollided);
	// Inclination angle of launch
	Jet::binary_write(coll_out,tmp=10);		
	Jet::binary_write(coll_out,inclination);
	// Length of data 
	Jet::binary_write(coll_out,tmp=11);		
	Jet::binary_write(coll_out,rows);
	Jet::binary_write(coll_out,cols);
	// Collision profile indices
	Jet::binary_write(coll_out,tmp=12);	
	for(int i=0; i<rows; i++) {
		for(int j=0; j<cols; j++) {
			float val = collisions[i][j];
			Jet::binary_write(coll_out,val);		
		}
	}
	coll_out.close();	
}

/* Function that takes the azimuth and inclination to perpendicular angles of ejection  */
/* in a cone with center direction given by m_jetLong, m_jetLat (degrees). Returns unit */
/* direction vector in the inertial Saturnian frame.									*/
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

/* Function to calcuate CDF of initial speed distribution for a given velocity v and  */
/* particle radius partRad. 														  */
double Jet::SpeedDistribution(const double & v, const float & partRad)
{
	return 1. - ( (v/m_v_gas)*(partRad/m_critRad) + 1 ) * pow(1-v/m_v_gas, partRad/m_critRad);
}

/* Function that takes in a particle radius and uniform sample in (0,1) and inverts */
/* the particle CDF, returning a sampled inital speed. 							    */
float Jet::SampleInitSpeed(const float & partRad, const double & sample)
{
	double a   = 0,
		   f_a = -sample,
		   b   = m_v_gas,
		   f_b = 1-sample,
		   err = 1,
		   tol = 1e-6,
		   c, 
		   f_c;

	// Simple bisection method to invert speed CDF and sample a speed, given a uniform
	// probability sample in (0,1) 
	while(err > tol) {
		c = (a+b)/2.;
		f_c = Jet::SpeedDistribution(c,partRad) - sample;

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
	return c;
}

/* Function to set normalization for cos^2 angular distrivbution over [inc0,inc1]. Note, */
/* the distribution is over [0,m_maxAngle], and this normalies over a subinterval. 		 */
void Jet::SetAngDistConst(const float & inc0, const float & inc1)
{
	m_angDistConst = 1. / ( m_maxAngleRad*( sin(PI*inc1*DEG2RAD/m_maxAngleRad) - 
						sin(PI*inc0*DEG2RAD/m_maxAngleRad) ) + PI*DEG2RAD*(inc1 - inc0) );
}

/* Function to calcuate CDF of cos^2 angular distribution for a given angle ang. Speed */
/* distribution is defined over [0,m_maxAngle], and normalized over [inc0,inc1]. Input */
/* and output are in radians. 														   */
double Jet::AngleDistribution(const double & inc0, const double & inc1, const double & ang)
{
	return m_angDistConst * ( m_maxAngleRad*( sin(PI*ang/m_maxAngleRad) - 
				sin(PI*inc0/m_maxAngleRad) ) + PI*(ang - inc0) );
}

/* Function to sample a cos^2 angular distribution over [inc0,inc1] by taking the inverse */
/* given some uniform sample in (0,1). Input and output are given in degrees.  			  */
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
/* to add into a second unordered map of aggregate residence times and charges.			  */ 
void Jet::UpdateDensity(unordered_map<long int,pair<float,float> > & local,
	unordered_map<long int,pair<float,float> > & aggregate)
{
	for(auto const & data : local) {
		aggregate[data.first].first  += (data.second).first;
		aggregate[data.first].second += (data.second).second;
	}
}

/* Simulate specific vector of speeds and azimuth angles, for a fixed inclination to the jet */
/* directional vector. Save binary output. Uses Open MP Parallelization. 			         */
void Jet::SpecSimOMP(Solver & systemSolver, const vector<float> & partSpeeds,
	const vector<double> & weights, const vector<float> & angleData, const float & totalTime,
	const float & volume, const float & partRad, const int & partRadId)
{
	// Initialize variables for computing and storing density profile and collision locations. 
	systemSolver.SetSize(partRad);
	unordered_map<long int,pair<float,float> > densCount;
	vector<vector<float> >   collisions;
	float inclination = angleData[0],
		  numAzimuth  = angleData[1],
		  dphi 		  = 360./angleData[1];
	int particlesLaunched = 0;

	// Determine particle speeds with nonzero weights. 
	int numSpeeds = partSpeeds.size(),
		bottom_ind = 0,
		top_ind = numSpeeds;
	for(int i=0; i<numSpeeds; i++) {
		if(weights[i] > 0) {
			bottom_ind = i;
			break;
		}
	}
	for(int i=numSpeeds-1; i>0; i--) {
		if(weights[i] > 0) {
			top_ind = i;
			break;
		}
	}

	// Loop over particle speeds.
	#pragma omp parallel for reduction(+:particlesLaunched) schedule(dynamic)
	for(int i=top_ind; i>=bottom_ind; i--) {

		// Declare temporary variables
		double *y = new double[CONST_numVariables];
		unordered_map<long int,pair<float,float> > threadCount; 
		vector<vector<float> > threadCollision;
		Solver tempSolver = systemSolver;
		double initVel    = partSpeeds[i];

		// Define time step and number of steps as a function of initial velocity
		// and gridsize. Include time step in residence time weight. Let dt <=60. 
		// Normalize weight with number of azimuth samplings and grid volume.
		double dt = 0.5 * m_gridSize_dx / initVel;
		if(dt > 60) dt = 60.0;
		int simSteps      = ceil(totalTime/dt);
		double tempWeight = weights[i] * dt / (volume*numAzimuth);

		// Loop over azimuth ejection angle.
		for(int j=0; j<numAzimuth; j++) {
			// Update azimuth and inclination angle and get velocity direction vector.
			bool didCollide = false;
			float azimuth  = j*dphi;
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
			tempSolver.ParticleSim(simSteps,dt,y,tempWeight,threadCount,didCollide,colLoc);
			bool isNAN = false;
			if(didCollide == true) {
				colLoc[6] = weights[i] / numAzimuth;
				for (int zz=0; zz<7; zz++) {
					if (std::isnan(colLoc[zz]) ){
						isNAN = true;
					}
				}
				if (isNAN) {
					std::cout << "Warning - NAN encountered in collisions (ignored).\n";
				}
				else {
					threadCollision.push_back(colLoc);					
				}
			}
			if (!isNAN) {
				particlesLaunched = particlesLaunched + 1;
			}
		}
		// Update master density and collision profiles.
		#pragma omp critical
		{
			Jet::UpdateDensity(threadCount, densCount);
			threadCount.clear();
			collisions.insert(collisions.end(),threadCollision.begin(),threadCollision.end());
			threadCollision.clear();
		}
		delete [] y; 
	}

	// Save density profile
	Jet::DensityWrite(densCount,inclination,numAzimuth,partRadId);

	// Save collision data
	int particlesCollided = collisions.size();
	Jet::CollisionWrite(collisions,inclination,particlesLaunched,particlesCollided, partRadId);

	cout << particlesCollided << "/" << particlesLaunched << " of the particles collided -- size: " <<
		partRad << ", angle: " << inclination << endl;
}

/* Simulate specific vector of speeds and azimuth angles, for a fixed inclination to the jet */
/* directional vector. Save binary output. Uses Open MP Parallelization. 			         */
void Jet::MonteCarlo_Jet(Solver & systemSolver, const int & numSpeeds, const int & numAzimuth,
	const float & inclination, const float & d_inclination, const float & totalTime,
	const float & volume, const float & partRad, const int & partRadId)
{
	// Initialize random number generators for azimuth in [0,360) and inclination
	// in inclination (+-) dinclination/2.
	float inc_a, inc_b;
	if(inclination == 0) {
		inc_a = 0;
		inc_b = d_inclination;
	}
	else {
		inc_a = (inclination - d_inclination);
		inc_b = (inclination + d_inclination);
	}
	Jet::SetAngDistConst(inc_a,inc_b);

	// Random number generators.
	random_device 				randDev;
	mt19937 	  				randEngine(randDev());	// Mersenne Twister pseudo-random generator
	uniform_real_distribution<> getUniform(0.,1.);
	uniform_real_distribution<> getAzimuth(0.,360.);

	// Initialize variables for computing and storing density profile and collision locations. 
	systemSolver.SetSize(partRad);
	unordered_map<long int,pair<float,float> > densCount;
	vector<vector<float> >   collisions;
	int particlesLaunched = 0, 
		monteCarlo        = 1;

	// Loop over particle speeds.
	#pragma omp parallel for reduction(+:particlesLaunched) schedule(dynamic)
	for(int i=0; i<numSpeeds; i++) {
		// Declare temporary variables
		double *y = new double[CONST_numVariables];
		unordered_map<long int,pair<float,float> > threadCount; 
		vector<vector<float> > threadCollision;
		Solver tempSolver = systemSolver;

		// Loop over azimuth ejection angle.
		for(int j=0; j<numAzimuth; j++) {
			// Sample initial velocity of particle. 
			double shift   = getUniform(randEngine);
			double initVel = Jet::SampleInitSpeed(partRad, shift);

			// Sample azimiuth angle uniformly over [0,360) and inclination from
			// a cos^2 angular distribution over [inc_a,inc_b).
			float az  = getAzimuth(randEngine),
				  inc = getUniform(randEngine);	  
			inc = Jet::SampleIncAngle(inc_a,inc_b,inc);

			// Define time step and number of steps as a function of initial velocity
			// and gridsize. Include time step in residence time weight. Let dt <=60. 
			// Normalize weight with number of azimuth and speed samplings and grid volume.
			double dt = 0.5 * m_gridSize_dx / initVel;
			if(dt > 60) dt = 60.;
			double tempWeight = dt / (numSpeeds*numAzimuth*volume);
			int simSteps      = ceil(totalTime/dt);

			// Get velocity direction vector, declare variables in OMP thread.
			bool collision = false; 
			vector<double> ejecDir(3); 
			vector<float>  colLoc(7);
			ejecDir = Jet::GetEjecVelocity(az,inc);

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
			tempSolver.ParticleSim(simSteps,dt,y,tempWeight,threadCount,collision,colLoc);

			// Normalize impact weights to represent one particle, add to total count
			bool isNAN = false;
			if(collision == true) {
				colLoc[6] = 1. / (numSpeeds * numAzimuth);
				for (int zz=0; zz<7; zz++) {
					if (std::isnan(colLoc[zz]) ){
						isNAN = true;
					}
				}
				if (isNAN) {
					std::cout << "Warning - NAN encountered in collisions (ignored).\n";
				}
				else {
					threadCollision.push_back(colLoc);					
				}
			}
			if (!isNAN) {
				particlesLaunched = particlesLaunched + 1;
			}
		}

		// Update master density and collision profiles.
		#pragma omp critical
		{
			Jet::UpdateDensity(threadCount, densCount);
			collisions.insert(collisions.end(),threadCollision.begin(),threadCollision.end());
		}
		delete [] y; 
	}

	// Save density profile
	Jet::DensityWrite(densCount,inclination,numAzimuth,partRadId,monteCarlo);

	// Save collision data
	int particlesCollided = collisions.size();
	Jet::CollisionWrite(collisions,inclination,particlesLaunched,
			particlesCollided,partRadId,monteCarlo);

	cout << "Size: " << densCount.size() << endl;
	cout << particlesCollided << "/" << particlesLaunched << " of the particles collided -- size: " <<
		partRad << ", angle: " << inclination << endl;

}


/* Simulate specific vector of speeds and azimuth angles, for a fixed inclination to the jet  */
/* directional vector. Only checks for collsion and saves locaiton in binary output. DOES NOT */
/* RECORD DENSITY. Uses Open MP Parallelization. 			      						      */
void Jet::CollisionMap(Solver & systemSolver, const vector<float> & partSpeeds,
	const vector<double> & weights, const vector<float> & angleData, const float & totalTime,
	const float & volume, const float & partRad, const int & partRadId)
{
	// Initialize variables for computing and storing density profile and collision locations. 
	systemSolver.SetSize(partRad);
	vector<vector<float> >   collisions;
	float inclination = angleData[0],
		  numAzimuth  = angleData[1],
		  dphi 		  = 360./angleData[1];
	int particlesLaunched = 0;

	// Determine particle speeds with nonzero weights. 
	int numSpeeds = partSpeeds.size(),
		bottom_ind = 0,
		top_ind = numSpeeds;
	for(int i=0; i<numSpeeds; i++) {
		if(weights[i] > 0) {
			bottom_ind = i;
			break;
		}
	}
	for(int i=numSpeeds-1; i>0; i--) {
		if(weights[i] > 0) {
			top_ind = i;
			break;
		}
	}

	// Loop over particle speeds.
	#pragma omp parallel for reduction(+:particlesLaunched) schedule(dynamic)
	for(int i=top_ind; i>=bottom_ind; i--) {

		// Declare temporary variables
		double *y = new double[CONST_numVariables];
		vector<vector<float> > threadCollision;
		Solver tempSolver = systemSolver;
		double initVel    = partSpeeds[i];

		// Define time step and number of steps as a function of initial velocity
		// and gridsize. Include time step in residence time weight. Let dt <=60. 
		// Normalize weight with number of azimuth samplings and grid volume.
		double dt = 0.5 * m_gridSize_dx / initVel;
		if(dt > 60) dt = 60.0;
		int simSteps      = ceil(totalTime/dt);
		double tempWeight = weights[i] * dt / (volume*numAzimuth);

		// Loop over azimuth ejection angle.
		for(int j=0; j<numAzimuth; j++) {
			// Update azimuth and inclination angle and get velocity direction vector.
			bool didCollide = false;
			float azimuth  = j*dphi;
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
			tempSolver.CheckCollision(simSteps,dt,y,tempWeight,didCollide,colLoc);
			if(didCollide == true) {
				colLoc[6] = weights[i] / numAzimuth;
				threadCollision.push_back(colLoc);
			}
			particlesLaunched = particlesLaunched + 1;

		}
		// Update master density and collision profiles.
		#pragma omp critical
		{
			collisions.insert(collisions.end(),threadCollision.begin(),threadCollision.end());
			threadCollision.clear();
		}
		delete [] y; 
	}

	// Save (empty) density profile for IDL purposes
	unordered_map<long int,pair<float,float> > densCount;
	densCount[0].first = 0.0;
	densCount[0].second = 0.0;
	Jet::DensityWrite(densCount,inclination,numAzimuth,partRadId);

	// Save collision data
	int particlesCollided = collisions.size();
	Jet::CollisionWrite(collisions,inclination,particlesLaunched,particlesCollided, partRadId);

	cout << particlesCollided << "/" << particlesLaunched << " of the particles collided -- size: " <<
		partRad << ", angle: " << inclination << endl;
}







