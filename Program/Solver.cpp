#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <new>
#include "genFunctions.h"
#include "Solver.h"
#include <gsl/gsl_integration.h>
using namespace std;
using namespace genFunctions;

/* Default Constructor */
Solver::Solver()
{
	m_change = 0;
	// Initialize plasma. 
	m_plasma.partRad  = 0.;
    m_plasma.Tec      = 0.;
    m_plasma.Teh      = 0.;
    m_plasma.Tp       = 0.;
    m_plasma.Tw       = 0.;
    m_plasma.Tnu      = 2.5;
    m_plasma.nec      = 0.;
    m_plasma.neh      = 0.;
    m_plasma.np       = 0.;
    m_plasma.nw       = 0.;
    m_plasma.vf       = 0.;
    m_plasma.j0_ec    = 0.;
    m_plasma.j0_eh    = 0.;
    m_plasma.j0_p     = 0.;
    m_plasma.j0_w     = 0.;
    m_plasma.j0_nu    = 0.;
    m_plasma.Xec      = 0.;
    m_plasma.Xeh      = 0.;
    m_plasma.Xp       = 0.;
    m_plasma.Xw       = 0.;
    m_plasma.Xnu      = 0.;
    m_plasma.Mp       = 0.;
    m_plasma.Mw       = 0.;
    m_plasma.Mec      = 0.;
    m_plasma.Meh      = 0.;
    m_plasma.Q_TO_PHI = 0.;
    m_plasma.d_m      = 2.0;
    m_plasma.E_m      = 500.;
    m_plasma.kTs      = 2.0;
    m_plasma.kappa    = 0.1;
    // Set Saturn Bfield coefficients - Connerney (1993)
	CONST_g10 = .21535;
	CONST_g20 = .01642;
	CONST_g30 = .02743;
	// Saturn Bfield coonstants - Simon (2013) - see Enceladus_Bfield_poynting.nb.
	CONST_mu0  = 4.*PI*1e-7,
	CONST_u0   = 26400.,
	CONST_B0   = 325.*1e-9,				// Assume constant strength background field, CONST_B0
	CONST_m0   = 17.6*1.672621*1e-27,
	CONST_n0   = 70.*1e6,
	CONST_a1nN = 0.609798,
	CONST_b1nN = 0.,
	CONST_g1eN = -0.390202,
	CONST_h1eN = 0.,
	CONST_a1iN = 1.3902,
	CONST_b1iN = 0.,
	CONST_g1iN = -0.780404,
	CONST_h1iN = 0.,
	CONST_a1sN = 0.609798,
	CONST_b1sN = 0.;
	// Set Saturn integration constant
	CONST_planIntegrate = 1.5*GLOBAL_J2*GLOBAL_radiusPlan*GLOBAL_radiusPlan;
}

/* Deconstructor */
Solver::~Solver()
{

}

/* Set solver to not include charging in particle equations of motion. */
void Solver::SetNoCharging()
{
	m_charging 		   = 0;
	CONST_numVariables = 12;
	Evaluate 		   = &Solver::EvaluateDeriv_NoCharge;
}

/* Set solver to use fixed charge to mass ratio in particle equations of motion. Input */
/* integer bfield decides whether to use Connery Bfield (bfield=0), or Simon Bfield    */
/* (bfield=1). Default Bfield is Connerney.   									       */
void Solver::SetConstCharge()
{
	m_charging 		   = 1;
	CONST_numVariables = 12;
	Evaluate 		   = &Solver::EvaluateDeriv_ConstCharge;
}

/* Set solver to include electric field, Bfield and plasma in particle equatons of motion. */
/* Input integer bfield decides whether to use Connery Bfield (bfield=0), or Simon Bfield  */
/* (bfield=1). Default Bfield is Connerney. 											   */
void Solver::SetCharging()
{
	m_charging 		   = 2;
	CONST_numVariables = 13;
	Evaluate 		   = &Solver::EvaluateDeriv;
}

/* Set planetary magnetic field to Connery (1993) model - default, input 1 - or Simon (2013) */
/* model - input 2. Note, 0 implies no plasma background.									 */
void Solver::SetBfield(const int & bfield)
{
	if(bfield == 1) {
		m_bfield = 1;
		Bfield   = &Solver::Bfield_Connerney;
	}
	else if(bfield == 2) {
		m_bfield = 2;
		Bfield   = &Solver::Bfield_Simon;
	}
}

/* Function to set the maximum number of extrapolations per iteration and error tolerance. */
/* Also allocates space in vectors for integration. 									   */
void Solver::SetIntegrator(const int & extrapolations, const double & error)
{
	CONST_maxDiv = extrapolations;
	CONST_errTol = error;
	// Set size of vectors to use in modified-Midpoint method.
	m_yTemp1.resize(CONST_numVariables);
	m_yTemp2.resize(CONST_numVariables);
	m_dTemp.resize(CONST_numVariables);

	// Vector of step sizes for extrapolation of modified-Midpoint.
	CONST_extrap.resize(CONST_maxDiv);
	for(int i=0; i<CONST_maxDiv; i++) {
		CONST_extrap[i] = 2.*(i+1.);
	}
	// Vector of vector of vectors to store tableau of data for extrapolation.
	m_extrapTableau.resize(CONST_maxDiv);
	for(int i=0; i<CONST_maxDiv; i++) {
		m_extrapTableau[i].resize(i+1);
		for(int j=0; j<=i; j++) {
			m_extrapTableau[i][j].resize(CONST_numVariables);
		}
	}
}

/* Function to set pole axis of Saturn. */
void Solver::SetPole(const double & x, const double & y, const double & z)
{
	CONST_poleX	= x;
	CONST_poleY	= y;
	CONST_poleZ	= z;
}

/* Function to set initial particle charge. */
void Solver::SetCharge(const float & c)
{
	CONST_potential = c;
}

/* Function to set initial particle radius and mass, given radius and density. */
/* Calls SetCharge2Mass function to set the appropriate charge to mass ratio.  */
void Solver::SetSize(const float & r, const float & rho)
{
	// If we have a large enough particle, it will charge fast enough that we 
	// may assume the charging to happen instantaneously. In this case it attains
	// equilibirum and we use constant charging equations. 
	if(r > 3) {
		CONST_partRad  	   = r;
		CONST_numVariables = 12;
		if(m_charging != 0) {
			Evaluate = &Solver::EvaluateDeriv_ConstCharge;			
		}
	}
	if(CONST_numVariables == 12) {
		CONST_partRad  = r;
		CONST_partMass = 4./3.*PI* pow(CONST_partRad,3) * rho;
		Solver::SetCharge2Mass();
	}
	else {
		CONST_partRad  = r*1e-6;
		CONST_partMass = 4./3.*PI * pow(CONST_partRad,3) * rho * 1000.;
	}
}

/* Function to set charge to mass ratio, after the charge and particle size have been set. */
void Solver::SetCharge2Mass()
{	
	CONST_Qd = .026562*CONST_potential/( CONST_partRad*CONST_partRad );
}

/* Function to set class vectors that were used to change from latitude/longitude moon */
/* system to inertial planetary frame. 												   */
void Solver::SetChangeBasisVec(const vector<double> & ex, const vector<double> & ey,
	const vector<double> & ez)
{
	m_ex = ex;
	m_ey = ey;
	m_ez = ez;
}

/* Function to create a spatial grid centered about moon with cube side length */
/* 'gridSize' and total grid size (+/- sizeX, +/- sizeY, [minZ,maxZ]). 		   */
void Solver::CreateDensityGrid(const int & minX, const int & maxX, 
	const int & minY, const int & maxY, const int & minZ, const int & maxZ,
	const float & gridSize)
{
	CONST_min_x    = minX;
	CONST_max_x    = maxX;
	CONST_min_y    = minY;
	CONST_max_y    = maxY;
	CONST_min_z    = minZ;
	CONST_max_z    = maxZ;
	CONST_gridSize = gridSize;
	m_Nx = (maxX - minX)/CONST_gridSize;
	m_Ny = (maxY - minY)/CONST_gridSize;
}

/* Function to set the number grid points we will use to discretize the surface of */
/* the moon. Note, we let dx = dy = dz, i.e. grid areas are square. Sets size of   */
/* grid based on size of the moon and number of desired grid points. 	    	   */
void Solver::CreateSurfaceGrid(const int & numGrids)
{
	CONST_surfN = numGrids;
	CONST_surfD = 2.*GLOBAL_radiusMoon/numGrids;
}

bool Solver::Abs_compare(const double & a, const double & b)
{
    return (abs(a) < abs(b));
}

/* Function to save simulation data. */
void Solver::SaveData(vector<vector<double> > & data, const int & rows, const int & columns)
{
	ofstream output("../Results/Data.csv");
	for(int p=0; p<rows; p++) {
		for(int q=0; q<columns; q++) {
			if(q == columns-1) output << data[p][q] << endl;
			else output << data[p][q] << ", ";
		}
	}
	output.close();	
}

/* Function to transfer particle coordinates (Px,Py,Pz) centered at planet to coordinates */
/* centered at moon, with (+y)-axis pointing the direction of rotation and (-x)-axis 	  */
/* pointing towards Saturn.																  */
void Solver::Transform(float & Px, float & Py, float & Pz, const double & Mx,
	const double & My, const double & Mz)
{
	// Change of basis vectors
	vector<double> ex,ey,ez, moonPos;
	ez 		= {0,0,1};
	moonPos = {Mx,My,Mz};
	ey 		= Cross(ez,moonPos);
	Normalize(ey);
	ex 		= Cross(ey,ez);
	Normalize(ex);

	// Transform particle
	double tempX = Px, 
		   tempY = Py, 
		   tempZ = Pz;
	Px = tempX * ex[0] + tempY * ex[1];
	Py = tempX * ey[0] + tempY * ey[1];
	Pz = tempZ;
}

/* Function to take in a vector of vectors of particle collision locations on 	   */
/* Enceladus in cartesian coordinates and return a vector of vectors containing    */
/* the associated geographical coordinates (latitude, longitude) in degrees. Note, */
/* this is using the assumption that Enceladus is spherical. 					   */
void Solver::Cartesian2Geog(vector<double> & collisionLocation)
{
	vector<double> loc(2);
	loc[0] = 90. - RAD2DEG * acos( collisionLocation[2]/GLOBAL_radiusMoon );
	loc[1] = RAD2DEG * atan2( collisionLocation[1], collisionLocation[0] );

	if(loc[1] < 0) loc[1] += 360.;
	collisionLocation = loc;
}

/* Function to transform physical coordinates (x,y,z) into spatial indices (i,j,k) for a 	   */
/* gridded tensor, and then transfer gridded indices (i,j,k) into an index for a 1-dimensional */ 
/* array. NOTE: we will initially  map (x,y,z) to indices (i,j,k), indexed s.t. i=0 is the 	   */
/* smallest value of x, j=0 the smallest value of y, etc. Tensor indexing (i,j,k) is done in   */
/* standard left to right, top to bottom, front to back. 									   */
long int Solver::GetDensityIndex(const double & x, const double & y, const double & z)
{	
	long int x_ind = floor( (x-CONST_min_x)/CONST_gridSize ),
		     y_ind = floor( (y-CONST_min_y)/CONST_gridSize ),
		     z_ind = floor( (z-CONST_min_z)/CONST_gridSize );

	long int ind   = (m_Nx*m_Ny*z_ind + m_Nx*y_ind + x_ind);
	return ind;	
}

/* Function to update the surface collision density profile with a new location. Takes input:   */
/* new impact location in Euclidean coordinates to include, and the current surface density     */
/* profile as an unordered map, where the density profile is a cell index along with the number */
/* of collisions in that cell. Function converts Euclidean coordinates to cell index, and       */
/* updates unordered map of density profile. 					     						    */
int Solver::GetCollisionIndex(vector<double> & newColl)
{
	// Project onto cube, translate face point lies in to a positive plane.  
	auto biggest = max_element(newColl.begin(), newColl.end(), Abs_compare);
	double maxEl = *biggest; 
	int maxInd 	 = distance(newColl.begin(), biggest);
	newColl.erase(biggest);
	for(int i=0; i<2; i++) {
		newColl[i] = GLOBAL_radiusMoon * (1. + newColl[i]/maxEl); 
	}

	// Determine which face of cube point lies on (see notes).  
	int faceNum, localInd;
	if(maxEl < 0) {
		if(maxInd == 2) faceNum = 0; 
		else if(maxInd == 1) faceNum = 2; 
		else faceNum = 4;
	}
	else {
		if(maxInd == 2) faceNum = 1;
		else if(maxInd == 1) faceNum = 3;
		else faceNum = 5;
	}
	// Get final index. 
	localInd = (CONST_surfN - 1)*floor(newColl[0] / CONST_surfD) + floor(newColl[1] / CONST_surfD);
	return faceNum*CONST_surfN*CONST_surfN + localInd;
}

/* Returns the dipole field L value at the observer location x. */
double Solver::X2L(const double & x, const double & y, const double & z)
{
	double r = sqrt( x*x + y*y + z*z ) / GLOBAL_radiusPlan;
	return r / (1. - pow(z*r/GLOBAL_radiusPlan, 2) );  
}

/* Function to compute an objects's inclination w.r.t. the plane perpendicular to */
/* the pole of Saturn. */ 
double Solver::Declination(const double & x, const double & y, const double & z,
	const double & vx, const double & vy, const double & vz)
{		
	// Get angular moment vector L.
	vector<double> r, v, L;
	r = {x,y,z};
		Normalize(r);
	v = {vx,vy,vz};
		Normalize(v);
	L = Cross(r,v);
		Normalize(L);

	// Compute declination.
	return acos( L[0]*CONST_poleX + L[1]*CONST_poleY + L[2]*CONST_poleZ );
}

/* Function to compute the euclidean B-field components */ 
vector<double> Solver::Bfield_Connerney(const double & x, const double & y, const double & z, 
	const vector<double> & moonPos, const double & r)
{
	vector<double> B(3);
	double cost    = z/r;
	double rRatio  = GLOBAL_radiusPlan/r;
	double rRatio3 = rRatio*rRatio*rRatio;
	double rRatio4 = rRatio3*rRatio;
	double rRatio5 = rRatio4*rRatio;

	// Transform Connerey (1993) equations to Cartesian coordinates. See Beckmann. Note, B-field 
	// is computed in Gauss not Tesla ( 1G = 10^(-4)T ) and we want in Tesla. 
	if(abs(cost) < 1) {	
		double B0 = 3.*rRatio3*CONST_g10*cost + 1.5*rRatio4*CONST_g20*( 5.*cost*cost - 1. ) + 
						.5*rRatio5*CONST_g30*cost*( 23.*cost*cost - 15. );
		double Bz = rRatio3*CONST_g10 + 3.*rRatio4*CONST_g20*cost + 1.5*rRatio5*CONST_g30*
			( cost*cost - 1. );
		B[0] = 1e-4*B0*x/r;
		B[1] = 1e-4*B0*y/r;
		B[2] = 1e-4*(B0*z/r - Bz);
	}
	else {
		double Br = 2.*rRatio3*CONST_g10*cost + 1.5*rRatio4*CONST_g20*( 3.*cost*cost - 1.) + 
						2.*rRatio5*CONST_g30*cost*( 5.*cost*cost - 3. );
		B[0] = 0.;
		B[1] = 0.;
		B[2] = 1e-4*Br*z/r;
	}
	return B;
}

/* Function to compute the euclidean B-field components from Simon (2013). */
vector<double> Solver::Bfield_Simon(const double & x, const double & y, const double & z, 
	const vector<double> & moonPos, const double & r)
{
	// Declare variables
	vector<double> B(3);
	double E0  = CONST_u0*CONST_B0,
		   va  = CONST_B0 / sqrt(CONST_m0*CONST_n0*CONST_mu0);
	double Ma  = CONST_u0/va;
	double sa  = 1. / (CONST_mu0*va*sqrt(1+pow(Ma,2)) ),
		   fak = 1. / sqrt(1 + Ma*Ma);
	double s0  = fak / (CONST_mu0*va);

	// Transform coordinates to be centered at Enceladus, with (+x)-axis pointing the direction 
	// of rotation and (+y)-axis pointing towards Saturn.calculuate particle distance.
	// Change of basis vectors
	vector<double> ex,ey,ez;
	ez 		= {0,0,1};
	ex 		= Cross(ez,moonPos);
	Normalize(ex);
	ey 		= Cross(ez,ex);
	Normalize(ex);

	// Transform particle, declare other variables
	double Px = x*ex[0] + y*ex[1],
		   Py = x*ey[0] + y*ey[1],
		   Pz = z,
		   dxPhi,
		   dyPhi, 	
		   partR;
	partR = sqrt(Px*Px + Py*Py);

	// Check if particle is in northern wing (z >0).
	if(Pz > 0) {
		// Check if particle is outside Enceladus tube.
		if(partR > GLOBAL_radiusMoon) {
			// Partial derivatives of potential for particle outside of Enceladus tube.
			dxPhi = E0*pow(GLOBAL_radiusMoon,2)*( CONST_g1eN*(-2.*Px*Py) / pow(partR,4) +
						CONST_h1eN*(pow(Py,2) - pow(Px,2)) / pow(partR,4) );
			dyPhi = E0*pow(GLOBAL_radiusMoon,2)*( CONST_g1eN*(pow(Px,2) - pow(Py,2)) / pow(partR,4) +
						CONST_h1eN*(-2.*Px*Py) / pow(partR,4) ) + E0;

			// B-field components for particle outside Enceladus tube and in northern wing. 
			B[0] = fak * ( CONST_mu0*sa*dyPhi - Ma*sqrt( pow(CONST_B0,2) - pow(CONST_mu0,2)*pow(sa,2) * 
						( pow(dxPhi,2)/(pow(Ma,2)+1) + pow(dyPhi,2) ) ) );
			B[1] = -fak * CONST_mu0 * sa * dxPhi;
			B[2] = -fak * ( CONST_mu0*sa*Ma*dyPhi + sqrt( pow(CONST_B0,2) - pow(CONST_mu0,2)*pow(sa,2) * 
						( pow(dxPhi,2)/(pow(Ma,2)+1) + pow(dyPhi,2) ) ) );
		}
		// If not, it is inside. 
		else {
			// Partial derivatives of potential for particle inside northern wing of Enceladus tube.
			dxPhi = E0*CONST_b1nN;
			dyPhi = E0*CONST_a1nN;

			// B-field components for particle inside northern wing of Enceladus tube.
			B[0] = fak * ( CONST_mu0*sa*dyPhi - Ma*sqrt( pow(CONST_B0,2) - pow(CONST_mu0,2)*pow(sa,2) * 
						( pow(dxPhi,2)/(pow(Ma,2)+1) + pow(dyPhi,2) ) ) );
			B[1] = -fak * CONST_mu0 * sa * dxPhi;
			B[2] = -fak * ( CONST_mu0*sa*Ma*dyPhi + sqrt( pow(CONST_B0,2) - pow(CONST_mu0,2)*pow(sa,2) * 
						( pow(dxPhi,2)/(pow(Ma,2)+1) + pow(dyPhi,2) ) ) );
		}
	}
	// Otherwise particle is in southern wing.
	else {
		// Check if particle is outside Enceladus tube.
		if(partR > GLOBAL_radiusMoon) {
			// Partial derivatives of potential if particle is outside of Enceladus tube.
			dxPhi = E0*pow(GLOBAL_radiusMoon,2)*( CONST_g1eN*(-2.*Px*Py) / pow(partR,4) +
						CONST_h1eN*(pow(Py,2) - pow(Px,2)) / pow(partR,4) );
			dyPhi = E0*pow(GLOBAL_radiusMoon,2)*( CONST_g1eN*(pow(Px,2) - pow(Py,2)) / pow(partR,4) +
						CONST_h1eN*(-2.*Px*Py) / pow(partR,4) ) + E0;

			// B-field components for a particle outside the Enceladus tube in the southern wing. 
			B[0] = -fak * ( CONST_mu0*sa*dyPhi - Ma*sqrt( pow(CONST_B0,2) - pow(CONST_mu0,2)*pow(sa,2) * 
						( pow(dxPhi,2)/(pow(Ma,2)+1) + pow(dyPhi,2) ) ) );
			B[1] = fak * CONST_mu0 * sa * dxPhi;
			B[2] = -fak * ( CONST_mu0*sa*Ma*dyPhi + sqrt( pow(CONST_B0,2) - pow(CONST_mu0,2)*pow(sa,2) * 
						( pow(dxPhi,2)/(pow(Ma,2)+1) + pow(dyPhi,2) ) ) );
		}
		// Check if particle is outside plume atmosphere.
		else if(partR > GLOBAL_radiusPlume) {
			// Partial derivatives of potential for particle inside southern wing of Enceladus tube,
			// but outside plume atmosphere.
			dxPhi = E0*CONST_b1iN + E0*pow(GLOBAL_radiusMoon,2)*( CONST_g1iN*(-2.*Px*Py) / pow(partR,4) +
						CONST_h1iN*(pow(Py,2) - pow(Px,2)) / pow(partR,4) );
			dyPhi = E0*CONST_a1iN + E0*pow(GLOBAL_radiusMoon,2)*( CONST_g1iN*(pow(Px,2) - pow(Py,2)) / 
						pow(partR,4) + CONST_h1iN*(-2.*Px*Py) / pow(partR,4) );

			// B-field components for particle inside southern wing of Enceladus tube, but outside 
			// plume atmosphere.
			B[0] = -fak * ( CONST_mu0*sa*dyPhi - Ma*sqrt( pow(CONST_B0,2) - pow(CONST_mu0,2)*pow(sa,2) * 
						( pow(dxPhi,2)/(pow(Ma,2)+1) + pow(dyPhi,2) ) ) );
			B[1] = fak * CONST_mu0 * sa * dxPhi;
			B[2] = -fak * ( CONST_mu0*sa*Ma*dyPhi + sqrt( pow(CONST_B0,2) - pow(CONST_mu0,2)*pow(sa,2) * 
						( pow(dxPhi,2)/(pow(Ma,2)+1) + pow(dyPhi,2) ) ) );
		}
		// If neither of the above, particle is inside plume atmosphere.
		else {
			// Partial derivatives of potential for particle inside southern wing of Enceladus tube
			// and within plume atmosphere.
			dxPhi = E0*CONST_b1sN;
			dyPhi = E0*CONST_a1sN;

			// B-field components for particle inside southern wing of Enceladus tube and within 
			// plume atmosphere.
			B[0] = -fak * ( CONST_mu0*sa*dyPhi - Ma*sqrt( pow(CONST_B0,2) - pow(CONST_mu0,2)*pow(sa,2) * 
						( pow(dxPhi,2)/(pow(Ma,2)+1) + pow(dyPhi,2) ) ) );
			B[1] = fak * CONST_mu0 * sa * dxPhi;
			B[2] = -fak * ( CONST_mu0*sa*Ma*dyPhi + sqrt( pow(CONST_B0,2) - pow(CONST_mu0,2)*pow(sa,2) * 
						( pow(dxPhi,2)/(pow(Ma,2)+1) + pow(dyPhi,2) ) ) );
		}
	}
	// Transform Bfield components from Enceladus system to Saturn system. Use Px,Py as temp variaxbles. 
	Px = B[0]*ex[0] + B[1]*ey[0];
	Py = B[0]*ex[1] + B[1]*ey[1];
	B[0] = Px;
	B[1] = Py;

	return B;
}

/* Get particle potential phi from charge Q. */ 
double Solver::PhiFromQ(const double & partQ)
{
	double epsilon0 = 8.85418781762e-12;
	return partQ / ( 4.*PI*epsilon0*CONST_partRad );
}

/* Solve special function via recursion relation given by ________ ? */
double Solver::SpecialF(const double & x, const int & n)
{
	if(n == 0) {
		return ( x*x * sqrt(PI/x) * exp(.25/x) * (1.-erf(0.5/sqrt(x))) ) / 2.;
	}
	else if(n == 1) {
		return ( x*x - Solver::SpecialF(x,0) ) / (2.*x);
	}
	else {
		double f;
		double fac = 2.*x;
		double f1  = Solver::SpecialF(x,1);
		double f2  = ( Solver::SpecialF(x,0) - f1 ) / fac;
		for(int i=2; i<=n-1; i++) {
			f  = (i*f1 - f2) / fac;
			f1 = f2;
			f2 = f;
		}
		return f;
	}
}

/* Function defined in Horyani (1996) to be integrated. */
double FB(double u, void * params)
{
	double x = *(double *) params;
	return pow(u,5) * exp( -u*(u*x + 1) );
}

/* Solve for special function integrated from B to oo by taking the difference of the */ 
/* function integrated from 0 to oo and 0 to B.   									  */
double Solver::SpecialFB(double x, double B, const int & n)
{
	// Generate structure for integration with room for maxInt intervals.
	int maxInt = 1000;
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(maxInt);
	// Create gsl compatible function, referencing 'FB' function defined in this code. 
	gsl_function F;
    F.function = &(FB);
    F.params = &x;	
    // Integer key based on kind of integration that can take values {1,2,...,6}, which
    // correspond to Gauss-Kronrod integration with 15,21,31,...,61 points. 
    int intKey = 5;
    // Initialize variables for results and error. Call integrator.
	double diff, err;	
	gsl_integration_qag(&F, 0., B, CONST_errTol, CONST_errTol, maxInt, intKey, workspace, &diff, &err);
	gsl_integration_workspace_free(workspace);
	
	// Return integral from B to oo. 
	return Solver::SpecialF(x,n) - x*x*diff; 
}

/* Initialize plasma variables in data structure, given input of Enceladus' euclidean */
/* coordinates. 																	  */
void Solver::SetPlasma(const double & x, const double & y, const double & z)
{
	double L = Solver::X2L(x,y,z);
	if(L < 2.5) L = 2.5;
    m_plasma.partRad = CONST_partRad;

    // Saturnian magnetosphere
	if(L < 20.)  {
		// Temperature cold e- [eV] 
		m_plasma.Tec = L < 9. ? .7*pow(6./0.7, (L-2.5)/6.5) : 6.;
		// Temperature hot e- [eV] 
		m_plasma.Teh = L <= 8. ? 800. :
			( L < 10. ? 800.*pow(0.125, (L-8.)/2.) :
				100. );
		// Temperature protons [eV]         
		m_plasma.Tp  = L < 3. ? 10. :
			(L < 3.5 ? 10.*pow(0.23, (L-3.)/0.5) : 
				(L < 4. ? 2.3 :
					(L < 6.5 ? 2.3*pow(4.34783, (L-4)/2.5) :
						10. )));
		// Temperature water group ions [eV] 
		m_plasma.Tw  = L < 3. ? 100 :
			(L < 3.5 ? 100.*pow(0.4, (L-3.)/0.5) :
				(L < 4. ? 40.*pow(0.75, (L-3.5)/0.5) :
					(L < 6.5 ? 30.*pow(0.3, (L-4.)/2.5) :
						100. )));
		// Number density cold e- [cm^-3]               
		m_plasma.nec = L < 5. ? 100.*pow(0.21, (L-2.5)/2.5) :
			(L < 10. ? 21.*pow(12./21., (L-5.)/5.) : 
				(L < 13. ? 12.*pow(0.7/12., (L-10.)/3.) :
					0.7 ));
		// Number density hot e- [cm^-3]            
		m_plasma.neh = L < 8. ? 0.2*pow(0.1/0.2, (L-2.5)/5.5) :
			0.1*pow(0.01/0.1, (L-8.)/12.);
		// Number density protons [cm^-3]   
		m_plasma.np  = L < 4. ? 0.2*pow(4./0.2, (L-2.5)/1.5) :
			(L < 4.5 ? 4.*pow(6./4., (L-4.)/0.5) :
				(L < 6.3 ? 6.*pow(3./6., (L-4.5)/1.8) :
					(L < 6.5 ? 3.*pow(1.2/3., (L-6.3)/0.2) :
						1.2*pow(0.7/1.2, (L-8.)/1.5) )));
		// Number density water group ions [cm^-3]                 
		m_plasma.nw  = L < 3.5 ? 0.2*pow(30./0.2, (L-2.5)/1.) :
			(L < 4.5 ? 30.*pow(22./30, (L-3.5)/1.) :
				(L < 6.3 ? 22.*pow(10./22., (L-4.5)/1.8) :
					(L < 6.5 ? 10.*pow(5./10., (L-6.3)/0.2) :
						5.*pow(1./5., (L-8.)/1.5) )));

		double density = (m_plasma.nec + m_plasma.neh) / (m_plasma.np + m_plasma.nw);
		m_plasma.np   *= density;
		m_plasma.nw   *= density;

		// subcorotation at L
		m_plasma.vf = L <= 3. ? 0.9 :
			(L <= 10. ? 0.9 - 0.4/7.*(L-3.) :
				(L <= 11.5 ? 0.5 + 0.1*(L-10.) :
					0.65 ));  
	}
	// Solar wind (not currently accounted for).
	else
	{

	}

	double m_p = 1.672621777e-27,	 // proton mass [kg] 
		   m_e = 9.1093829140e-31,    // electron mass [kg]
		   e   = 1.602e-19;           // elementary charge [C]
	double m_w = 16. * m_p;           // water group ion mass [kg]

	// Electron and ion charged currents are described by the OML 
	// approximation as given in Horanyi, Annu. Rev. Astron. Astrophys.
	// 1996, 34, 383-416 (H96).  

	double vec = sqrt( 8. * e * m_plasma.Tec/ (PI*m_e) ),    // cold electron thermal speed
		   veh = sqrt( 8. * e * m_plasma.Teh/ (PI*m_e) ),    // hot electron thermal speed
		   vp  = sqrt( 8. * e * m_plasma.Tp / (PI*m_p) ),    // proton thermal speed
	       vw  = sqrt( 8. * e * m_plasma.Tw / (PI*m_w) ),	 // water group ions thermal speed 
		   f   = pow(m_plasma.partRad,2) * e * PI;

	// OLM current at phi = 0 eV 
	m_plasma.j0_ec = - m_plasma.nec * 1e+6 * vec * f;
	m_plasma.j0_eh = - m_plasma.neh * 1e+6 * veh * f;
	m_plasma.j0_p  =   m_plasma.np  * 1e+6 * vp  * f;
	m_plasma.j0_w  =   m_plasma.nw  * 1e+6 * vw  * f;
	m_plasma.j0_nu =   2.*2.5e+14*f / (GLOBAL_plan2sun*GLOBAL_plan2sun*1e+6);

	// Precompute chi parameters @ phi = 1 eV
	m_plasma.Xec   = -1./m_plasma.Tec;
	m_plasma.Xeh   = -1./m_plasma.Teh;
	m_plasma.Xp    =  1./m_plasma.Tp;
	m_plasma.Xw    =  1./m_plasma.Tw;
	m_plasma.Xnu   =  1./m_plasma.Tnu;

	// Precompute mach numbers @ v = 1 km/s
	m_plasma.Mp    = 1e+3 / sqrt( 2. * e * m_plasma.Tp/m_p );
	m_plasma.Mw    = 1e+3 / sqrt( 2. * e * m_plasma.Tw/m_w );
	m_plasma.Mec   = 1e+3 / sqrt( 2. * e * m_plasma.Tec/m_e );
	m_plasma.Meh   = 1e+3 / sqrt( 2. * e * m_plasma.Teh/m_e );

	// Precompute special functions for efficiency
	double spFct = m_plasma.E_m/4.;	
	CONST_SpecialF5_Xec = Solver::SpecialFB(-m_plasma.Xec*spFct, spFct, 5);
	CONST_SpecialF5_Xeh = Solver::SpecialFB(-m_plasma.Xeh*spFct, spFct, 5);
	CONST_SpecialF_Xec  = Solver::SpecialF(-m_plasma.Xec*spFct, 5);
	CONST_SpecialF_Xeh  = Solver::SpecialF(-m_plasma.Xeh*spFct, 5);
}

/* Compute change in particle charge given current potential and relative velocity. */
double Solver::phi2j(const double & phi, const double & v)
{
   // mach numbers @ v
   double Mp  = m_plasma.Mp  * v;
   double Mw  = m_plasma.Mw  * v;

	// spherical particle emits 2 e_sec
   double fes   = 3.7*2.*m_plasma.d_m;
   
   // Declare current variables.
   double J_ec, J_eh, J_nu, J_p, J_w, J_sec, J_seh;

   // Current onto positive grain
   if(phi > 0.) {
		// Cold e- collection current (Horanyi 1996 eq. 4)
		J_ec = m_plasma.j0_ec * ( 1-m_plasma.Xec*phi );
		// Hot e- collection current (Horanyi 1996 eq. 4)
		J_eh = m_plasma.j0_eh * ( 1-m_plasma.Xeh*phi ); 
		// Photo current (Horanyi 1996 eq. 13)
		J_nu = m_plasma.j0_nu * m_plasma.kappa * exp( -m_plasma.Xnu*phi );
		// Proton collection j (Horanyi 1996 eq. 15)
		J_p  = m_plasma.j0_p/4. * ( (pow(Mp,2) + 0.5 - m_plasma.Xp*phi) * sqrt(PI)/Mp *
			( erf( Mp + sqrt(m_plasma.Xp*phi) ) + erf( Mp-sqrt(m_plasma.Xp*phi) ) ) + 
			( sqrt(m_plasma.Xp*phi)/Mp+1 ) * exp( -pow(Mp-sqrt(m_plasma.Xp*phi),2) ) -    
			( sqrt(m_plasma.Xp*phi)/Mp-1 ) * exp( -pow(Mp+sqrt(m_plasma.Xp*phi),2) ) );
		// Water group ions collection j (Horanyi 1996 eq. 15)          
		J_w  = m_plasma.j0_w/4. * ( (pow(Mw,2) + 0.5 - m_plasma.Xw*phi)*sqrt(PI)/Mw *  
			( erf( Mw + sqrt(m_plasma.Xw*phi) ) + erf( Mw-sqrt(m_plasma.Xw*phi) ) ) +
			( sqrt(m_plasma.Xw*phi)/Mw+1 ) * exp( -pow(Mw-sqrt(m_plasma.Xw*phi),2) ) -
			( sqrt(m_plasma.Xw*phi)/Mw-1 ) * exp( -pow(Mw+sqrt(m_plasma.Xw*phi),2) ) );
		// Secondary e- due to cold e- (Horanyi 1996 eq. 12)       
		J_sec = fes * ( -m_plasma.j0_ec ) * ( 1.+1./m_plasma.kTs*phi ) * exp( phi*(-1./m_plasma.kTs-m_plasma.Xec) ) * 
			CONST_SpecialF5_Xec;
		// Secondary e- due to hot e- (Horanyi 1996 eq. 12) 
		J_seh = exp( phi*( -1./m_plasma.kTs-m_plasma.Xeh ) ) * CONST_SpecialF5_Xeh;
	}
	// Current onto negative grain
	else {
		// cold e- collection current (Horanyi 1996 eq. 4)      
		J_ec = m_plasma.j0_ec * exp( -m_plasma.Xec*phi );
		// hot e- collection current (Horanyi 1996 eq. 4)
		J_eh = m_plasma.j0_eh * exp( -m_plasma.Xeh*phi ); 
		// photo current (Horanyi 1996 eq. 13)
		J_nu = m_plasma.j0_nu * m_plasma.kappa;  
		// proton collection j (Horanyi 1996 eq. 14)
		J_p  = m_plasma.j0_p/2. * ( sqrt(PI)/Mp * (pow(Mp,2) + 0.5 - m_plasma.Xp*phi) * erf(Mp) + exp(-pow(Mp,2)) );
		// water group ions collection j (Horanyi 1996 eq. 14)          
		J_w  = m_plasma.j0_w/2. * ( sqrt(PI)/Mw * (pow(Mw,2) + 0.5 - m_plasma.Xw*phi) * erf(Mw) + exp(-pow(Mw,2)) );
		// secondary e- due to cold e-  (Horanyi 1996 eq. 10)      
		J_sec = fes*(-J_ec) * CONST_SpecialF_Xec;
		// secondary e- due to hot e- (Horanyi 1996 eq. 10)             
		J_seh = fes*(-J_eh) * CONST_SpecialF_Xeh;
	}
   
	// SAVING CURRENTS TO PLOT
  	// vector<double> currents = {J_ec*1e19,J_eh*1e19,J_nu*1e19,J_p*1e19,J_w*1e19,J_sec*1e19,J_seh*1e19,phi,Mp,Mw};
  	// m_currents.push_back(currents);

	return J_ec + J_eh + J_nu + J_p + J_w + J_sec + J_seh;
}

/* System of equations without charging equations. */
void Solver::EvaluateDeriv_NoCharge(vector<double> & Derivs)
{
	// Take previous values from pointer Derivs.
	vector<double> P  = { Derivs[0], Derivs[1],  Derivs[2] }, 
				   Pv = { Derivs[3], Derivs[4],  Derivs[5] },
				   E  = { Derivs[6], Derivs[7],  Derivs[8] },
				   Ev = { Derivs[9], Derivs[10], Derivs[11] };

	// Particle Constants
	double part2plan   = Norm(P),									 	  // particle distance to planet
		   part2plan_2 = part2plan*part2plan,						      // particle distance squared
		   part2moon   = Distance(P,E),							     	  // particle distance to moon
		   P_moon      = GLOBAL_muMoon / (part2moon*part2moon*part2moon), // Particle/moon equation of motion constant
		   P_plan      = GLOBAL_muPlan / pow(part2plan,5), 		     	  // Particle/planet equation of motion constant
		   sinP_2  	   = (P[2]/part2plan) * (P[2]/part2plan);		      // sin^2(declination) = (z/r)^2
	
	// Moon Constants
	double moon2plan   = Norm(E),							   // Moon distance to planet
		   moon2plan_2 = moon2plan * moon2plan,				   // Moon distance to planet squared
		   M_plan  	   = GLOBAL_muPlan / pow(moon2plan,5),	   // Planet/moon equation of motion constant
		   sinE_2  	   = (E[2]/moon2plan) * (E[2]/moon2plan);  // sin^2(declination) = (z/r)^2

	// Particle Equations of motion.
	Derivs[0] = Pv[0];
	Derivs[1] = Pv[1];	
	Derivs[2] = Pv[2];
	Derivs[3] = -P_plan*( part2plan_2 - CONST_planIntegrate*(5.*sinP_2 - 1.) )*P[0] - P_moon*(P[0]-E[0]);
	Derivs[4] = -P_plan*( part2plan_2 - CONST_planIntegrate*(5.*sinP_2 - 1.) )*P[1] - P_moon*(P[1]-E[1]);
	Derivs[5] = -P_plan*( part2plan_2 - CONST_planIntegrate*(5.*sinP_2 - 3.) )*P[2] - P_moon*(P[2]-E[2]);

	// Enceladus equations of motion.
	Derivs[6]  = Ev[0];
	Derivs[7]  = Ev[1];	
	Derivs[8]  = Ev[2];
	Derivs[9]  = -M_plan*( moon2plan_2 - CONST_planIntegrate*(5.*sinE_2 - 1.) )*E[0];
	Derivs[10] = -M_plan*( moon2plan_2 - CONST_planIntegrate*(5.*sinE_2 - 1.) )*E[1];
	Derivs[11] = -M_plan*( moon2plan_2 - CONST_planIntegrate*(5.*sinE_2 - 3.) )*E[2];
}


/* System of equations with a constant charge. Make sure variable CONST_Qd is set. */
void Solver::EvaluateDeriv_ConstCharge(vector<double> & Derivs)
{
	// Take previous values from pointer Derivs.
	vector<double> P  = { Derivs[0], Derivs[1],  Derivs[2] }, 
				   Pv = { Derivs[3], Derivs[4],  Derivs[5] },
				   E  = { Derivs[6], Derivs[7],  Derivs[8] },
				   Ev = { Derivs[9], Derivs[10], Derivs[11] };

	// Particle Constants
	double part2plan   = Norm(P),									 	  // particle distance to planet
		   part2plan_2 = part2plan*part2plan,						      // particle distance squared
		   part2moon   = Distance(P,E),							     	  // particle distance to moon
		   P_moon      = GLOBAL_muMoon / (part2moon*part2moon*part2moon), // Particle/moon equation of motion constant
		   P_plan      = GLOBAL_muPlan / pow(part2plan,5), 		     	  // Particle/planet equation of motion constant
		   sinP_2  	   = (P[2]/part2plan) * (P[2]/part2plan);		      // sin^2(declination) = (z/r)^2
	
	// Moon Constants
	double moon2plan   = Norm(E),							   // Moon distance to planet
		   moon2plan_2 = moon2plan * moon2plan,				   // Moon distance to planet squared
		   M_plan  	   = GLOBAL_muPlan / pow(moon2plan,5),	   // Planet/moon equation of motion constant
		   sinE_2  	   = (E[2]/moon2plan) * (E[2]/moon2plan);  // sin^2(declination) = (z/r)^2

	// Charge Equation constants. Temp vector refers to (v - Om x r), for particle velocity
	// v, position r and magnetic field angular velocity Om. Om has orbital frequency of GLOBAL_Omega_p
	// hours, and momentum vector parallel to the pole, due to rigid corotation.  
	vector<double> field, temp;
	double Om_x = CONST_poleX * 2.*PI/( 60.*60.*GLOBAL_Omega_p ),
		   Om_y = CONST_poleY * 2.*PI/( 60.*60.*GLOBAL_Omega_p ),
		   Om_z = CONST_poleZ * 2.*PI/( 60.*60.*GLOBAL_Omega_p );
	field = (this->*Bfield)(P[0],P[1],P[2],E,part2plan); 
	temp  = { Pv[0]- Om_y*P[2] + Om_z*P[1],
			  Pv[1]- Om_z*P[0] + Om_x*P[2],
			  Pv[2]- Om_x*P[1] + Om_y*P[0] };
	field = Cross(temp,field);

	// Particle Equations of motion.
	Derivs[0] = Pv[0];
	Derivs[1] = Pv[1];	
	Derivs[2] = Pv[2];
	Derivs[3] = -P_plan*( part2plan_2 - CONST_planIntegrate*(5.*sinP_2 - 1.) )*P[0] -
					P_moon*(P[0]-E[0]) + CONST_Qd*field[0];
	Derivs[4] = -P_plan*( part2plan_2 - CONST_planIntegrate*(5.*sinP_2 - 1.) )*P[1] -
					P_moon*(P[1]-E[1]) + CONST_Qd*field[1];
	Derivs[5] = -P_plan*( part2plan_2 - CONST_planIntegrate*(5.*sinP_2 - 3.) )*P[2] -
					P_moon*(P[2]-E[2]) + CONST_Qd*field[2];

	// Enceladus equations of motion.
	Derivs[6]  = Ev[0];
	Derivs[7]  = Ev[1];	
	Derivs[8]  = Ev[2];
	Derivs[9]  = -M_plan*( moon2plan_2 - CONST_planIntegrate*(5.*sinE_2 - 1.) )*E[0];
	Derivs[10] = -M_plan*( moon2plan_2 - CONST_planIntegrate*(5.*sinE_2 - 1.) )*E[1];
	Derivs[11] = -M_plan*( moon2plan_2 - CONST_planIntegrate*(5.*sinE_2 - 3.) )*E[2];
}

/* System of equations */
void Solver::EvaluateDeriv(vector<double> & Derivs)
{
	// Take previous values from pointer Derivs.
	vector<double> P  = { Derivs[0], Derivs[1],  Derivs[2] }, 
				   Pv = { Derivs[3], Derivs[4],  Derivs[5] },
				   E  = { Derivs[6], Derivs[7],  Derivs[8] },
				   Ev = { Derivs[9], Derivs[10], Derivs[11] };
	double Q = Derivs[12];

	// Particle Constants
	double part2plan   = Norm(P),									 	  // particle distance to planet
		   part2plan_2 = part2plan*part2plan,						      // particle distance squared
		   part2moon   = Distance(P,E),							     	  // particle distance to moon
		   P_moon      = GLOBAL_muMoon / (part2moon*part2moon*part2moon), // Particle/moon equation of motion constant
		   P_plan      = GLOBAL_muPlan / pow(part2plan,5), 		     	  // Particle/planet equation of motion constant
		   sinP_2  	   = (P[2]/part2plan) * (P[2]/part2plan);		      // sin^2(declination) = (z/r)^2
	
	// Moon Constants
	double moon2plan   = Norm(E),							   // Moon distance to planet
		   moon2plan_2 = moon2plan * moon2plan,				   // Moon distance to planet squared
		   M_plan  	   = GLOBAL_muPlan / pow(moon2plan,5),	   // Planet/moon equation of motion constant
		   sinE_2  	   = (E[2]/moon2plan) * (E[2]/moon2plan);  // sin^2(declination) = (z/r)^2

	// Charge Equation constants. Temp vector refers to (v - Om x r), for particle velocity
	// v, position r and magnetic field angular velocity Om. Om has orbital frequency of GLOBAL_Omega_p
	// hours, and momentum vector parallel to the pole, due to rigid corotation.
	vector<double> field, temp;
	double omega = 2.*PI / (60.*60.*GLOBAL_Omega_p),
		   Om_x  = CONST_poleX * omega,
		   Om_y  = CONST_poleY * omega,
		   Om_z  = CONST_poleZ * omega;
	field = (this->*Bfield)(P[0],P[1],P[2],E,part2plan);
	temp  = { Pv[0]- Om_y*P[2] + Om_z*P[1],
			  Pv[1]- Om_z*P[0] + Om_x*P[2],
			  Pv[2]- Om_x*P[1] + Om_y*P[0] };
	field = Cross(temp,field);
	double phi = Solver::PhiFromQ(Q); 
	double vp  = sqrt( (Pv[0] + m_plasma.vf*omega*P[1])*(Pv[0] + m_plasma.vf*omega*P[1])  
		 			+ (Pv[1] - m_plasma.vf*omega*P[0])*(Pv[1] - m_plasma.vf*omega*P[0])  
		 			+ Pv[2]*Pv[2]);

	// Particle Equations of motion.
	Derivs[0] = Pv[0];
	Derivs[1] = Pv[1];	
	Derivs[2] = Pv[2];
	Derivs[3] = -P_plan*( part2plan_2 - CONST_planIntegrate*(5.*sinP_2 - 1.) )*P[0] -
					P_moon*(P[0]-E[0]) + Q/CONST_partMass*field[0];
	Derivs[4] = -P_plan*( part2plan_2 - CONST_planIntegrate*(5.*sinP_2 - 1.) )*P[1] - 
					P_moon*(P[1]-E[1]) + Q/CONST_partMass*field[1];
	Derivs[5] = -P_plan*( part2plan_2 - CONST_planIntegrate*(5.*sinP_2 - 3.) )*P[2] - 
					P_moon*(P[2]-E[2]) + Q/CONST_partMass*field[2];

	// Enceladus equations of motion.
	Derivs[6]  = Ev[0];
	Derivs[7]  = Ev[1];	
	Derivs[8]  = Ev[2];
	Derivs[9]  = -M_plan*( moon2plan_2 - CONST_planIntegrate*(5.*sinE_2 - 1.) )*E[0];
	Derivs[10] = -M_plan*( moon2plan_2 - CONST_planIntegrate*(5.*sinE_2 - 1.) )*E[1];
	Derivs[11] = -M_plan*( moon2plan_2 - CONST_planIntegrate*(5.*sinE_2 - 3.) )*E[2];

	// Update charging equations, determine if equilibrium charge is reached. 
	Derivs[12] = Solver::phi2j(phi, vp);
}

/* Function to run modified midpoint method, to be called by Bulirsch-Stoer. */
void Solver::ModMidpoint(const double & dt, vector<double> & y_n, const int & numSteps)
{
	double h = dt/numSteps, temp_var;

	// Evaluate first step using forward Euler. 
	for(int j=0; j<CONST_numVariables; j++) { 				
		m_yTemp1[j] = y_n[j];
		m_dTemp[j]  = y_n[j];
	}
	// Call Evaluate pointer to evaluate function derivatives. Note, 'this' refers to the current
	// instance of the class Solver.
	(this->*Evaluate)(m_dTemp);

	for(int j=0; j<CONST_numVariables; j++) {
		m_yTemp2[j] = m_yTemp1[j] + h*m_dTemp[j];
		m_dTemp[j]  = m_yTemp2[j];
	}
	(this->*Evaluate)(m_dTemp);

	// Loop to run modified midpoint method with numSteps steps.
	for(int i=1; i<numSteps; i++) {
		// Update y_(i+1) = y_(i-1) + 2hf(y_i). Replace y_i with y_(i+1) and y_(i-1) with y_i.
		for(int j=0; j<CONST_numVariables; j++) { 
			temp_var    = m_yTemp1[j];
			m_yTemp1[j] = m_yTemp2[j];
			m_yTemp2[j] = temp_var + 2.*h*m_dTemp[j];
		}
		// Set d_temp = f(y_i)
		for(int j=0; j<CONST_numVariables; j++)	m_dTemp[j] = m_yTemp2[j];
		(this->*Evaluate)(m_dTemp);
	}
	// Solve for u_(t+1) = 1/2( y_n + y_(n-1) + hf(y)_n ).
	for(int j=0; j<CONST_numVariables; j++) {
		y_n[j] = .5*( m_yTemp2[j] + m_yTemp1[j] + h*m_dTemp[j] );
	}
}

/* Function to take one time step using a Bulirsch Stoer integrator. */
int Solver::BulirschStoerStep(const double & dt, double *y_n, const double & tol, bool & collideMoon)
{
	int count  = 0;
	double err = 1., Q;
	if(CONST_numVariables == 13) Q = y_n[12];

	// Initialize tableau for this step with current state. 
	for(int i=0; i<CONST_numVariables; i++) {
		m_extrapTableau[0][0][i] = y_n[i];
	}
	Solver::ModMidpoint(dt, m_extrapTableau[0][0], CONST_extrap[0]);

	while(err > tol) {
		// Increase counter
		count = count+1;

		// Break Bulirsch-Stoer iteration, reduce next step size and try again if max number of
		// extrapolation iterations has been reached. 
		if(count == CONST_maxDiv) return 0;

		// Contruct y_n copy and call ModMidpoint to operate with count steps.
		for(int i=0; i<CONST_numVariables; i++) {
			m_extrapTableau[count][0][i] = y_n[i];
		}
		Solver::ModMidpoint(dt, m_extrapTableau[count][0], CONST_extrap[count]);

		// Loop to extrapolate results. 
		for(int i=1; i<=count; i++) {
			for(int j=0; j<CONST_numVariables; j++) {
				m_extrapTableau[count][i][j] = m_extrapTableau[count][i-1][j] +
					( m_extrapTableau[count][i-1][j] - m_extrapTableau[count-1][i-1][j] ) /
					( (CONST_extrap[count]*CONST_extrap[count]) / (CONST_extrap[count-i]*CONST_extrap[count-i]) - 1. );
			}
		}
		// Save relative error vector in trivial location in tableau and compute norm.  
		for(int i=0; i<CONST_numVariables; i++) {
			m_extrapTableau[0][0][i] = ( m_extrapTableau[count][count][i] - m_extrapTableau[count][count-1][i] ) / 
				m_extrapTableau[count][count-1][i];
		}
		err = Norm(m_extrapTableau[0][0],CONST_numVariables);
	}
	// Change values in y_n to extrapolated vector y_(n+1).
	for(int i=0; i<CONST_numVariables; i++) {
		y_n[i] = m_extrapTableau[count][count][i];
	}

	// Check if equilibrium charge has been reached. If it has, switch to constant charge to
	// mass ratio. NOTE - this is only for Connerney Bfield model, not Simon. 
	if(m_bfield == 1 && CONST_numVariables == 13) {
		err = abs(Q - y_n[12]) / abs(y_n[12]);
		if(err < 1e-12) {
			CONST_numVariables = 12;
			Evaluate 		   = &Solver::EvaluateDeriv_ConstCharge;
			CONST_partRad 	  *= 1e6;
			CONST_partMass 	   = 4./3.*PI* pow(CONST_partRad,3);
			CONST_potential    = y_n[12] / (4.*PI*CONST_partRad*1e-6*8.854187*1e-12);
			Solver::SetCharge2Mass();
			m_change = 1;
		}
	}

	// Transform if y_(n+1) collided with moon.
	double dist = sqrt( (y_n[0] - y_n[6])*(y_n[0] - y_n[6]) +
						(y_n[1] - y_n[7])*(y_n[1] - y_n[7]) + 
						(y_n[2] - y_n[8])*(y_n[2] - y_n[8]) );
	if(dist < GLOBAL_radiusMoon) collideMoon = true;	
	else collideMoon = false;

	// Return if the step was successful.
	return 1;
}

/* Function to create and save a particle object for given IC, step-size and steps */
void Solver::CreateParticle(const int & steps, const double & dt, double *y, const int & data)
{	
	// Initialize Simulation 
	vector<int> locationInd;			 		   		  
	bool   collideMoon = false;	
	float  partPos_x,
		   partPos_y,
		   partPos_z;
	double dt_next,
		   curr_step;
	int    step_success;
			 		   	     
	/*-------------------------------------------------------------------------------------*/
	/*------------------------------ Saving Various Results -------------------------------*/
	/*-------------------------------------------------------------------------------------*/

	// Array to store data every numData steps. Iniitalize elements and save starting point. 
	int numData 	 = 1; 			// Store current state every numData iterations.
	int dataIndex 	 = 1;			// Index to keep track of number of states stored.
	int dataCols     = 0;
	vector<vector<double> > storeData;
	vector<double> nextState;
	vector<float> initLoc = {float(y[0]-y[6]),float(y[1]-y[7]),float(y[2]-y[8])};
	Solver::Transform(initLoc[0],initLoc[1],initLoc[2],y[6],y[7],y[8]);

	// Compute and save particle with respect to moon.
	if(data == 1) {
		dataCols = 3;
		partPos_x = y[0] - y[6];
		partPos_y = y[1] - y[7];
		partPos_z = y[2] - y[8];
		Solver::Transform(partPos_x,partPos_y,partPos_z,y[6],y[7],y[8]);
		nextState.push_back(partPos_x);
		nextState.push_back(partPos_y);
		nextState.push_back(partPos_z);	
		storeData.push_back(nextState);
	}
	// Save particle and moon locations with respect to planet.
	else if(data == 2) {
		dataCols = 6;
		nextState.push_back(y[0]);
		nextState.push_back(y[1]);
		nextState.push_back(y[2]);
		nextState.push_back(y[6]);
		nextState.push_back(y[7]);
		nextState.push_back(y[8]);
		storeData.push_back(nextState);
	}
	// Save particle charge.
	else if(data == 3) { 
		dataCols = 1;
		nextState.push_back(y[12]);
		storeData.push_back(nextState);
	}

	/*-------------------------------------------------------------------------------------*/
	/*------------------------------ Loop To Run Simulation -------------------------------*/
	/*-------------------------------------------------------------------------------------*/
	for(int i=0; i<steps; i++)
	{	
		// Take one time step of size dt. If desired error has not been reached, reduce step
		// size and repeat until we have completed one step of size dt.
		dt_next   = dt;
		curr_step = 0.;
		while( abs(dt - curr_step) > 1e-10 ) {
			step_success = Solver::BulirschStoerStep(dt_next, y, 1e-12, collideMoon);
			if(step_success > 1e-10) {
				curr_step += dt_next;
				dt_next    = dt - curr_step;
			}
			else {
				dt_next /= 2.;
			}
		}
		// Compute particle position with respect to Enceladus (in initial lat/long frame).
		partPos_x = y[0] - y[6];
		partPos_y = y[1] - y[7];
		partPos_z = y[2] - y[8];
		Solver::Transform(partPos_x,partPos_y,partPos_z,y[6],y[7],y[8]);

		// Store desired data every numData steps
		if(i % numData == 0) {
			if(data == 1) {
				nextState[0] = partPos_x;
				nextState[1] = partPos_y;
				nextState[2] = partPos_z;
				storeData.push_back(nextState);
				dataIndex = dataIndex + 1;
			}
			else if(data == 2) {
				nextState[0] = y[0];
				nextState[1] = y[1];
				nextState[2] = y[2];
				nextState[3] = y[6];
				nextState[4] = y[7];
				nextState[5] = y[8];
				storeData.push_back(nextState);
				dataIndex = dataIndex + 1;
			}
			else if(data == 3) { 
				nextState[0] = y[12];
				storeData.push_back(nextState);
				dataIndex = dataIndex + 1;
			}
		}
		// If particle collided save impact location and stop iterating.
		if( collideMoon == true ) {					
			vector<double> whereCollide = {partPos_x,partPos_y,partPos_z};
			// Normalize(whereCollide, 3, GLOBAL_radiusMoon);
			// Solver::Cartesian2Geog(whereCollide); 
			vector<float> speedCollide = {float(y[3]-y[9]),float(y[4]-y[10]),float(y[5]-y[11])};
			Solver::Transform(speedCollide[0],speedCollide[1],speedCollide[2],y[9],y[10],y[11]);
			// cout << "Particle Collided at: (" << whereCollide[0] << ", " << whereCollide[1] << ")." << endl;
			// cout << "Particle Velocity:    (" << speedCollide[0] << ", " << speedCollide[1] << ", " << speedCollide[2] << ")." << endl;
			cout << "Particle Collided at: (" << whereCollide[0] << ", " << whereCollide[1] << ", " << whereCollide[2] << ")." << endl;
			cout << "Distance from start: (" << abs(initLoc[0]-partPos_x) << "," << abs(initLoc[1]-partPos_y) << "," << abs(initLoc[2]-partPos_z) << ")" << endl;
			break;
		}
		// If particle has not collided and is within area of interest, save its location:
		if( partPos_x > CONST_min_x && partPos_x < CONST_max_x && 
			partPos_y > CONST_min_y && partPos_y < CONST_max_y &&
			partPos_z > CONST_min_z && partPos_z < CONST_max_z )
		{
			long int ind = Solver::GetDensityIndex(partPos_x,partPos_y,partPos_z);	
			locationInd.push_back(ind);
		}
	}
	double dec1 = Solver::Declination(y[6],y[7],y[8],y[9],y[10],y[11])*RAD2DEG;
	double dec2 = Solver::Declination(y[0],y[1],y[2],y[3],y[4],y[5])*RAD2DEG;
	cout << "Enc. Declination: " << dec1 << endl << "Part. Declination: " << dec2 << endl;

	// Save simulation data. 
	if(data == 1 || data == 2 || data == 3) {
		Solver::SaveData(storeData,dataIndex,dataCols);
		// Solver::SaveData(m_currents,m_currents.size(),10);
	}
	// Set CONST_numVariables and associated variables to inital status if equilibrium
	//charge was reached.
	if(m_change == 1) {
		CONST_numVariables += 1;
		CONST_partRad      *= 1e-6;
		CONST_partMass      = 4./3.*PI * pow(CONST_partRad,3) * 1000.;
		Evaluate 	   		= &Solver::EvaluateDeriv;
		m_change 			= 0;
	}	
}

/* Function to simulate a particle and return if it collided with the moon. */
void Solver::CheckCollision(const int & steps, const double & dt, double *y, const double & weight, 
							bool & collideMoon, vector<float> & whereCollide)
{
	// Initialize variables. 
	float  partPos_x,
		   partPos_y,
		   partPos_z;
	double dt_next,
		   curr_step;
	int    step_success;

	// Loop To Run Simulation.
	for(int i=0; i<steps; i++)
	{
		// Take one time step of size dt. If desired error has not been reached, reduce step
		// size and repeat until we have completed one step of size dt.
		dt_next   = dt;
		curr_step = 0.;
		while( abs(dt - curr_step) > 1e-10 ) {
			step_success = Solver::BulirschStoerStep(dt_next, y, 1e-12, collideMoon);
			if(step_success > 1e-10) {
				curr_step += dt_next;
				dt_next    = dt - curr_step;
			}
			else {
				dt_next /= 2.;
			}
		}

		// If particle collided save impact location and stop iterating.
		if( collideMoon == true ) {
			// Compute particle position with respect to Enceladus (in initial lat/long frame).
			partPos_x = y[0] - y[6];
			partPos_y = y[1] - y[7];
			partPos_z = y[2] - y[8];
			Solver::Transform(partPos_x,partPos_y,partPos_z,y[6],y[7],y[8]);

			// Save location in whereCollide vector
			whereCollide[0] = partPos_x;
			whereCollide[1] = partPos_y;
			whereCollide[2] = partPos_z;
			whereCollide[3] = y[3] - y[9];
			whereCollide[4] = y[4] - y[10];
			whereCollide[5] = y[5] - y[11];
			Solver::Transform(whereCollide[3],whereCollide[4],whereCollide[5],y[9],y[10],y[11]);
			break;
		}
	}
	// Set CONST_numVariables and associated variables to inital status if equilibrium
	//charge was reached.
	if(m_change == 1) {
		CONST_numVariables += 1;
		CONST_partRad      *= 1e-6;
		CONST_partMass      = 4./3.*PI * pow(CONST_partRad,3) * 1000.;
		Evaluate 	   		= &Solver::EvaluateDeriv;
		m_change 				= 0;
	}
}

/* Function to simulate a particle and return its declination after the given time period. */
double Solver::FinalDeclination(const int & steps, const double & dt, double *y)
{
	// Initialize variables. 
	bool collideMoon = false;	
	double dt_next,
		   curr_step;
	int step_success;

	// Loop To Run Simulation.
	for(int i=0; i<steps; i++)
	{
		// Take one time step of size dt. If desired error has not been reached, reduce step
		// size and repeat until we have completed one step of size dt.
		dt_next   = dt;
		curr_step = 0.;
		while( abs(dt - curr_step) > 1e-10 ) {
			step_success = Solver::BulirschStoerStep(dt_next, y, 1e-12, collideMoon);
			if(step_success > 1e-10) {
				curr_step += dt_next;
				dt_next    = dt - curr_step;
			}
			else {
				dt_next /= 2.;
			}
		}
		// If particle collided stop iterating.
		if( collideMoon == true ) break;
	}
	// Set CONST_numVariables and associated variables to inital status if equilibrium
	//charge was reached.
	if(m_change == 1) {
		CONST_numVariables += 1;
		CONST_partRad      *= 1e-6;
		CONST_partMass      = 4./3.*PI * pow(CONST_partRad,3) * 1000.;
		Evaluate 	   		= &Solver::EvaluateDeriv;
		m_change 				= 0;
	}
	// Compute and return particle declination.
	return Solver::Declination(y[0],y[1],y[2],y[3],y[4],y[5])*RAD2DEG;
}

/* Function to simulate a particle and track where it has been in a spatial grid about the moon.  */
/* Also saves in input reference variables a boolean telling if the particle collided and a state */
/* vector at that point in time. 																  */
void Solver::ParticleSim(const int & steps, const double & dt, double *y, const double & weight, 
	unordered_map<long int,pair<float,float> > & locationInd, bool & collideMoon,
	vector<float> & whereCollide)
{
	// Initialize variables. 
	float  partPos_x,
		   partPos_y,
		   partPos_z;
	double dt_next,
		   curr_step;
	int    step_success;

	// Loop To Run Simulation.
	for(int i=0; i<steps; i++)
	{
		// Take one time step of size dt. If desired error has not been reached, reduce step
		// size and repeat until we have completed one step of size dt.
		dt_next   = dt;
		curr_step = 0.;
		while( abs(dt - curr_step) > 1e-10 ) {
			step_success = Solver::BulirschStoerStep(dt_next, y, 1e-12, collideMoon);
			if(step_success > 1e-10) {
				curr_step += dt_next;
				dt_next    = dt - curr_step;
			}
			else {
				dt_next /= 2.;
			}
		}
		// Compute particle position with respect to Enceladus (in initial lat/long frame).
		partPos_x = y[0] - y[6];
		partPos_y = y[1] - y[7];
		partPos_z = y[2] - y[8];
		Solver::Transform(partPos_x,partPos_y,partPos_z,y[6],y[7],y[8]);

		// If particle collided save impact location and stop iterating.
		if( collideMoon == true ) {
			whereCollide[0] = partPos_x;
			whereCollide[1] = partPos_y;
			whereCollide[2] = partPos_z;
			whereCollide[3] = y[3] - y[9];
			whereCollide[4] = y[4] - y[10];
			whereCollide[5] = y[5] - y[11];
			Solver::Transform(whereCollide[3],whereCollide[4],whereCollide[5],y[9],y[10],y[11]);
			break;
		}
		// If particle has not collided and is within area of interest, save its location:
		if( partPos_x > CONST_min_x && partPos_x < CONST_max_x && 
			partPos_y > CONST_min_y && partPos_y < CONST_max_y &&
			partPos_z > CONST_min_z && partPos_z < CONST_max_z )
		{
			long int ind = Solver::GetDensityIndex(partPos_x,partPos_y,partPos_z);
			locationInd[ind].first  += weight;
			// Save charge as well if charges are being tracked. 
			if(m_charging == 2) {
				locationInd[ind].second += weight*y[12];
			}
			else if(m_charging == 1) { 
				locationInd[ind].second += weight*CONST_potential;
			}
		}
	}
	// Set CONST_numVariables and associated variables to inital status if equilibrium
	// charge was reached.
	if(m_change == 1) {
		CONST_numVariables += 1;
		CONST_partRad      *= 1e-6;
		CONST_partMass      = 4./3.*PI * pow(CONST_partRad,3) * 1000.;
		Evaluate 	   		= &Solver::EvaluateDeriv;
		m_change 			= 0;
	}		
}

