#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>
#include <vector>
#include <random>
#include <cmath>
#include <new>
#include "genFunctions.h"
#include "Solver.h"
#include "Jet.h"
#ifdef _OPENMP
	#include <omp.h>
#endif
using namespace std;
using namespace genFunctions;

/* Function to take in vectors of Enceladus position and velocity, save them as member */
/* then compute the initial particle location and velocity due to Enceladus' velocity and */
/* Enceladus' rotation. Also stores member unit vector initPos with the directional */
/* vector of the particle launch location. */
void GetInitCond(vector<double> &moonPos, vector<double> &moonVel,
	vector<double> &parPos, vector<double> &dirVec, vector<double> &parNonEjVel, 
	double &lat, double &lon, const float & partAzimuth, const float & partInc)
{
	double vRot, phi, theta, jetZenith = 0.0, jetAzimuth = 0.0;

	vector<double> ex, ey, ez, vx, vy, vz, jetPos, jetDir;

	// Convert latitude/longitude to radians for this code.
	lon = DEG2RAD*lon;
	lat = DEG2RAD*lat;
	
	// Put longitude and latitude into Cartesian unit direction vector. Note, these
	// basis vectors rotate from planet centered frame to moon centered frame. 
	ez = {0., 0., 1.};
	ey = Cross(ez, moonPos);
	Normalize(ey);
	ex = Cross(ey, ez);
	Normalize(ex);
	// Since transformation is orthogonal, taking transpose we can trasnfer from
	// moon frame to planet fram - Unit vector of jet's location in inertial frame.
	jetPos = {cos(lon)*cos(lat)*ex[0] + sin(lon)*cos(lat)*ey[0],
	          cos(lon)*cos(lat)*ex[1] + sin(lon)*cos(lat)*ey[1],
			  sin(lat)};
	Normalize(jetPos);

	// Set initial dust position.
	parPos[0] = moonPos[0] + GLOBAL_radiusMoon*jetPos[0];
	parPos[1] = moonPos[1] + GLOBAL_radiusMoon*jetPos[1];
	parPos[2] = moonPos[2] + GLOBAL_radiusMoon*jetPos[2];

	// Account for Enceladus rotational speed. 
	vector<double> rotationDir = Cross(ez, jetPos);
	Normalize(rotationDir);
	vRot = 2.*PI*cos(lat)*GLOBAL_radiusMoon / GLOBAL_periodMoon;

	// Set initial dust velocity.
	parNonEjVel[0] = moonVel[0] + vRot*rotationDir[0];
	parNonEjVel[1] = moonVel[1] + vRot*rotationDir[1];
	parNonEjVel[2] = moonVel[2] + vRot*rotationDir[2];

    // Ejection direction of the jet (including tilt)
    vz = jetPos;
    vx = {0., 0., 1.};
    vy = Cross(vz, vx);
    	Normalize(vy);
    vx = Cross(vy, vz);
    	Normalize(vx);
    phi   = -jetAzimuth*DEG2RAD;        // because it's western longitude
    theta = PI/2. - jetZenith*DEG2RAD;  // zenith -> theta
    jetDir = { cos(phi)*cos(theta)*vx[0] + sin(phi)*cos(theta)*vy[0] + sin(theta)*vz[0],
    		     cos(phi)*cos(theta)*vx[1] + sin(phi)*cos(theta)*vy[1] + sin(theta)*vz[1],
    		     cos(phi)*cos(theta)*vx[2] + sin(phi)*cos(theta)*vy[2] + sin(theta)*vz[2] };

	// build orthonormal system for the jet in inertial frame, where jetDir 
	// is parallel to +z and +x is in apex direction.
    vz = jetDir;
    if(Dot(ex, vz) < 1.) vx = ex;
    else vx = ez;
    vy = Cross(vz, vx);
    	Normalize(vy);
    vx = Cross(vy, vz);
   		Normalize(vx);

   	// Get particle ejection direction. 
	double partTheta = PI/2. - partInc*DEG2RAD; // Because we must have longitude in [-180,180]
	double partPhi   = -partAzimuth*DEG2RAD;

	dirVec[0] = cos(partPhi)*cos(partTheta)*vx[0] + sin(partPhi)*cos(partTheta)*vy[0] + sin(partTheta)*vz[0];
	dirVec[1] = cos(partPhi)*cos(partTheta)*vx[1] + sin(partPhi)*cos(partTheta)*vy[1] + sin(partTheta)*vz[1];
	dirVec[2] = cos(partPhi)*cos(partTheta)*vx[2] + sin(partPhi)*cos(partTheta)*vy[2] + sin(partTheta)*vz[2];
}


/* Function to calculate the declination of a particle launched from the south pole of */
/* Enceladus for initial velocities minVel:velStep:maxVel. Saves declinations in .csv file. */
void DeclinationSpeed(Solver & systemSolver, const int & numVariables, 
	const vector<double> & moonPos, const vector<double> & moonVel, const int & simSteps,
	const double & dt, const float & minVel, const float & maxVel, const float & velStep,
	const float & partRad)
{
	int loops 			   = ceil((maxVel - minVel)/velStep);
	float *partDeclination = new float[loops+1];
	float *velocities 	   = new float[loops+1];

	#pragma omp parallel for
	for(int i=0;i<=loops;i++)
	{
		// Set current velocity
		velocities[i] = minVel + i*velStep;

		Solver tempSolver = systemSolver;
		double *y = new double[numVariables];
		// Reset initial conditions
		y[0]  = moonPos[0];
		y[1]  = moonPos[1];
		y[2]  = moonPos[2] - GLOBAL_radiusMoon;
		y[3]  = moonVel[0];
		y[4]  = moonVel[1];
		y[5]  = moonVel[2] - velocities[i];
		y[6]  = moonPos[0];
		y[7]  = moonPos[1];
		y[8]  = moonPos[2];
		y[9]  = moonVel[0];
		y[10] = moonVel[1];
		y[11] = moonVel[2];
		if(numVariables == 13)
		{
			y[12] = 0.;
		}

		// Have systemSolver object create and simulate a particle.
		partDeclination[i] = tempSolver.FinalDeclination(simSteps, dt, y);
		delete [] y;
	}

	// Save declination and velocity as .csv file.
	stringstream filename;
	if(GLOBAL_bodyID == 1) { 
		filename << "../Results/Declination/Enc_Declination_r" << partRad << ".csv";
	}
	else if(GLOBAL_bodyID == 2) { 
		filename << "../Results/Declination/Eur_Declination_r" << partRad << ".csv";
	}
	ofstream output(filename.str());
	for(int p=0; p<loops+1; p++)
	{
		output << velocities[p] << ", " << partDeclination[p] << endl;  
	}
	output.close();

	// Free pointers.
	delete [] velocities;
	delete [] partDeclination;
}

/* Find the escape speed from particle launched near the south pole. Finds escape speed for */
/* numPoints launch locations around numRings different rings centered on the south pole,   */
/* equispaced from the south pole to maxInc degrees from the pole. 							*/
void EscapeSpeed(Solver & systemSolver, const int & numVariables,
	const vector<double> & moonPos, const vector<double> & moonVel,
	const int & simSteps, const double & dt, const int & numRings, const int & numPoints,
	float & maxInc, const float & partRad)
{
	// Number of spatial rings centered at south pole as a list of latitudes. 
	maxInc  	     = maxInc*DEG2RAD;
	float *latitude  = new float[numRings];
	double latSep    = maxInc/numRings;
	for(int i=0; i<numRings; i++){ latitude[i] = -PI/2. + latSep*(i+1); }
	
	// Create double pointer to store data, vector to keep current angles. Note that we increase 
	// the number of launch sites by 7 for each ring further from the moon. Thus we can allocate 
	// space as  numLaunch = numRings*numPoints + 7*SUM_(i=0)^(numRings-1) i, where
	// 7*SUM_(i=0)^(numRings-1) i = numRings*(numRings-1)*7/2 
	int numLaunch = numRings*numPoints + 7./2.*numRings*(numRings - 1) + 1;
	vector<vector<double> > escVel(numLaunch);

	// Solve for change of basis vectors.
	vector<double> ex(3);
	vector<double> ey(3);
	vector<double> ez(3);
	SetChangeBasis(moonPos, ex, ey, ez);

	// Loops to find particle escape speed along rings.
	#pragma omp parallel for
	for(int i=0; i<numRings; i++)
	{
		cout << i << endl;

		// Number of points around each spatial ring as a list of longitudes. Increase number
		// of points as rings get larger (try to maintain density). 
		Solver tempSolver = systemSolver;
		double *y = new double[numVariables];
		vector<double> currentAngle(3);
		int numPole       = numPoints+7*i;
		float *longitude  = new float[numPole];
		double longSep    = 2.*PI/numPole;
		for(int k=0; k<numPole; k++) longitude[k] = longSep*k;

		for(int j=0; j<numPole; j++)
		{
			// Set initial collision boolean to true and initial minimum particle velocity.
			bool collision = true; 
			float initVel   = .145;
			
			// Loop over increasing velocities until particle escapes. 
			while(collision == true)
			{
				initVel += .005;
				// Reset Enceladus initial conditions, particle launch location, and increase velocity.
				y[0]  = moonPos[0];
				y[1]  = moonPos[1];
				y[2]  = moonPos[2];
				y[3]  = moonVel[0];
				y[4]  = moonVel[1];
				y[5]  = moonVel[2];
				y[6]  = moonPos[0];
				y[7]  = moonPos[1];
				y[8]  = moonPos[2];
				y[9]  = moonVel[0];
				y[10] = moonVel[1];
				y[11] = moonVel[2];
				if(numVariables == 13)
				{
					y[12] = 0.;
				}
				InitializeDust(y, ex, ey, ez, initVel, latitude[i], longitude[j]);

				// Simulate particle 
				collision = tempSolver.CheckCollision(simSteps, dt, y);
			}
			
			// Save current launch location (x,y) and launch velocity in 3 element vector.
			int ind = i*numPoints + 7./2.*i*(i - 1) + j;
			currentAngle[0] = cos(latitude[i])*cos(longitude[j]);
			currentAngle[1] = cos(latitude[i])*sin(longitude[j]);
			currentAngle[2] = initVel;
			escVel[ind] = currentAngle;
		}
		delete [] longitude;
		delete [] y;
	}

	// Test the South Pole as well.
	bool collision = true;
	float initVel = .145;
	double *y = new double[numVariables];
	vector<double> currentAngle(3);

	// Loop over increasing velocities until particle escapes.
	while(collision == true)
	{
		// Reset Enceladus initial conditions, particle launch location, and increase velocity.
		initVel = initVel + .005;
		y[0]  = moonPos[0];
		y[1]  = moonPos[1];
		y[2]  = moonPos[2];
		y[3]  = moonVel[0];
		y[4]  = moonVel[1];
		y[5]  = moonVel[2];
		y[6]  = moonPos[0];
		y[7]  = moonPos[1];
		y[8]  = moonPos[2];
		y[9]  = moonVel[0];
		y[10] = moonVel[1];
		y[11] = moonVel[2];
		if(numVariables == 13)
		{
			y[12] = 0.;
		}		
		InitializeDust(y, ex, ey, ez, initVel, -PI/2., 0.);
		// Simulate particle 
		collision = systemSolver.CheckCollision(simSteps, dt, y);
	}
	currentAngle[0] = 0.;
	currentAngle[1] = 0.;
	currentAngle[2] = initVel;
	escVel[numLaunch-1] = currentAngle;
	
	// Save escape velocity pointer as .csv file. 
	stringstream filename;
	if(GLOBAL_bodyID == 1) { 
		filename << "../Results/EscapeSpeed/Enc_EscapeSpeed_r" << partRad << ".csv";
	}
	else if(GLOBAL_bodyID == 2) { 
		filename << "../Results/EscapeSpeed/Eur_EscapeSpeed_r" << partRad << ".csv";
	}
	ofstream output(filename.str());
	for(int i=0; i<escVel.size(); i++)
	{
		output << escVel[i][0] << ", " << escVel[i][1] << ", " << escVel[i][2] << endl;  
	}
	output.close();

	// Free pointers.
	delete [] latitude;
	delete [] y;
}

int main(int argc, char *argv[])
{
	if(argc == 1)
	{
		cout << "\nPlease include an optional initial flag specifying particle charging: \n" <<
				"\t -nocharge \n \t -constcharge \n \t -charge (default) \n" <<
				"followed by a mandatory flag specifying simulation: " << endl << endl;
		cout << "-sim: Launch one particle orthogonal to surface and save data. Takes optional flags: \n" <<
			" \t -size (particle size in um; default 1) \n" <<
			" \t -orb  (number of orbits to simulate; default 2) \n" << 
			" \t -vel  (initial particle velocity in km/s; default .25) \n" <<  
			" \t -long (launch longitude on Enceladus in degree; default 0) \n" <<
			" \t -ang  (angle of launch from Enceladus pole in degrees; default 0) \n" <<
			" \t -save (1: save particle w.r.t. Enceladus, 2: save particle and Enceladus \n" <<
			" \t        w.r.t. Saturn, 0: do not save data; default 1)" << endl << endl;
		cout << "-dec: Save particle declination vs. initial velocity. Takes optional flags: \n" <<
			" \t -size   (particle size in um; default 1) \n" <<
			" \t -orb    (number of orbits to simulate; default 2) \n" << 
			" \t -minvel (minimum initial velocity in km/s; default .2) \n" <<  
			" \t -maxvel (maximum initial velocity in km/s; default .35) \n" <<
			" \t -step   (step size of velocities over [minvel,maxvel]; default .005)" << endl << endl;
		cout << "-esc: Save particle escape speeds for particles launched orthogonal to the surface \n" << 
			"in rings centered about Enceladus' south pole. Takes optional flags: \n" <<
			" \t -size   (particle size in um; default 1) \n" <<
			" \t -orb    (number of orbits to simulate; default 2) \n" << 
			" \t -maxang (maximum angle in degrees from pole to launch particles; default 60) \n" <<  
			" \t -rings  (number of rings to launch particles from; default 14) \n" <<
			" \t -points (number of points in first ring; default 40)" << endl << endl;
		cout << "-jet: Simulate one inclination of a jet for a given particle size. Saves density profile\n" << 
			" in grid about Enceladus and state vectors of particles before collision. Takes optional flags: \n" <<
			" \t -jetid      (jet id to set jet location from Porco et al.; default 1) \n" <<  
			" \t -sizeid     (index for particle size from bins; default 300) \n" <<
			" \t -inc 	     (inclination of particle launch in degrees; default 0) \n" <<
			" \t -dphi 	     (dphi degrees for phi in [0,360]; default 5) \n" <<
			" \t -orthogonal (simulates jet orthogonal to surface) \n" << 
			" \t -orb        (number of orbits to simulate; default 2)" << endl << endl;
			return 0;
	}

	// --------------------------------------------------------------------------- //
	// ---------------------------- Set up simulation ---------------------------- //
	// --------------------------------------------------------------------------- //

	// Set Enceladus global constants.
	SetEnceladus();

	// Define variables for simulation.	
	// Integer values		
	int	codeVersion = 10,		bFieldModel  = 1,		jetRef 	    = 1,
		numCore     = 1,		charging 	 = 2,		jetNum      = 1,
		data 	    = 1,		xmin 		= -250,	
		xmax 		= 250,		ymin		 = -250,	ymax 	 	= 250,
		zmin 		= -1000,	zmax 		 = 500,		numRings    = 14,
		numPoints    = 40,		sizeIndex    = 300, 	orthogonal  = 0,
		extrapolate  = 15, 		choice       = 0,  		numVariables,
		jet_ind 	= -1, 		simSteps;
	// Float values
	float gridSize    = 2.5,		volume       = pow(gridSize,3)*1e+3,
		  partPot 	  = -1.49,	inclination  = 0.,		dphi 	    = 5.,
		  maxInc 	  = 15.,	orbits  	 = 2.,		initVel 	= .25,
		  angle 	  = 0., 	longitude    = 0.,		minVel 	    = .15,
		  maxVel 	  = .45,	velStep 	 = .005,	maxAng 	    = 60.,
		  partRad 	  = 1,		totalTime;
	 // Double values
	double pole_RA    = 40.66*DEG2RAD, 		errorTol = 1e-12,
		   pole_DEC   = 83.52*DEG2RAD,		dt 		 = 5.;	

	// Enceladus initial conditions with respect to J2000 at time Jan 1. 2014, 00:00.
	vector<double> moonPos = {13225.931, -236286.61,  16298.732};
	vector<double> moonVel = {12.612619, 0.58868893, -1.1307378};

	// Pole of Enceladus w.r.t. J2000 at time Jan 1. 2014, 00:00.
	vector<double> pole = { cos(pole_DEC)*cos(pole_RA),
				 	 		cos(pole_DEC)*sin(pole_RA),
			     		    sin(pole_DEC) };

	// Take input on charging equations. If no input provided, use full charging equations.
    int firstInput = 2;
	if(strcmp(argv[1],"-constcharge") == 0)		charging = 1;
	else if(strcmp(argv[1],"-nocharge") == 0)	charging = 0;
	else if(strcmp(argv[1],"-charge") == 0)		charging = 2;
	else 
	{
		charging   = 2; 
		firstInput = 1;
	}

	// Transform coordinate system. 
	TransformSystem(moonPos,moonVel,pole);

	// Solve for change of basis vectors.
	vector<double> ex(3), ey(3), ez(3);
	SetChangeBasis(moonPos, ex, ey, ez);

	// Create solver.
	Solver systemSolver;
	systemSolver.SetCharge(partPot);
	if(charging == 0) {
		numVariables = 12;
		systemSolver.SetNoCharging();
	}	   
	else if(charging == 1) {
		numVariables = 12;
		systemSolver.SetConstCharge();
	}
	else {
		numVariables = 13;
		systemSolver.SetCharging();
	}  
	systemSolver.SetBfield(bFieldModel);
	systemSolver.SetIntegrator(extrapolate, errorTol);
	systemSolver.SetPole(pole[0],pole[1],pole[2]);
	systemSolver.SetChangeBasisVec(ex,ey,ez);
	systemSolver.CreateDensityGrid(xmin,xmax,ymin,ymax,zmin,zmax,gridSize);	

	// ---------------------------------------------------------------------------- //
	// --------------------------- Read input arguments --------------------------- //
	// ---------------------------------------------------------------------------- //

	if(strcmp(argv[firstInput],"-sim") == 0)
	{
		choice    = 1;
		for(int i=2; i<argc; i++)
		{
			if(strcmp(argv[i],"-size") == 0) {
				i += 1;
				partRad = atof(argv[i]);
			}
			if(strcmp(argv[i],"-bfield") == 0) {
				i += 1;
				bFieldModel = atoi(argv[i]);
				systemSolver.SetBfield(bFieldModel);
			}
			else if(strcmp(argv[i],"-orb") == 0) {
				i += 1;
				orbits = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-vel") == 0) {
				i += 1;
				initVel = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-long") == 0) {
				i += 1;
				longitude = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-ang") == 0) {
				i += 1;
				angle = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-save") == 0) {
				i += 1;
				data = atoi(argv[i]);
			}
		}
	}
	else if(strcmp(argv[firstInput],"-dec") == 0)
	{
		choice  = 2;
		for(int i=2; i<argc; i++)
		{
			if(strcmp(argv[i],"-size") == 0) {
				i += 1;
				partRad = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-orb") == 0) {
				i += 1;
				orbits = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-maxvel") == 0) {
				i += 1;
				maxVel = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-minvel") == 0) {
				i += 1;
				minVel = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-step") == 0) {
				i += 1;
				velStep = atof(argv[i]);
			}
		}

	}
	else if(strcmp(argv[firstInput],"-esc") == 0)
	{
		choice    = 3;
		for(int i=2; i<argc; i++)
		{
			if(strcmp(argv[i],"-size") == 0) {
				i += 1;
				partRad = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-maxang") == 0) {
				i += 1;
				maxAng = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-orb") == 0) {
				i += 1;
				orbits = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-rings") == 0) {
				i += 1;
				numRings = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-points") == 0) {
				i += 1;
				numPoints = atoi(argv[i]);
			}
		}
	}
	else if(strcmp(argv[firstInput],"-jet") == 0)
	{
		choice      = 4;	
		for(int i=2; i<argc; i++)
		{
			if(strcmp(argv[i],"-sizeid") == 0) {
				i += 1;
				sizeIndex = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-orb") == 0) {
				i += 1;
				orbits = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-jetid") == 0) {
				i += 1;
				jetNum = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-orthogonal") == 0) {
				orthogonal = 1;
			}
			else if(strcmp(argv[i],"-inc") == 0) {
				i += 1;
				inclination = atof(argv[i]);;
			}
			else if(strcmp(argv[i],"-dphi") == 0) {
				i += 1;
				dphi = atof(argv[i]);;
			}								
		}
	}
	// Compute total timesteps given dt and number of orbits.
	totalTime = orbits*GLOBAL_periodMoon;
	simSteps  = ceil(totalTime/dt);

	//---------------------------------------------------------------------------------------------//
	//--------------------------------------- Single solve ----------------------------------------//
	//---------------------------------------------------------------------------------------------//
	if(choice == 1)
	{
		// Set particle size and plasma environment if charging equations are used. 
		systemSolver.SetSize(partRad);
		if(charging == 2) {
			systemSolver.SetPlasma(moonPos[0], moonPos[1], moonPos[2]);
		}


		double latitude = -90.0,
			   longitude = 0.0,
			   partAzimuth = 10.0,
			   partInc     = 10.0;
		vector<double> ejecDir(3), parPos(3), parNonEjVel(3);
		GetInitCond(moonPos, moonVel, parPos, ejecDir, parNonEjVel,
					latitude, longitude, partAzimuth, partInc);

		double *y = new double[numVariables];
		y[0]  = parPos[0];
		y[1]  = parPos[1];
		y[2]  = parPos[2];
		y[3]  = parNonEjVel[0] + initVel*ejecDir[0];
		y[4]  = parNonEjVel[1] + initVel*ejecDir[1];
		y[5]  = parNonEjVel[2] + initVel*ejecDir[2];
		y[6]  = moonPos[0];
		y[7]  = moonPos[1];
		y[8]  = moonPos[2];
		y[9]  = moonVel[0];
		y[10] = moonVel[1];
		y[11] = moonVel[2];
		if(numVariables == 13) {
			y[12] = 0;
		}

		// cout <<y[0]<<", "<<y[1]<<", "<<y[2]<<", "<<y[3]<<", "<<y[4]<<", "<<y[5]<< endl; 
		// cout <<y[6]<<", "<<y[7]<<", "<<y[8]<<", "<<y[9]<<", "<<y[10]<<", "<<y[11]<< endl << endl; 

		// vector<double> M = {y[6],y[7],y[8]}, B0, B1;
		// double rTemp = sqrt(pow(y[0],2)+pow(y[1],2)+pow(y[2],2));
		// B0 = systemSolver.Bfield_Connerney(y[0],y[1],y[2],M,rTemp); 
		// B1 = systemSolver.Bfield_Simon(y[0],y[1],y[2],M,rTemp); 
		// cout << "Connerney: " << B0[0] << ", " << B0[1] << ", " << B0[2] << endl;
		// cout << "Simon: 	" << B1[0] << ", " << B1[1] << ", " << B1[2] << endl;

		// systemSolver.CreateParticle(simSteps, dt, y, data);
		// double potential = y[12] / (4.*PI*partRad*1e-6*8.854187*1e-12);
		// cout << "Potential: " << potential << ", Charge: " << y[12] << endl;

		systemSolver.CheckCollision(simSteps, dt, y);

		delete [] y;

	}
	//---------------------------------------------------------------------------------------------//
	//---------------------------------- Declination vs. speed ------------------------------------//
	//---------------------------------------------------------------------------------------------//
	else if(choice == 2)
	{
		// Set particle size and plasma environment if charging equations are used. 
		systemSolver.SetSize(partRad);
		if(charging == 2) {
			systemSolver.SetPlasma(moonPos[0],moonPos[1],moonPos[2]);
		}
		DeclinationSpeed(systemSolver,numVariables,moonPos,moonVel,
				simSteps,dt,minVel,maxVel,velStep,partRad);
	}
	//---------------------------------------------------------------------------------------------//
	//--------------------------------------- Escape Speed ----------------------------------------//
	//---------------------------------------------------------------------------------------------//
	else if(choice == 3)
	{
		// Set particle size and plasma environment if charging equations are used. 
		systemSolver.SetSize(partRad);
		if(charging == 2) {
			systemSolver.SetPlasma(moonPos[0],moonPos[1],moonPos[2]);
		}
		EscapeSpeed(systemSolver,numVariables,moonPos,moonVel,
				simSteps,dt,numRings,numPoints,maxAng,partRad);
	}
	//---------------------------------------------------------------------------------------------//
	//--------------------------------------- Jet Density -----------------------------------------//
	//---------------------------------------------------------------------------------------------//
	else if(choice == 4)
	{
		//------------------------------------ Load Data ------------------------------------//

		// Load .csv of jet locations and associated numeric IDs.
		char delim;	// dummy variable to absord commas from .csv
		vector<vector<float> > jetLoc;
		vector<int> jet_IDs;
		ifstream getJet("../Data/EncJetSources.csv");
		assert(getJet.is_open());
		while(!getJet.eof()) {
			float tempD;
			vector<float> tempV(4,0);
			getJet >> tempD >> delim >> tempV[0] >> delim >> tempV[1] >> delim >> tempV[2] >> delim >> tempV[3];
			jet_IDs.push_back(tempD);
			jetLoc.push_back(tempV);
		}
		getJet.close();

		// Load .csv of particle speeds.
		vector<float> partSpeeds;
		ifstream getSpeed("../Data/Speed_ms.csv");
		assert(getSpeed.is_open());
		while(!getSpeed.eof()) {
			float temp;
			getSpeed >> temp; 						 // Load speed bins in m/s
			partSpeeds.push_back( temp / 1000.);	 // Convert to km/s
		}
		getSpeed.close();

		// Load .csv of particle sizes.
		vector<float> partSizes;
		ifstream getSize("../Data/Size_m_reduced.csv");
		assert(getSize.is_open());
		while(!getSize.eof()) {
			float temp;
			getSize >> temp;					 // Load size bins in m.
			partSizes.push_back(temp*1000000.);	 // Convert to um
		}
		getSize.close();

		// Load .csv of weights for all particle speeds and sizes into a vector of vectors.
		vector<vector<double> > partWeights;
		ifstream getWeight("../Data/Weights_Normalized_1e-3.csv");
		assert(getWeight.is_open());
		while(!getWeight.eof()) {
			// Define temporary variables, get next line in .csv array.
			vector<double> temp;
			string line, tempS;
			double tempD;
			size_t pos;
			getline(getWeight,line);	

			// Pick out numeric substrings between commas, convert to double and save in vector.
			while( (pos=line.find(",")) != string::npos ) {	
				tempS = line.substr(0,pos);
				tempD = atof(tempS.c_str());
				temp.push_back(tempD);
				line.erase(0,pos+1);
			}
			// Get last element of line, pass current line into vector of vectors. 
			tempD = atof(tempS.c_str());
			temp.push_back(tempD);
			partWeights.push_back(temp);
		}
		getWeight.close();

		//----------------------------------------------------------------------------------//

		// Set particle size and plasma environment if charging equations are used. 
		partRad   = partSizes[sizeIndex];
		systemSolver.SetSize(partRad);
		if(charging == 2) {
			systemSolver.SetPlasma(moonPos[0], moonPos[1], moonPos[2]);
		}
		
		// Test Jet location
		for(int i=0; i<jet_IDs.size(); i++) {
			if(jetNum == jet_IDs[i]) {
				jet_ind = i;
			}
		}
		if(jet_ind < 0) {
			cout << "Error - jet not found. \n";
			return -1; 
		}



		vector<float> angleData   = {inclination, float(360./dphi)};
		vector<float> jetLocation = jetLoc[jet_ind];

		cout  << "jet @ " << jetLocation << endl;

		if(orthogonal == 1) {
			jetLocation[2] = 0;
			jetLocation[3] = 0;
		}

		// Create Jet object.
		Jet Eruptor;
		Eruptor.SetLocation(jetLocation[0],jetLocation[1],jetLocation[2],jetLocation[3]);
		Eruptor.SetNumVariables(numVariables);
		Eruptor.SetInitCond(moonPos,moonVel);
		Eruptor.SetMaxAngle(maxInc);
		Eruptor.SetGridData(xmin,xmax,ymin,ymax,zmin,zmax,gridSize);
		Eruptor.SetSimulationData(jetNum,codeVersion,bFieldModel,charging,jetRef,numCore);

		// Call to simulate 
		// Eruptor.SpecSimOMP(systemSolver, partSpeeds, partWeights[sizeIndex], angleData,
		// 		totalTime, volume, partSizes[sizeIndex], sizeIndex);



	}



	cout << endl;
	return 0;
}


















