#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>
#include <vector>
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


void GetGridDims(const float &gridHeight, const float &gridWidth, const float &lat,
	const float &lon, const vector<double> &moonPos, int &xmin, int &xmax, int &ymin,
	int &ymax, int &zmin, int &zmax)
{
	float gridHeight_orth = GLOBAL_radiusMoon - sqrt(abs(pow(GLOBAL_radiusMoon,2)-pow(gridWidth,2)));

	// Convert latitude/longitude to radians.
	float lon_r = DEG2RAD*(lon-180.),
		  lat_r = DEG2RAD*lat;
	
	// Put longitude and latitude into Cartesian unit direction vector. Note, these
	// basis vectors rotate from planet centered frame to moon centered frame. 
	vector<double> ex,ey,ez, jetPos;
	ez = {0., 0., 1.};
	ey = Cross(ez, moonPos);
	Normalize(ey);
	ex = Cross(ey, ez);
	Normalize(ex);

	// Since transformation is orthogonal, taking transpose we can trasnfer from
	// moon frame to planet frame - Unit vector of jet's location in inertial frame.
	jetPos = {cos(lon_r)*cos(lat_r)*ex[0] + sin(lon_r)*cos(lat_r)*ey[0],
	          cos(lon_r)*cos(lat_r)*ex[1] + sin(lon_r)*cos(lat_r)*ey[1],
			  sin(lat_r)};
	Normalize(jetPos);

	// Set initial dust position in Jovian frame.
	float x_cent = GLOBAL_radiusMoon*jetPos[0],
		  y_cent = GLOBAL_radiusMoon*jetPos[1],
		  z_cent = GLOBAL_radiusMoon*jetPos[2];

cout << gridHeight_orth << endl;

	// Maximum change in x (and equivalently y) and z in the positive and negative
	// directions. Note, use a positive angle for these calculations.
	lat_r = abs(lat_r);
	float dxmin = gridWidth*sin(lat_r) + gridHeight_orth*cos(lat_r),
		  dxmax = max(gridHeight*cos(lat_r),gridWidth*sin(lat_r)),
		  dzmin = max(gridWidth*cos(lat_r),gridHeight*sin(lat_r)),
		  dzmax = gridWidth*cos(lat_r) + gridHeight_orth*sin(lat_r);

cout << "(" << x_cent << "," << y_cent << "," << z_cent << ")" << endl;
cout << dxmin << ", " << dxmax << ", " << dzmin << ", " << dzmax << endl;

	// Set min, max of x,y,z w.r.t. the maximum change in x,y,z and the location of jet
	xmin = floor(x_cent - dxmin);
	xmax = ceil(x_cent + dxmax);
	ymin = floor(y_cent - dxmin);
	ymax = ceil(y_cent + dxmax);
	zmin = floor(z_cent - dzmin);
	zmax = ceil(z_cent + dzmax);
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
	float *latitude  = new float[numRings];
	double latSep    = maxInc/numRings;
	for(int i=0; i<numRings; i++){ latitude[i] = -90. + latSep*(i+1); }
	
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
		double longSep    = 360./numPole;
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
			currentAngle[0] = cos(DEG2RAD*latitude[i])*cos(DEG2RAD*longitude[j]);
			currentAngle[1] = cos(DEG2RAD*latitude[i])*sin(DEG2RAD*longitude[j]);
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
		InitializeDust(y, ex, ey, ez, initVel, -90., 0.);
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
				"\t -nocharge (only option right now)\n \t -constcharge \n \t -charge \n" <<
				"followed by a mandatory flag specifying simulation: " << endl << endl;
		cout << "-sim: Launch one particle orthogonal to surface and save data. Takes optional flags: \n" <<
			" \t -size     (particle size in um; default 1) \n" <<
			" \t -orb 	   (number of orbits to simulate; default 2) \n" << 
			" \t -vel  	   (initial particle velocity in km/s; default 0.25) \n" <<  
			" \t -long 	   (launch longitude on Europa in degree; default 0) \n" <<
			" \t -ang 	   (angle of launch from Europa pole in degrees; default 0) \n" <<
			" \t -gridsize (size (km) of single cube in density cube) \n" <<
			" \t -save 	   (1: save particle w.r.t. Europa, 2: save particle and Europa \n" <<
			" \t        w.r.t. Saturn, 0: do not save data; default 1)" << endl << endl;
		cout << "-dec: Save particle declination vs. initial velocity. Takes optional flags: \n" <<
			" \t -size   (particle size in um; default 1) \n" <<
			" \t -orb    (number of orbits to simulate; default 2) \n" << 
			" \t -minvel (minimum initial velocity in km/s; default .2) \n" <<  
			" \t -maxvel (maximum initial velocity in km/s; default .35) \n" <<
			" \t -step   (step size of velocities over [minvel,maxvel]; default .005)" << endl << endl;
		cout << "-esc: Save particle escape speeds for particles launched orthogonal to the surface \n" << 
			"in rings centered about Europa's south pole. Takes optional flags: \n" <<
			" \t -size   (particle size in um; default 1) \n" <<
			" \t -orb    (number of orbits to simulate; default 2) \n" << 
			" \t -maxang (maximum angle in degrees from pole to launch particles; default 60) \n" <<  
			" \t -rings  (number of rings to launch particles from; default 14) \n" <<
			" \t -points (number of points in first ring; default 40)" << endl << endl;
			return 0;
	}
	// FIX DEFAULT VALUES


	// --------------------------------------------------------------------------- //
	// ---------------------------- Set up simulation ---------------------------- //
	// --------------------------------------------------------------------------- //

	// Set Europa global constants.
	SetEuropa(); 

	// Define variables for simulation.	
	// Integer values
	int	codeVersion = 10,		charging    = 0,		data 	    = 1,
		xmin 		= -4000,	xmax 		= 4000, 	ymin 	    = -4000,	
		ymax 		= 4000,	 	zmin 	    = -7500,	zmax 		= 0,
		numRings    = 14,		numPoints   = 40,		orthogonal  = 0,
		extrapolate = 15,		choice      = 0,		jetNum 		= 1,
		sizeIndex   = 273, 		angIndex    = 0, 		bFieldModel = 0,
		jetRef 		= 1,		numCore 	= 1,		numVariables,
		simSteps;
	// Float values
	float gridSize    = 1.5,	orbits      = 2.,		initVel    = 0.5,
		  angle 	  = 0., 	longitude   = 0.,		maxAng 	   = 60,
		  minVel 	  = 0,		maxVel 	    = 0,		velStep    = 0,
		  partRad	  = 1, 		maxInc 		= 15, 		gridHeight = 300,
		  gridWidth	 = 300,		inc 		= 0,		azimuth    = 0,
		  volume,				totalTime;
	// Double values
	double pole_RA  = 268.08*DEG2RAD,	errorTol = 1e-12,
	   	   pole_DEC = 64.51*DEG2RAD, 	dt       = 2.5;

	// Europa initial conditions with respect to J2000 at time Jan 1. 2014, 00:00.
	vector<double> moonPos = {-659245.02020193, -133050.42521924, -73672.92361506};
	vector<double> moonVel = {2.99536494, -11.97251661, -5.78613310};

	// Pole of Europa w.r.t. J2000 at time Jan 1. 2014, 00:00.
	vector<double> pole = { cos(pole_DEC)*cos(pole_RA),
			     			cos(pole_DEC)*sin(pole_RA),
			     			sin(pole_DEC) };

	// Take input on charging equations. If no input provided, use no particle charge.
    int firstInput = 2;
	if(strcmp(argv[1],"-constcharge") == 0)
	{
		cout << "Only no charge for Europa." << endl;
		return 0;
	}
	else if(strcmp(argv[1],"-charge") == 0)	
	{
		cout << "Only no charge for Europa." << endl;
		return 0;
	}
	else if(strcmp(argv[1],"-nocharge") == 0)	charging = 0;
	else 
	{
		charging   = 0; 
		firstInput = 1;
	}

	// Transform coordinate system. 
	TransformSystem(moonPos,moonVel,pole);

	// Solve for change of basis vectors.
	vector<double> ex(3), ey(3), ez(3);
	SetChangeBasis(moonPos, ex, ey, ez);

	// Create solver.
	Solver systemSolver;
	if(charging == 0)
	{
		numVariables = 12;
		systemSolver.SetNoCharging();
	}	   
	else if(charging == 1)
	{
		// numVariables = 12;
		// systemSolver.SetConstCharge();
		// systemSolver.SetCharge(partPot);
	}
	else
	{
		// numVariables = 13;
		// systemSolver.SetCharging();
	}  
	systemSolver.SetIntegrator(extrapolate, errorTol);
	systemSolver.SetPole(pole[0],pole[1],pole[2]);

	xmin = -1000;
	xmax = -800;
	ymin = -100;
	ymax = 100;
	zmin = -2200;
	zmax = -1200;

	systemSolver.CreateDensityGrid(xmin,xmax,ymin,ymax,zmin,zmax,gridSize);	

	// ---------------------------------------------------------------------------- //
	// --------------------------- Read input arguments --------------------------- //
	// ---------------------------------------------------------------------------- //

	if(strcmp(argv[firstInput],"-sim") == 0)
	{
		choice    = 1;
		for(int i=2; i<argc; i++)
		{
			if(strcmp(argv[i],"-size") == 0)
			{
				i += 1;
				partRad = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-orb") == 0)
			{
				i += 1;
				orbits = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-vel") == 0)
			{
				i += 1;
				initVel = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-long") == 0)
			{
				i += 1;
				longitude = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-ang") == 0)
			{
				i += 1;
				angle = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-inc") == 0)
			{
				i += 1;
				inc = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-azimuth") == 0)
			{
				i += 1;
				azimuth = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-save") == 0)
			{
				i += 1;
				data = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-gridsize") == 0)
			{
				i += 1;
				gridSize = atof(argv[i]);
			}
		}
	}
	else if(strcmp(argv[firstInput],"-dec") == 0)
	{
		choice  = 2;
		for(int i=2; i<argc; i++)
		{
			if(strcmp(argv[i],"-size") == 0)
			{
				i += 1;
				partRad = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-orb") == 0)
			{
				i += 1;
				orbits = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-maxvel") == 0)
			{
				i += 1;
				maxVel = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-minvel") == 0)
			{
				i += 1;
				minVel = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-step") == 0)
			{
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
			if(strcmp(argv[i],"-size") == 0)
			{
				i += 1;
				partRad = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-maxang") == 0)
			{
				i += 1;
				maxAng = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-orb") == 0)
			{
				i += 1;
				orbits = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-rings") == 0)
			{
				i += 1;
				numRings = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-points") == 0)
			{
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
			if(strcmp(argv[i],"-sizeid") == 0)
			{
				i += 1;
				sizeIndex = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-orb") == 0)
			{
				i += 1;
				orbits = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-jetid") == 0)
			{
				i += 1;
				jetNum = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-orthogonal") == 0)
			{
				orthogonal = 1;
			}
			else if(strcmp(argv[i],"-angid") == 0)
			{
				i += 1;
				angIndex = atoi(argv[i]);;
			}
			else if(strcmp(argv[i],"-gridsize") == 0)
			{
				i += 1;
				gridSize = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-height") == 0)
			{
				i += 1;
				gridHeight = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-width") == 0)
			{
				i += 1;
				gridWidth = atof(argv[i]);
			}							
		}
	}
	// Compute total time given number of orbits.
	totalTime = orbits*GLOBAL_periodMoon;
	simSteps  = ceil(totalTime/dt);
	volume    = pow(gridSize*1e+3,3);
	dt = 0.5 * gridSize / initVel;

	//---------------------------------------------------------------------------------------------//
	//--------------------------------------- Single solve ----------------------------------------//
	//---------------------------------------------------------------------------------------------//
	if(choice == 1)
	{
		// Set particle size and plasma environment if charging equations are used. 
		systemSolver.SetSize(partRad);
		if(charging == 2)
		{
			systemSolver.SetPlasma(moonPos[0], moonPos[1], moonPos[2]);
		}
		
		double *y = new double[numVariables];
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
			y[12] = 0;
		}
		InitializeDust(y,ex,ey,ez,initVel,angle - 90.,longitude,inc,azimuth);
		systemSolver.CreateParticle(simSteps, dt, y, data);
		delete [] y;

	}
	//---------------------------------------------------------------------------------------------//
	//---------------------------------- Declination vs. speed ------------------------------------//
	//---------------------------------------------------------------------------------------------//
	else if(choice == 2)
	{
		// Set particle size and plasma environment if charging equations are used. 
		systemSolver.SetSize(partRad);
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
		EscapeSpeed(systemSolver,numVariables,moonPos,moonVel,
				simSteps,dt,numRings,numPoints,maxAng,partRad);
	}
	//---------------------------------------------------------------------------------------------//
	//--------------------------------------- Jet Density -----------------------------------------//
	//---------------------------------------------------------------------------------------------//
	else if(choice == 4)
	{
		// Load .csv of jet locations and associated numeric IDs.
		char delim;	// dummy variable to absord commas from .csv
		vector<vector<float> > jetLoc;
		vector<int> jetID;
		ifstream getJet("../Data/EurJetSources.csv");
		assert(getJet.is_open());
		while(!getJet.eof()) {
			float tempD;
			vector<float> tempV(4,0);
			getJet >> tempD >> delim >> tempV[0] >> delim >> tempV[1] >> delim >> tempV[2] >> delim >> tempV[3];
			jetID.push_back(tempD);
			jetLoc.push_back(tempV);
		}
		getJet.close();

		// Load .csv of cone ejection angles and number of points per circle.
		vector<vector<float> > coneData;
		ifstream getCone("../Data/ConeData.csv");
		// ifstream getCone("../Data/ConeData_Const.csv");
		assert(getCone.is_open());
		while(!getCone.eof())
		{
			vector<float> temp(2,0);
			getCone >> temp[0] >> delim >> temp[1];
			coneData.push_back(temp);
		}
		getCone.close();

		// Load .csv of particle sizes.
		vector<float> partSizes;
		ifstream getSize("../Data/Size_m_reduced.csv");
		assert(getSize.is_open());
		while(!getSize.eof())
		{
			float temp;
			getSize >> temp;					 // Load size bins in m.
			partSizes.push_back(temp*1000000.);	 // Convert to um
		}
		getSize.close();
		float d_inclination = coneData[1][0] - coneData[0][0];

		// Set Jet location (2 options from Roth et al).
		vector<float> jetLocation;
		jetLocation = jetLoc[jetNum];

		// Get optimal grid size as a function of jet location. 
		GetGridDims(gridHeight, gridWidth, jetLocation[0],jetLocation[1],moonPos,xmin,xmax,ymin,ymax,zmin,zmax);
		systemSolver.CreateDensityGrid(xmin,xmax,ymin,ymax,zmin,zmax,gridSize);	

		cout << "(" << xmin << "," << xmax << ") , (" << ymin << "," << ymax << "), (" << zmin << "," << zmax << ") - " << gridSize << endl;

		// Create Jet.
		float v_gas   = .7,
			  critRad = .8;
		Jet Eruptor;
		Eruptor.SetSpeedDist(v_gas,critRad);
		Eruptor.SetNumVariables(numVariables);
		Eruptor.SetLocation(jetLocation[0],jetLocation[1],jetLocation[2],jetLocation[3]);
		Eruptor.SetInitCond(moonPos,moonVel);
		Eruptor.SetMaxAngle(maxInc);
		Eruptor.SetGridData(xmin,xmax,ymin,ymax,zmin,zmax,gridSize);
		Eruptor.SetSimulationData(jetNum,codeVersion,bFieldModel,charging,jetRef,numCore);
	}	
	

	cout << endl;
	return 0;
}


















