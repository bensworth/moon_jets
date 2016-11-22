#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>
#include <new>
#include "mpi.h"
#include "genFunctions.h"
#include "Solver.h"
#include "Jet.h"
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

	// Maximum change in x (and equivalently y) and z in the positive and negative
	// directions. Note, use a positive angle for these calculations.
	lat_r = abs(lat_r);
	float dxmin = gridWidth*sin(lat_r) + gridHeight_orth*cos(lat_r),
		  dxmax = max(gridHeight*cos(lat_r),gridWidth*sin(lat_r)),
		  dzmin = max(gridWidth*cos(lat_r),gridHeight*sin(lat_r)),
		  dzmax = gridWidth*cos(lat_r) + gridHeight_orth*sin(lat_r);

	// Set min, max of x,y,z w.r.t. the maximum change in x,y,z and the location of jet
	xmin = floor(x_cent - dxmin);
	xmax = ceil(x_cent + dxmax);
	ymin = floor(y_cent - dxmin);
	ymax = ceil(y_cent + dxmax);
	zmin = floor(z_cent - dzmin);
	zmax = ceil(z_cent + dzmax);
}

int main(int argc, char *argv[])
{
	//---------------------------------- Initialize parallel ----------------------------------//

	int rank,
		numProcess, 
		requestTag = 0,
		radTag 	   = 1,
		angTag 	   = 2;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcess);
	int numCore = 12 * numProcess;

			cout << rank << "/" << numProcess << endl;

	// Set system parameters for Enceladus/Saturn
	SetEuropa();

	//-------------------------------------- Master node --------------------------------------//
	if(rank == 0)
	{
		// Load .csv of cone ejection angles and number of points per circle.
		char delim;	// dummy variable to absord commas from .csv
		vector<vector<double> > coneData;
		ifstream getCone("../Data/EurPlume_Cone.csv");
		// ifstream getCone("../Data/ConeData_Const.csv");
		assert(getCone.is_open());
		while(!getCone.eof())
		{
			vector<double> temp(2,0);
			getCone >> temp[0] >> delim >> temp[1];
			coneData.push_back(temp);
		}
		getCone.close();

		// List of particle radii to simulate. 
		// vector<int> sizes = {246,261,273,284,292,300,307,313,319,329,
		// 					 339,346,353,368,380,390,399,414,460};
		// vector<int> sizes = {273,307,334,353,380,408};		
		// vector<int> sizes = {438,461,488,508};		
		vector<int> sizes = {322,368};		

		// Number of tasks to perform and number of tasks complete. 
		int N 	   = sizes.size(),
			numAng = coneData.size();

		// Hand out tasks to perform to workers until there are none left. 
		for(int i=0; i<N; i++)
		// for(int j=0; j<numAng; j++)
		{
			for(int j=numAng-1; j>=0; j--)
			// for(int i=0; i<N; i++)
			{
				int nextProc;
				MPI_Recv(&nextProc, 1, MPI_INT, MPI_ANY_SOURCE, requestTag, MPI_COMM_WORLD, &status);
				MPI_Send(&sizes[i], 1, MPI_INT, nextProc, radTag, MPI_COMM_WORLD);
				MPI_Send(&j, 1, MPI_INT, nextProc, angTag, MPI_COMM_WORLD);
				cout << "Work distributed to worker " << nextProc << " -- " << sizes[i] <<
					": " << coneData[j][0] << endl;
			}
		}

		// After all tasks have been handed out, call each worker to retire. 
		cout << "All work distributed, beginning shutdown." << endl;
		for(int i=1; i<numProcess; i++)
		{
			int nextProc, end = -1;
			MPI_Recv(&nextProc, 1, MPI_INT, MPI_ANY_SOURCE, requestTag, MPI_COMM_WORLD, &status);
			MPI_Send(&end, 1, MPI_INT, nextProc, radTag, MPI_COMM_WORLD);
		}
	}

	//-------------------------------------- Worker nodes -------------------------------------//
	else
	{
		//---------------------------- Load input/set parameters ----------------------------//

		// Define variables for simulation.	
		// Integer values
		int	codeVersion  = 8,
			jetNum 	     = 1,
			jetRef 	     = 1,
			charging 	 = 0,
			bFieldModel  = 0,
			xmin 		 = -1565,
			xmax 	     = 1565,
			ymin		 = -1565,
			ymax   		 = 1565,
			zmin 	     = -3000,
			zmax 		 = 0,
			orthogonal   = 1,
			extrapolate  = 15,
			numVariables = 12,
			numSpeeds    = 100;
		// Float values
		float orbits     = 2,
			  maxInc     = 45,
			  gridSize   = 1,
			  critRad    = .8,
	  		  gridHeight = 300,
			  gridWidth	 = 300,
			  v_gas      = .7,
	   		  volume,
			  totalTime;
		// Double values
		double errorTol = 1e-12,
		   	   pole_DEC = 64.51*DEG2RAD,
		   	   pole_RA  = 268.08*DEG2RAD;

		for(int i=1; i<argc; i++)
		{
			if(strcmp(argv[i],"-vgas") == 0)
			{
				i += 1;
				v_gas = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-RC") == 0)
			{
				i += 1;
				critRad = atof(argv[i]);
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
			else if(strcmp(argv[i],"-maxinc") == 0)
			{
				i += 1;
				maxInc = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-orb") == 0)
			{
				i += 1;
				orbits = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-orthogonal") == 0)
			{
				orthogonal = 1;
			}	
			else if(strcmp(argv[i],"-jetid") == 0)
			{
				i += 1;
				jetNum = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-numspeeds") == 0)
			{
				i += 1;
				numSpeeds = atoi(argv[i]);
			}
		}
		// Compute total time given number of orbits.
		totalTime = orbits*GLOBAL_periodMoon;
		volume    = pow(gridSize*1e+3,3);

		//------------------------------------ Load Data ------------------------------------//

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
		ifstream getCone("../Data/EurPlume_Cone.csv");
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
		ifstream getSize("../Data/Size_m.csv");
		assert(getSize.is_open());
		while(!getSize.eof())
		{
			float temp;
			getSize >> temp;					 // Load size bins in m.
			partSizes.push_back(temp*1000000.);	 // Convert to um
		}
		getSize.close();
		float d_inclination = coneData[1][0] - coneData[0][0];

		//-------------------------------- Initialize system --------------------------------//
		
		// Europa initial conditions with respect to J2000 at time Jan 1. 2014, 00:00.
		vector<double> moonPos = {-659245.02020193, -133050.42521924, -73672.92361506};
		vector<double> moonVel = {2.99536494, -11.97251661, -5.78613310};

		// Pole of Europa w.r.t. J2000 at time Jan 1. 2014, 00:00.
		vector<double> pole = { cos(pole_DEC)*cos(pole_RA),
						    	cos(pole_DEC)*sin(pole_RA),
						    	sin(pole_DEC) };

		// Transform coordinate system. 
		TransformSystem(moonPos,moonVel,pole);

		// Solve for change of basis vectors.
		vector<double> ex(3), ey(3), ez(3);
		SetChangeBasis(moonPos,ex,ey,ez);

		// Create solver.
		Solver systemSolver;
		systemSolver.SetNoCharging();							// No charging for Europa
		systemSolver.SetIntegrator(extrapolate,errorTol);
		systemSolver.SetPole(pole[0],pole[1],pole[2]);
		systemSolver.SetChangeBasisVec(ex,ey,ez);

		// Set Jet location (2 options from Roth et al).
		vector<float> jetLocation;
		jetLocation = jetLoc[jetNum];

		// Get optimal grid size as a function of jet location. 
		GetGridDims(gridHeight,gridWidth,jetLocation[0],jetLocation[1],moonPos,xmin,xmax,ymin,ymax,zmin,zmax);
		systemSolver.CreateDensityGrid(xmin,xmax,ymin,ymax,zmin,zmax,gridSize);	

		// Create Jet.
		Jet Eruptor;
		Eruptor.SetSpeedDist(v_gas,critRad);
		Eruptor.SetNumVariables(numVariables);
		Eruptor.SetLocation(jetLocation[0],jetLocation[1],jetLocation[2],jetLocation[3]);
		Eruptor.SetInitCond(moonPos,moonVel);
		Eruptor.SetMaxAngle(maxInc);
		Eruptor.SetGridData(xmin,xmax,ymin,ymax,zmin,zmax,gridSize);
		Eruptor.SetSimulationData(jetNum,codeVersion,bFieldModel,charging,jetRef,numCore);

		//-------------------------------- Update and simulate jet --------------------------------//
		int finish = 0; 
		while(finish == 0)
		{
			// Request new job from master.
			int sizeIndex, angIndex;
			MPI_Send(&rank, 1, MPI_INT, 0, requestTag, MPI_COMM_WORLD);
			MPI_Recv(&sizeIndex, 1, MPI_INT, 0, radTag, MPI_COMM_WORLD, &status);

			// If recieved work call, do said work. 
			if(sizeIndex >= 0) 
			{
				MPI_Recv(&angIndex, 1, MPI_INT, 0, angTag, MPI_COMM_WORLD, &status);
				cout << "Received work order -- " << rank << endl;

				// Simulate Jet for given particle size. 
				float inclination = coneData[angIndex][0];
				int   numAzimuth  = coneData[angIndex][1];
				Eruptor.MonteCarlo_Jet(systemSolver,numSpeeds,numAzimuth,inclination,
						d_inclination,totalTime,volume,partSizes[sizeIndex],sizeIndex);
			}
			// If recieved retire call, break loop, finalize MPI.
			else
			{
				finish = 1;
				cout << "Bedtime -- " << rank << endl;
				break;
			}
		}
	}

	if(rank == 0)
	{
		cout << "Goodnight world." << endl;
	}
	MPI_Finalize();
	return 0;
}











