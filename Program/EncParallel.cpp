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
#include <cassert>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>
#include <new>
#include "mpi.h"
#include <omp.h>
#include "genFunctions.h"
#include "Solver.h"
#include "Jet.h"
using namespace std;
using namespace genFunctions;

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
	int num_threads = omp_get_max_threads();

			cout << rank << "/" << numProcess << ", OMP = " << num_threads << endl;

	// Set system parameters for Enceladus/Saturn
	SetEnceladus();

	//-------------------------------------- Master node --------------------------------------//
	if(rank == 0) {
		// Load .csv of cone ejection angles and number of points per circle.
		char delim;	// dummy variable to absord commas from .csv
		vector<vector<double> > coneData;
		ifstream getCone("../Data/EncConeData.csv");
		// ifstream getCone("../Data/EncConeData_Const.csv");
		assert(getCone.is_open());

		while(!getCone.eof()) {
			vector<double> temp(2,0);
			getCone >> temp[0] >> delim >> temp[1];
			coneData.push_back(temp);
		}
		getCone.close();

		// List of particle radii to simulate. 
		// vector<int> sizes = {273,313,343,353,381,438,462};
		// vector<int> sizes = {339,348,361,373,380,391,400,414,427,432,437};
		vector<int> sizes = {372,395,411,427,432};

		// Number of tasks to perform and number of tasks complete. 
		int N = sizes.size();
		int numAng = coneData.size();

		// Hand out tasks to perform to workers until there are none left. 
		for(int i=0; i<N; i++)
		{
			for(int j=numAng-1; j>=0; j--)
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
		for(int i=1; i<numProcess; i++) {
			int nextProc, end = -1;
			MPI_Recv(&nextProc, 1, MPI_INT, MPI_ANY_SOURCE, requestTag, MPI_COMM_WORLD, &status);
			MPI_Send(&end, 1, MPI_INT, nextProc, radTag, MPI_COMM_WORLD);
		}
	}

	//-------------------------------------- Worker nodes -------------------------------------//
	else {
		//---------------------------- Load input/set parameters ----------------------------//
		
		// Constant parameters for the system.
		double errorTol    = 1e-12;

		// Parameters controlling simulation. 
		float partPot  = -1.49,		orbits = 2.,	maxInc = 15.,
			  gridSize = 5.,		critRad = 0.2,  v_gas = 0.5,
			  volume,	  			totalTime;

		int   oneOrbit = 118442,	extrapolate  = 15, 
			  xmin	   = -250,		xmax		 = 250, 	ymin		 = -250,
			  ymax	   = 250,		zmin 		 = -1000, 	zmax 		 = 500,
			  charging = 2,			numVariables = 13,		numSpeeds 	 = 100;

		// Information about current simulation for output purposes.
		int   orthogonal = 0,		jetNum 	    = 1,	bFieldModel  = 1,
			  jetRef 	 = 1,		codeVersion = 10,	jet_ind 	 = -1;

		bool collisionOnly = false, montecarlo = false;

		for(int i=1; i<argc; i++) {

			if(strcmp(argv[i],"-nocharge") == 0) {
				charging     = 0;
				numVariables = 12;
			}
			else if(strcmp(argv[i],"-gridsize") == 0) {
				i += 1;
				gridSize = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-xmin") == 0) {
				i += 1;
				xmin = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-xmax") == 0) {
				i += 1;
				xmax = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-ymin") == 0) {
				i += 1;
				ymin = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-ymax") == 0) {
				i += 1;
				ymax = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-zmin") == 0) {
				i += 1;
				zmin = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-zmax") == 0) {
				i += 1;
				zmax = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-maxinc") == 0) {
				i += 1;
				maxInc = atof(argv[i]);
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
			else if(strcmp(argv[i],"-plasma") == 0) {
				i += 1;
				bFieldModel = atoi(argv[i]);
			}	
			else if(strcmp(argv[i],"-jetref") == 0) {
				i += 1;
				jetRef = atoi(argv[i]);
			}
			else if(strcmp(argv[i],"-collisions") == 0) {
				collisionOnly = true;
				if(rank == 1) {
					cout << "Collisions only\n";
				}
			}
			else if(strcmp(argv[i],"-montecarlo") == 0) {
				montecarlo = true;
			}
			else if(strcmp(argv[i],"-vgas") == 0)
			{
				i += 1;
				v_gas = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-RC") == 0)
			{
				i += 1;
				critRad = atof(argv[i]);
			}
			else if(strcmp(argv[i],"-numspeeds") == 0)
			{
				i += 1;
				numSpeeds = atoi(argv[i]);
			}
		}
		totalTime = ceil(orbits*oneOrbit);
		volume    = pow(gridSize*1e+3,3);

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

		// Load .csv of cone ejection angles and number of points per circle.
		vector<vector<float> > coneData;
		ifstream getCone("../Data/EncConeData.csv");
		// ifstream getCone("../Data/EncConeData_Const.csv");
		assert(getCone.is_open());

		while(!getCone.eof()) {
			vector<float> temp(2,0);
			getCone >> temp[0] >> delim >> temp[1];
			coneData.push_back(temp);
		}
		getCone.close();
		float d_inclination = coneData[1][0] - coneData[0][0];

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

		//-------------------------------- Initialize system --------------------------------//
		
		// Enceladus initial conditions with respect to J2000 at time Jan 1. 2014, 00:00.
		vector<double> moonPos = {13225.931, -236286.61, 16298.732};
		vector<double> moonVel = {12.612619, 0.58868893, -1.1307378};

		// Pole of Enceladus w.r.t. J2000 at time Jan 1. 2014, 00:00.
		double pole_RA  	= 40.66*DEG2RAD;
		double pole_DEC 	= (83.52)*DEG2RAD; 
		vector<double> pole = { cos(pole_DEC)*cos(pole_RA),
						    	cos(pole_DEC)*sin(pole_RA),
						    	sin(pole_DEC) };

		// Transform coordinate system. 
		TransformSystem(moonPos,moonVel,pole);

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
		systemSolver.SetIntegrator(extrapolate,errorTol);
		systemSolver.SetPole(pole[0],pole[1],pole[2]);
		systemSolver.CreateDensityGrid(xmin,xmax,ymin,ymax,zmin,zmax,gridSize);	

		// Create Jet.
		Jet Eruptor;
		for(int i=0; i<jet_IDs.size(); i++) {
			if(jetNum == jet_IDs[i]) {
				jet_ind = i;
			}
		}
		if(jet_ind < 0) {
			cout << "Error - jet not found. \n";
			return -1; 
		}

		vector<float> jetLocation = jetLoc[jet_ind];
		if(orthogonal == 1) {
			jetLocation[2] = 0;
			jetLocation[3] = 0;
		}

		if(rank == 1) {
			std::cout << "Jet ID - " << jetNum << ", loc - (" << jetLocation[0] << ", " << jetLocation[1] << ", " << jetLocation[2] << ", " << jetLocation[3] << ")\n";
		}

		// Set speed distribution if running monte Carlo simulation
		if (montecarlo) {
			Eruptor.SetSpeedDist(v_gas,critRad);	
		}
		Eruptor.SetNumVariables(numVariables);
		Eruptor.SetLocation(jetLocation[0],jetLocation[1],jetLocation[2],jetLocation[3]);
		Eruptor.SetInitCond(moonPos,moonVel);
		Eruptor.SetMaxAngle(maxInc);
		Eruptor.SetGridData(xmin,xmax,ymin,ymax,zmin,zmax,gridSize);
		Eruptor.SetSimulationData(jetNum,codeVersion,bFieldModel,charging,jetRef,numCore);

		//-------------------------------- Update and simulate jet --------------------------------//
		int finish = 0; 
		while(finish == 0) {
			// Request new job from master.
			int sizeIndex, angIndex;
			MPI_Send(&rank, 1, MPI_INT, 0, requestTag, MPI_COMM_WORLD);
			MPI_Recv(&sizeIndex, 1, MPI_INT, 0, radTag, MPI_COMM_WORLD, &status);

			// If recieved work call, do said work. 
			if(sizeIndex >= 0)  {
				MPI_Recv(&angIndex, 1, MPI_INT, 0, angTag, MPI_COMM_WORLD, &status);
			
				// Finalize solver initialization. 
				systemSolver.SetSize(partSizes[sizeIndex]);
				systemSolver.SetPlasma(moonPos[0], moonPos[1], moonPos[2]);

				// Simulate Jet for given particle size. 
				if (collisionOnly) {
					cout << "Received work order -- " << rank << ": size ind = " <<
							sizeIndex << ", ang ind = " << angIndex << std::endl;
					Eruptor.CollisionMap(systemSolver, partSpeeds, partWeights[sizeIndex], coneData[angIndex],
										 totalTime, volume, partSizes[sizeIndex], sizeIndex);
				}
				else if (montecarlo) {
					cout << "Received work order -- " << rank << ": rc = " << critRad << ", vgas = " << vgas <<
							", size ind = " << sizeIndex << ", ang ind = " << angIndex << std::endl;
					float inclination = coneData[angIndex][0];
					int   numAzimuth  = coneData[angIndex][1];
					Eruptor.MonteCarlo_Jet(systemSolver, numSpeeds, numAzimuth, inclination, d_inclination,
										   totalTime, volume, partSizes[sizeIndex], sizeIndex);
				} 
				else {
					cout << "Received work order -- " << rank << ": size ind = " <<
							sizeIndex << ", ang ind = " << angIndex << std::endl;
					Eruptor.SpecSimOMP(systemSolver, partSpeeds, partWeights[sizeIndex], coneData[angIndex],
									   totalTime, volume, partSizes[sizeIndex], sizeIndex);
				}
			}
			// If recieved retire call, break loop, finalize MPI.
			else {
				finish = 1;
				cout << "Bedtime -- " << rank << endl;
				break;
			}
		}
	}
	if(rank == 0) {
		cout << "Goodnight world." << endl;
	}
	MPI_Finalize();
	return 0;
}











