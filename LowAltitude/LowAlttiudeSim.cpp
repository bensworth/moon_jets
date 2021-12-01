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
#include <omp.h>
#include "genFunctions.h"
#include "Solver.h"
#include "Jet.h"
using namespace std;
using namespace genFunctions;

int main(int argc, char *argv[])
{
    //---------------------------------- Initialize parallel ----------------------------------//
    int num_threads = omp_get_max_threads();
    cout << num_threads << " OMP threads.\n";

    // Set system parameters for Enceladus/Saturn
    SetEnceladus();

    //---------------------------- Load input/set parameters ----------------------------//
    
    double errorTol = 1e-12;
    float partPot = -1.49,
          initVel,
          partRad;
    int   extrapolate = 15, 
          charging = 2,     // Full particle charging
          numVariables,
          numAzimuth = 100,
          partRad_ind = 0,
          initVel_ind = 0,
          num_inner_inc = 4,    // Number of particles simulated/inclination bin
          numSpeeds = 28,
          numRadii = 30;
    int   orthogonal = 1,
          bFieldModel  = 1, // Connerey charging model
          jet_ind = -1;

    for(int i=1; i<argc; i++) {

        if(strcmp(argv[i],"-nocharge") == 0) {
            charging     = 0;
            numVariables = 12;
        }
        // else if(strcmp(argv[i],"-gridsize") == 0) {
        //     i += 1;
        //     gridSize = atof(argv[i]);
        // }
        else if(strcmp(argv[i],"-naz") == 0)
        {
            i += 1;
            numAzimuth = atoi(argv[i]);
        }
        else if(strcmp(argv[i],"-r") == 0)
        {
            i += 1;
            partRad_ind = atoi(argv[i]);
        }
        else if(strcmp(argv[i],"-s") == 0)
        {
            i += 1;
            initVel_ind = atoi(argv[i]);
        }
    }

    // Check for valid radii and speed indices
    if (partRad_ind > (numRadii-1)) {
        std::cout << "Invalid particle radius, must be between [0," << numRadii-1 << "].\n";
        return -1;
    }
    if (initVel_ind > (numSpeeds-1)) {
        std::cout << "Invalid particle speed, must be between [0," << numSpeeds-1 << "].\n";
        return -1;
    }

    //------------------------------------ Load Data ------------------------------------//

    // Create vector of particle speeds in km/s between [25,50,...,700]
    float dvel = 0.7/numSpeeds;
    vector<float> partSpeeds(28);
    for (int vv=0; vv<partSpeeds.size(); vv++) {
        partSpeeds[vv] = (vv+1)*dvel;
    }

    // Create vector of particle sizes in um between [0.5,1.0,...,15]
    float drad = 15.0 / numRadii;
    vector<float> partSizes(30);
    for (int vv=0; vv<partSizes.size(); vv++) {
        partSizes[vv] = (vv+1)*drad;
    }

    //-------------------------------- Initialize system --------------------------------//
    
    // Enceladus initial conditions with respect to J2000 at time Jan 1. 2014, 00:00.
    vector<double> moonPos = {13225.931, -236286.61, 16298.732};
    vector<double> moonVel = {12.612619, 0.58868893, -1.1307378};

    // Pole of Enceladus w.r.t. J2000 at time Jan 1. 2014, 00:00.
    double pole_RA      = 40.66*DEG2RAD;
    double pole_DEC     = (83.52)*DEG2RAD; 
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
    systemSolver.SetSize(partSizes[partRad_ind]);
    systemSolver.SetPlasma(moonPos[0], moonPos[1], moonPos[2]);

    // Create Jet.
    Jet Eruptor;
    vector<float> jetLocation{-90.0, 0, 0, 0};

    // Set jet parameters
    Eruptor.SetNumVariables(numVariables);
    Eruptor.SetLocation(jetLocation[0],jetLocation[1],jetLocation[2],jetLocation[3]);
    Eruptor.SetInitCond(moonPos,moonVel);

    std::cout << "Particle radius = " << partSizes[partRad_ind] << " um,\n"
        << "Initial particle speed = " << partSpeeds[initVel_ind] << " km/s\n";

    //-------------------------------- Simulate jet --------------------------------//
    Eruptor.HoverSimOMP(systemSolver, numAzimuth, partRad_ind,
        partSizes[partRad_ind], initVel_ind, partSpeeds[initVel_ind],
        num_inner_inc);

    return 0;
}











