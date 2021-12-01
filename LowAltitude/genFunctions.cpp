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

#include <vector>
#include <cmath>
#include "genFunctions.h"
using namespace std;

namespace genFunctions{

	double PI 	   = 3.141592653589793;
	double DEG2RAD = PI/180.0,
		   RAD2DEG = 180.0/PI;
	int    GLOBAL_bodyID     = 0;
	// Planetary constants
	double GLOBAL_radiusPlan = 0,
		   GLOBAL_muPlan     = 0,
		   GLOBAL_J2         = 0,
		   GLOBAL_plan2sun   = 0, 
		   GLOBAL_Omega_p    = 0;
	// Moon constants
	double GLOBAL_radiusMoon  = 0,
		   GLOBAL_moonSMA     = 0,
		   GLOBAL_muMoon      = 0,
		   GLOBAL_periodMoon  = 0,
		   GLOBAL_radiusPlume = 0;

	void SetEnceladus()
	{
		GLOBAL_bodyID = 1;

		// Planet (equatorial) radius, gravitational constant, J2 moment, distance to sun,  
		// planetary Bfield corotation angular rate (Omega_p).
		GLOBAL_radiusPlan = 60330.;
		GLOBAL_muPlan  	  = 3.7931*1e+7;
		GLOBAL_J2  		  = .016298;
		GLOBAL_plan2sun   = 1373880928.32861920;	// This is only for plasma -- Saturn specific?
		GLOBAL_Omega_p    = 10.785;

		// Moon (equatorial) radius, semi-major axis, gravitational constant and orbital period.
		GLOBAL_radiusMoon  = 249.4;
		GLOBAL_moonSMA	   = 237837.39;
		GLOBAL_muMoon  	   = 7.204291351372671;
		GLOBAL_periodMoon  = 2.*PI*sqrt( GLOBAL_moonSMA*GLOBAL_moonSMA*GLOBAL_moonSMA/GLOBAL_muPlan );
		GLOBAL_radiusPlume = .8*GLOBAL_radiusMoon;
	}

	void SetEuropa()
	{
		GLOBAL_bodyID = 2;

		// Planet (equatorial) radius, gravitational constant, J2 moment.
		GLOBAL_radiusPlan = 71492.;
		GLOBAL_muPlan  	  = 1.26712767881e+8;
		GLOBAL_J2  		  = .014736;

		// Moon (equatorial) radius, semi-major axis, gravitational constant and orbital period.
		GLOBAL_radiusMoon = 1564.13;
		GLOBAL_moonSMA	  = 670900;
		GLOBAL_muMoon  	  = 3202.935;
		GLOBAL_periodMoon = 2.*PI*sqrt( GLOBAL_moonSMA*GLOBAL_moonSMA*GLOBAL_moonSMA/GLOBAL_muPlan );
	}

	double NormP(const double *A, const int & size)
	{	
		double N = 0.;
		for(int i=0; i<size; i++)
		{
			N += A[i]*A[i];
		}
		return sqrt(N);
	}

	double Norm(const vector<double> & A, const int & size)
	{
		double N = 0.;
		for(int i=0; i<size; i++)
		{
			N += A[i]*A[i];
		}
		return sqrt(N);
	}

	void Normalize(vector<double> & A, const int & size, const double & r)
	{
		double N = 0.;
		for(int i=0; i<size; i++)
		{
			N += A[i]*A[i];
		}
		N = sqrt(N);
		for(int i=0; i<size; i++)
		{
			A[i] = A[i]*r/N;
		}
	}
	double Dot(const vector<double> & A, const vector<double> & B, const int & size)
	{
		double D = 0.;
		for(int i=0; i<size; i++)
		{
			D += A[i]*B[i];
		}
		return D;
	}

	vector<double> Cross(const vector<double> & A, const vector<double> & B)
	{
		vector<double> C(3);
		C[0] = A[1]*B[2] - A[2]*B[1];
		C[1] = A[2]*B[0] - A[0]*B[2];
		C[2] = A[0]*B[1] - A[1]*B[0];
		return C;
	}

	double Distance(const vector<double> & A, const vector<double> & B, const int & size)
	{
		double N =0.;
		for(int i=0; i<size; i++)
		{
			N += (A[i] - B[i])*(A[i] - B[i]);
		}
		return sqrt(N);
	}

	// Get 4d indices for flattened 1d array. Must have
	//		0<=i<n1, 0<=j<n2, 0<=k<n3, 0<=l<n4
	// with array dimensions n1 x n2 x n3 x n4 
	int get4dind(int i, int j, int k, int l, int n1, int n2, int n3, int n4)
	{
	    return i*n2*n3*n4 + j*n3*n4 + k*n4 + l;
	}

	// Get 4d indices for flattened 1d array. Must have
	//		0<=i<n1, 0<=j<n2,
	// with array dimensions n1 x n2
	int get2dind(int i, int j, int n1, int n2)
	{
	    return i*n2 + j;
	}

	/* Function to perform quaternion multiplication on quaternions Q1 = (w1,V1), */
	/* Q2 = (w2,V2), for angles w1,w2 and vector parts V1,V2. */
	vector<double> QuatMult(const vector<double> & Q1, const vector<double> & Q2)
	{
		// Take cross product of vector parts of quaternions. 
		vector<double> V1(3);
			V1[0] = Q1[1];
			V1[1] = Q1[2];
			V1[2] = Q1[3];
		vector<double> V2(3);
			V2[0] = Q2[1];
			V2[1] = Q2[2];
			V2[2] = Q2[3];
		V1 = Cross(V1,V2);
		
		// Form product Q = Q1*Q2.
		vector<double> Q(4);
		Q[0] = Q1[0]*Q2[0] - ( Q1[1]*Q2[1] + Q1[2]*Q2[2] + Q1[3]*Q2[3] );
		Q[1] = Q1[0]*Q2[1] + Q2[0]*Q1[1] + V1[0];
		Q[2] = Q1[0]*Q2[2] + Q2[0]*Q1[2] + V1[1];
		Q[3] = Q1[0]*Q2[3] + Q2[0]*Q1[3] + V1[2];
		return Q;
	}
	/* Transform system to be in a Saturn inertial frame, with Enceladus orbiting in the xy-plane. */
	void TransformSystem(vector<double> & pos, vector<double> & vel, vector<double> & pole)
	{
		// Angular momentum vector.
		vector<double> r = pos, 
					   v = vel,
					   angMom;
		Normalize(r);
		Normalize(v);
		angMom = Cross(r,v);

		// Form matrix rows to orthogonalize angular momentum vector.
		vector<double> row1, row2, row3;
		double n1 = sqrt( pow(angMom[0],2) + pow(angMom[1],2) ), 
			   n2 = Norm(angMom);
		row1 = { angMom[1]/n1,
				 -angMom[0]/n1,
				 0. };
	    row2 = { angMom[0]*angMom[2]/(n2*n1),
	   			 angMom[1]*angMom[2]/(n2*n1),
	   			 -n1/n2 };
	    row3 = { angMom[0]/n2,
	   	 		 angMom[1]/n2,
	   	 		 angMom[2]/n2 };

		// Transfer Enceladus system to have the xy-plane be the orbital plane
		// Note, this is just matrix multiplication.
		double x0, y0, z0;
		x0 		= Dot(pos,row1);
		y0 		= Dot(pos,row2);
		z0 		= Dot(pos,row3);
		pos[0] 	= x0; 
		pos[1] 	= y0;
		pos[2] 	= z0;
		x0 		= Dot(vel,row1);
		y0 		= Dot(vel,row2);
		z0 		= Dot(vel,row3);
		vel[0] 	= x0;
		vel[1] 	= y0;
		vel[2] 	= z0;	
		x0 		= Dot(pole,row1);
		y0 		= Dot(pole,row2);
		z0		= Dot(pole,row3);
		pole[0] = x0;
		pole[1] = y0;
		pole[2] = z0;

		// Rotate initial (x,y) coordinate of Enceladus in orbital plane to be (x,0).
		double denom, temp;
		denom   = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
		
		temp 	= (vel[0]*pos[0] + vel[1]*pos[1])/denom;
		vel[1]  = (-vel[0]*pos[1] + vel[1]*pos[0])/denom;
		vel[0]  = temp;
		temp 	= (pole[0]*pos[0] + pole[1]*pos[1])/denom;
		
		pole[1] = (-pole[0]*pos[1] + pole[1]*pos[0])/denom;
		pole[0] = temp;
		
		pos[0]  = (pos[0]*pos[0] + pos[1]*pos[1])/denom;
		pos[1]  = 0.;
	}

	/* Function to create columns of change of basis matrix from Enceladus latitude/longitude */
	/* frame to Saturn intertial frame. */
	void SetChangeBasis(const vector<double> moonPos, vector<double> & ex, vector<double> & ey,
		vector<double> & ez)
	{
		ez[0] = 0.;
		ez[1] = 0.;
		ez[2] = 1.;
		ey = Cross(ez,moonPos);
		Normalize(ey);
		ex = Cross(ey,ez);
		Normalize(ex);
	}

	/* Function that takes in a pointer to the current system state, initial particle velocity, and */ 
	/* longitude and latitude of the Jet. The first six elements of the pointer variable represent  */
	/* the particle position and velocity, which are initially equal to Enceladus' position and     */
	/* velocity. Function sets the initial particle position and velocity. 							*/ 
	void InitializeDust(double *currentState,  vector<double> ex,  vector<double> ey, 
		 vector<double> ez, const float & initVel, const float & latitude,
		const float & longitude, const float &inc, const float &azimuth)
	{
		vector<double> moonPos, jetPos, jetDir, rotationDir, vx, vy, vz;

		double lon   = DEG2RAD*longitude - PI,
			   lat 	 = DEG2RAD*latitude,
			   phi   = -azimuth*DEG2RAD,   // because it's western longitude
	    	   theta = PI/2.-inc*DEG2RAD;  // zenith -> theta

	    // Basis for moon position.   
		moonPos = { currentState[6], currentState[7], currentState[8] };
		ez[0] = 0.;
		ez[1] = 0.;
		ez[2] = 1.;
		ey = Cross(ez,moonPos);
		Normalize(ey);
		ex = Cross(ey,ez);
		Normalize(ex);
		
		// Transferring from longitude/latitude to Cartesian and then to the new Cartesian
		// system can be done as follows:
		jetPos = { cos(lon)*cos(lat)*ex[0] + sin(lon)*cos(lat)*ey[0],
				  	  cos(lon)*cos(lat)*ex[1] + sin(lon)*cos(lat)*ey[1],
			 	  	  sin(lat) };
		Normalize(jetPos);

		// Set initial dust position.
		currentState[0] += GLOBAL_radiusMoon*jetPos[0];
		currentState[1] += GLOBAL_radiusMoon*jetPos[1];
		currentState[2] += GLOBAL_radiusMoon*jetPos[2];

		// Account for Enceladus rotational speed. 
		rotationDir = Cross(ez,jetPos);
		Normalize(rotationDir);
		double vRot = 2.*PI*cos(lat)*GLOBAL_radiusMoon / GLOBAL_periodMoon;
		
	    // Ejection direction of the jet (including tilt)
	    vz = jetPos;
	    vx = {0.,0.,1.};
	    vy = Cross(vz, vx);
	    Normalize(vy);
	    vx = Cross(vy, vz);
	    Normalize(vx);
	    jetDir = { cos(phi)*cos(theta)*vx[0] + sin(phi)*cos(theta)*vy[0] + sin(theta)*vz[0],
	    		   cos(phi)*cos(theta)*vx[1] + sin(phi)*cos(theta)*vy[1] + sin(theta)*vz[1],
	    		   cos(phi)*cos(theta)*vx[2] + sin(phi)*cos(theta)*vy[2] + sin(theta)*vz[2] };

		// Set initial dust velocity.
		currentState[3] += initVel*jetDir[0] + vRot*rotationDir[0];
		currentState[4] += initVel*jetDir[1] + vRot*rotationDir[1];
		currentState[5] += initVel*jetDir[2] + vRot*rotationDir[2];
	}
}
