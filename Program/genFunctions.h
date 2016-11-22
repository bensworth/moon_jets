#ifndef GENFUNCTIONSHEADERDEF
#define GENFUNCTIONSHEADERDEF
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

namespace genFunctions
{
	/* --------------------------- Declare constants --------------------------- */
	extern double PI,
				  DEG2RAD, 
				  RAD2DEG;
	// Planet radius, gravitational constant, J2 moment, distance to sun, planetary 
	// Bfield corotation angular rate (Omega_p).
	extern double GLOBAL_radiusPlan, 
				  GLOBAL_muPlan, 
				  GLOBAL_J2, 
				  GLOBAL_plan2sun, 
				  GLOBAL_Omega_p;
	// Moon radius, semi-major axis, gravitational constant and orbital period.
	extern double GLOBAL_radiusMoon,
				  GLOBAL_moonSMA,
				  GLOBAL_muMoon, 
				  GLOBAL_periodMoon,
				  GLOBAL_radiusPlume;
	extern int 	  GLOBAL_bodyID;

	/* --------------------------- Declare functions --------------------------- */
	void 		   SetEnceladus();
	void 		   SetEuropa();
	double 		   NormP(const double *A, const int & size);
	double		   Norm(const vector<double> & A, const int & size = 3);
	void		   Normalize(vector<double> & A, const int & size = 3, const double & r = 1);
	double 		   Dot(const vector<double> & A, const vector<double> & B, const int & size = 3);
	vector<double> Cross(const vector<double> & A, const vector<double> & B);
	double 		   Distance(const vector<double> & A, const vector<double> & B, const int & size = 3);
	vector<double> QuatMult(const vector<double> & Q1, const vector<double> & Q2);
	void 		   TransformSystem(vector<double> & pos, vector<double> & vel, vector<double> & pole);
	void 		   SetChangeBasis(const vector<double> encPos, vector<double> & ex, vector<double> & ey,
					vector<double> & ez);
	void 		   InitializeDust(double *currentState,  vector<double> ex,  vector<double> ey, 
		 			vector<double> ez, const float & initVel, const float & latitude,
					const float & longitude, const float &inc = 0, const float &azimuth = 0);


	template<typename T>
	ostream& operator<< (ostream& out, const vector<T>& v) {
	    out << "[";
	    size_t last = v.size() - 1;
	    for(size_t i = 0; i < v.size(); ++i) {
	        out << v[i];
	        if (i != last)
	            out << ", ";
	    }
	    out << "]";
	    return out;
	}
}

#endif
