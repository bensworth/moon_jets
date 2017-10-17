#ifndef SOLVERHEADERDEF
#define SOLVERHEADERDEF
#include <unordered_map>
#include <vector>
using namespace std;
using namespace genFunctions;
class Solver
{
public:

	/*-------------------------------------------------------------------------------------*/
	/*-------------------------------------- Methods --------------------------------------*/
	/*-------------------------------------------------------------------------------------*/

	Solver();  
	~Solver();

	void   SetNoCharging();
	void   SetConstCharge();
	void   SetCharging();
	void   SetBfield(const int & bfield = 1);
	void   SetIntegrator(const int & extrapolations, const double & error);
	void   SetPole(const double & x, const double & y, const double & z);
	void   SetCharge(const float & c);
	void   SetSize(const float & r, const float & rho = 1.);
	void   SetCharge2Mass();
	void   SetPlasma(const double & x, const double & y, const double & z);
	void   CreateDensityGrid(const int & minX, const int & maxX, 
			const int & minY, const int & maxY, const int & minZ, const int & maxZ,
			const float & gridSize);
	void   CreateParticle(const int & steps, const double & dt, double *y, const int & data);
	void   CheckCollision(const int & steps, const double & dt, double *y, const double & weight, 
			 bool & collideMoon, vector<float> & whereCollide);
	double FinalDeclination(const int & steps, const double & dt, double *y);
	void   ParticleSim(const int & steps, const double & dt, double *y, const double & weight, 
			unordered_map<long int,pair<float,float> > & locationInd, bool & collideMoon,
			vector<float> & whereCollide);

	vector<double> Bfield_Connerney(const double & x, const double & y, const double & z, 
			const vector<double> & moonPos, const double & r);
	vector<double> Bfield_Simon(const double & x, const double & y, const double & z, 
			const vector<double> & moonPos, const double & r);

private:

	vector<vector<double> > m_currents;


	// Equation of motion constants
	double CONST_poleX,
		   CONST_poleY,
		   CONST_poleZ,
		   CONST_planIntegrate;
	int    CONST_numVariables,
		   m_change; 
	// Vectors to be used in solving steps of the equations. 
	int    CONST_maxDiv;
	double CONST_errTol;
	vector<vector<vector<double> > > m_extrapTableau;
	vector<double> m_yTemp1,
				   m_yTemp2,
				   m_dTemp,
				   CONST_extrap;
	// Fixed charge with assumed potential CONST_potential, and particle radius CONST_partRad.
	float  CONST_potential,
		   CONST_partRad,
		   CONST_partMass;
	double CONST_Qd;
	// Charging equations constants
	double m_bfield,
		   m_charging,
		   CONST_SpecialF_Xec,
		   CONST_SpecialF_Xeh,
		   CONST_SpecialF5_Xec,
		   CONST_SpecialF5_Xeh;
   	// Saturn B-field coefficients - Connerney (1993). 
	double CONST_g10,
		   CONST_g20,
		   CONST_g30;
	// Saturn Bfield coonstants - Simon (2013) - see Enceladus_Bfield_poynting.nb.
	double CONST_mu0,
		   CONST_u0,
		   CONST_B0,			// Assume constant strength background field, CONST_B0
		   CONST_m0,
		   CONST_n0,
		   CONST_a1nN,
		   CONST_b1nN,
		   CONST_g1eN,
		   CONST_h1eN,
		   CONST_a1iN,
		   CONST_b1iN,
		   CONST_g1iN,
		   CONST_h1iN,
		   CONST_a1sN,
		   CONST_b1sN;
    // Density spatial grid constants.
	float  CONST_gridSize;
	int	   CONST_min_x,
		   CONST_max_x,
		   CONST_min_y,
		   CONST_max_y,
   		   CONST_min_z,
		   CONST_max_z,
		   m_Nx,
		   m_Ny;
	// Surface discretization constants.
	int    CONST_surfN;
	double CONST_surfD;
	// Structure of plasma variables. 
	struct Plasma
	{
		double partRad;
	    double Tec;
		double Teh;
	    double Tp;
	    double Tw;
	    double Tnu;
	    double nec;
	    double neh;
	    double np;
	    double nw;
	    double vf;
	    double j0_ec;
	    double j0_eh;
	    double j0_p;
	    double j0_w;
	    double j0_nu;
	    double Xec;
	    double Xeh;
	    double Xp;
	    double Xw;
	    double Xnu;
	    double Mp;
	    double Mw;
	    double Mec;
	    double Meh;
	    double Q_TO_PHI;
	    double d_m;
	    double E_m;
	    double kTs;
	    double kappa;
	};
	Plasma m_plasma;

	// Pointer to function to evaluate derivatives of system. Chosen based on number of variables. 
	void (Solver::*Evaluate)(vector<double> &);
	// Pointer to function to evaluate Bfield.
	vector<double> (Solver::*Bfield)(const double &, const double &, const double &, 
		const vector<double> &, const double &);

	/*-------------------------------------------------------------------------------------*/
	/*-------------------------------------- Methods --------------------------------------*/
	/*-------------------------------------------------------------------------------------*/

	static bool Abs_compare(const double & a, const double & b);
	void   SaveData(vector<vector<double> > & data, const int & rows, const int & columns);
	void   Transform(float & Px, float & Py, float & Pz, const double & Ex,
			const double & Ey, const double & Ez);
	void   Cartesian2Geog(vector<double> & collisionLocation);
	long int GetDensityIndex(const double & x, const double & y, const double & z);
	double X2L(const double & x, const double & y, const double & z);
	double Declination(const double & x, const double & y, const double & z, const double & vx,
			const double & vy, const double & vz);

	// vector<double> Bfield_Connerney(const double & x, const double & y, const double & z, 
	// 		const vector<double> & moonPos, const double & r);
	// vector<double> Bfield_Simon(const double & x, const double & y, const double & z, 
	// 		const vector<double> & moonPos, const double & r);

	double PhiFromQ(const double & partQ);
	double SpecialF(const double & x, const int & n);
	double SpecialFB(double x, double B, const int & n);
	double phi2j(const double & phi, const double & v);
	void   EvaluateDeriv_NoCharge(vector<double> & Derivs);
	void   EvaluateDeriv_ConstCharge(vector<double> & Derivs);
	void   EvaluateDeriv(vector<double> & Derivs);
	void   ModMidpoint(const double & dt, vector<double> & y_n, const int & numSteps);
	int    BulirschStoerStep(const double & dt, double *y_n, const double & tol,
			bool & collideEnceladus);
};

#endif
