#ifndef JETHEADERDEF
#define JETHEADERDEF
class Solver;
#include <unordered_map>
using namespace std;

class Jet
{
public:

	/*-------------------------------------------------------------------------------------*/
	/*-------------------------------------- Methods --------------------------------------*/
	/*-------------------------------------------------------------------------------------*/

	Jet();
	~Jet();

	void SetNumVariables(const int & n);
	void SetLocation(const float & latitude, const float & longitude, const float & jetZen = 0,
		  const float & jetAz = 0);
	void SetSpeedDist(const float & v_gas, const float & critRad);
	void SetInitCond(const vector<double> & pos, const vector<double> & vel);
	void SetMaxAngle(const float & theta);
	void SetSimulationData(const int & jetId, const int & codeVersion, const int & bField, 
		  const int & charging, const int & jetRef, const int & numCore);
	void SetGridData(const int & gridMin_x, const int & gridMax_x, const int & gridMin_y,
		  const int & gridMax_y, const int & gridMin_z, const int & gridMax_z,
		  const float & gridSize_dx);
	void SpecSimOMP(Solver & systemSolver, const vector<float> & partSpeeds,
		  const vector<double> & weights, const vector<float> & angleData, const float & totalTime,
		  const float & volume, const float & partRad, const int & partRadId);
	void MonteCarlo_Jet(Solver & systemSolver, const int & numSpeeds, const int & numAzimuth,
		  const float & inclination, const float & d_inclination, const float & totalTime,
		  const float & volume, const float & partRad, const int & partRadId);

	void CollisionMap(Solver & systemSolver, const vector<float> & partSpeeds,
		const vector<double> & weights, const vector<float> & angleData, const float & totalTime,
		const float & volume, const float & partRad, const int & partRadId);

private:

	// Information about jet to save in binary file. 
	int 			m_jetId,			//
					m_codeVersion, 		//
					m_bFieldModel,	 	//
					m_jetRef, 			//
					m_numCore, 			//
					m_dataID,			//
					m_charging, 		//
					m_numGrid_x, m_numGrid_y, m_numGrid_z,
					m_minGrid_x, m_minGrid_y, m_minGrid_z;
	float      	   	m_gridSize_dx;     	//
	// Constant parameters/variables.
    int 		   	CONST_numVariables; //
    double 			m_angDistConst;		// normalization constant for angular distribution
	float         	m_maxAngle,			// maximum jet ejection angle (degrees)
					m_maxAngleRad,		// maximum jet ejection angle (radians)
					m_latitude,     	// latitude of the jet's location on moon
					m_longitude,    	// longitude of the jet's location on moon
					m_jetZenith,    	// zenith angle of the jet's tilt
					m_jetAzimuth,   	// azimuth of the jet's tilt
					m_v_gas,			// Plume gas speed for speed distribution
					m_critRad;			// Plume critical radius for speed distribution
	vector<double> 	m_jetPos;			// unit vector of jet's location in inertial frame
	vector<double> 	m_jetDir;			// jet's ejection direction in inertial frame
    vector<double> 	m_ex, m_ey, m_ez;	// Moon body-fixed orthonormal system
    									// (+z: angular momentum, +x: moon apex)
    vector<double> 	m_vx, m_vy, m_vz;	// Jet's orthonormal system
    									// (+z: jet's ejection dir., +x: moon apex)
	vector<double> 	m_moonPos,			// Moon position at start of simulation
					m_moonVel, 			// Moon velocity vector at start of simulation
					m_parPos, 			// Dust position at ejection
					m_parNonEjVel;		//

	/*-------------------------------------------------------------------------------------*/
	/*-------------------------------------- Methods --------------------------------------*/
	/*-------------------------------------------------------------------------------------*/

	template<typename T> ostream& binary_write(ostream& stream, const T& value)
	{
	    return stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
	}

	void   DensityWrite(const unordered_map<long int,pair<float,float> > & Density, const float & inclination,
		   const int & numAzimuth, const int & partRadId, const int & monteCarlo = 0);
	void   CollisionWrite(const vector<vector<float> > Collisions, const float & inclination,
		   const int & particlesLaunched, const int & particlesCollided, const int & partRadId,
		   const int & monteCarlo = 0);
	vector<double> GetEjecVelocity(const float & azimuth, const float & inc);
	double  SpeedDistribution(const double & v, const float & partRad);
	float  SampleInitSpeed(const float & partRad, const double & shift);	
	void   SetAngDistConst(const float & inc0, const float & inc1);
	double AngleDistribution(const double & inc0, const double & inc1, const double & ang);
	float  SampleIncAngle(const float & inc0, const float & inc1, const float & sample);
	void   UpdateDensity(unordered_map<long int,pair<float,float> > & local,
		   unordered_map<long int,pair<float,float> > & aggregate);

};

#endif
