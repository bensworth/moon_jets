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

#ifndef JETHEADERDEF
#define JETHEADERDEF
class Solver;
#include <unordered_map>
#include "H5Cpp.h"
using namespace H5;
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
    void SetInitCond(const vector<double> & pos, const vector<double> & vel);
    void SetMaxAngle(const float & theta);
    void SetSimulationData(const int & jetId, const int & codeVersion, const int & bField, 
          const int & charging, const int & jetRef, const int & numCore);
    void SetGridData(const int & gridMin_x, const int & gridMax_x, const int & gridMin_y,
          const int & gridMax_y, const int & gridMin_z, const int & gridMax_z,
          const float & gridSize_dx);
    void HoverSimOMP(Solver & systemSolver, const int &numAzimuth, const int &partRad_ind,
    const float &partRad, const int &initVel_ind, const float &initVel);

private:

    // Information about jet to save in binary file. 
    int             m_jetId,            //
                    m_codeVersion,      //
                    m_bFieldModel,      //
                    m_jetRef,           //
                    m_numCore,          //
                    m_dataID,           //
                    m_charging,         //
                    m_numGrid_x, m_numGrid_y, m_numGrid_z,
                    m_minGrid_x, m_minGrid_y, m_minGrid_z;
    float           m_gridSize_dx;      //
    static const int m_nr;
    static const int m_nphi;
    static const int m_nvr;
    static const int m_nvphi;
    static const float m_min_altitude;
    static const float m_max_altitude;
    static const float m_max_velocity;
    static const float m_max_rphi;
    static const float m_max_vphi;
    const float m_dr;
    const double m_dphi;
    const float m_dvr;
    const double m_dvphi;

    vector< vector<double> > m_velVolume;
    vector< vector<double> > m_locVolume;

    // Constant parameters/variables.
    int             CONST_numVariables; //
    double          m_angDistConst;     // normalization constant for angular distribution
    float           m_maxAngle,         // maximum jet ejection angle (degrees)
                    m_maxAngleRad,      // maximum jet ejection angle (radians)
                    m_latitude,         // latitude of the jet's location on moon
                    m_longitude,        // longitude of the jet's location on moon
                    m_jetZenith,        // zenith angle of the jet's tilt
                    m_jetAzimuth,       // azimuth of the jet's tilt
                    m_v_gas,            // Plume gas speed for speed distribution (km/s)
                    m_critRad;          // Plume critical radius for speed distribution (um)
    vector<double>  m_jetPos;           // unit vector of jet's location in inertial frame
    vector<double>  m_jetDir;           // jet's ejection direction in inertial frame
    vector<double>  m_ex, m_ey, m_ez;   // Moon body-fixed orthonormal system
                                        // (+z: angular momentum, +x: moon apex)
    vector<double>  m_vx, m_vy, m_vz;   // Jet's orthonormal system
                                        // (+z: jet's ejection dir., +x: moon apex)
    vector<double>  m_moonPos,          // Moon position at start of simulation
                    m_moonVel,          // Moon velocity vector at start of simulation
                    m_parPos,           // Dust position at ejection
                    m_parNonEjVel;      //

    /*-------------------------------------------------------------------------------------*/
    /*-------------------------------------- Methods --------------------------------------*/
    /*-------------------------------------------------------------------------------------*/

    template<typename T> ostream& binary_write(ostream& stream, const T& value)
    {
        return stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    void   DensityWrite(const unordered_map<long int,pair<float,float> > & Density, const float & inclination,
           const int & numAzimuth, const int & partRadId, const int & monteCarlo = 0);
    vector<double> GetEjecVelocity(const float & azimuth, const float & inc);
    void   SetAngDistConst(const float & inc0, const float & inc1);
    double AngleDistribution(const double & inc0, const double & inc1, const double & ang);
    float  SampleIncAngle(const float & inc0, const float & inc1, const float & sample);
    void   UpdateDensity(unordered_map<long int,pair<float,float> > & local,
           unordered_map<long int,pair<float,float> > & aggregate);

    // TODO : check ordering; why is dataset.write() different than dataspace.write?
    template <size_t nr, size_t nphi, size_t nvr, size_t nvphi>
    void HDF5DistWrite(float (&residenceTime)[nr][nphi][nvr][nvphi],
        const int &numAzimuth, const int &partRad_ind, const float &partRad,
        const int &initVel_ind, const float &initVel)
    {
        // Try block to detect exceptions raised by any of the calls inside it
        try
        {
            // Turn off the auto-printing when failure occurs so that we can
            // handle the errors appropriately
            Exception::dontPrint();

            // Create file; H5F_ACC_TRUNC means if the file exists, open as
            // read only, otherwise create new file
            std::string  file_name = "./data/EncResidenceTime_r" + std::to_string(partRad_ind) +
                "_s" + std::to_string(initVel_ind) + ".hdf5";
            H5File file(file_name, H5F_ACC_TRUNC);
            std::string dataset_name;
            FloatType dfloat(PredType::NATIVE_FLOAT);  // Native C++ Float32
            IntType dint(PredType::NATIVE_INT);        // Native C++ Int

            // 4d array storage
            int ndims = 4;
            hsize_t dim4[4] = {nr, nphi, nvr, nvphi};  // dataset dimensions
            DataSpace dataspace4d(ndims, dim4);                  // Create data space w/ given dims
            dataset_name = "residence_time";
            DataSet dataset2d = file.createDataSet(dataset_name, dfloat, dataspace4d);
            dataset2d.write(residenceTime, dfloat);    // Write data to dataset (works w/ array a[][][])

            // Particle radius attribute
            dataset_name = "particle_radius"; 
            hsize_t dim1[1] = {1};   // dataset dimensions
            DataSpace radius_dspace(1, dim1);  // Create data space w/ given dims
            Attribute att_radius = file.createAttribute (dataset_name, dfloat, radius_dspace);
            att_radius.write(dfloat, &partRad);    // Write data to dataset

            // Particle speed attribute
            dataset_name = "initial_speed"; 
            DataSpace speed_dspace(1, dim1);  // Create data space w/ given dims
            Attribute att_speed = file.createAttribute (dataset_name, dfloat, speed_dspace);
            att_speed.write(dfloat, &initVel);    // Write data to dataset

            // Number azimuth attribute
            dataset_name = "number_azimuth"; 
            DataSpace azimuth_dspace(1, dim1);  // Create data space w/ given dims
            Attribute att_azimuth = file.createAttribute (dataset_name, dint, azimuth_dspace);
            att_azimuth.write(dint, &numAzimuth);    // Write data to dataset
        }
        // catch failures
        catch( FileIException error ) { error.printErrorStack(); } // H5File operations error
        catch( DataSetIException error ) { error.printErrorStack(); } // DataSet operations error
        catch( DataSpaceIException error ) { error.printErrorStack(); } // DataSpace operations errpr
        catch( DataTypeIException error ) { error.printErrorStack(); } // DataSpace operations error
    }
};

#endif