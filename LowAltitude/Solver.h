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

#ifndef SOLVERHEADERDEF
#define SOLVERHEADERDEF
#include <unordered_map>
#include <memory>
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
    void   CreateDistributionGrid(const float &min_alt,
            const float &max_alt, const int &num_alt,
            const float &max_vel, const int &num_vel,
            const float &max_rphi, const int &num_rphi, 
            const float &max_vphi, const int &num_vphi);

    vector<double> Bfield_Connerney(const double & x, const double & y, const double & z, 
            const vector<double> & moonPos, const double & r);
    vector<double> Bfield_Simon(const double & x, const double & y, const double & z, 
            const vector<double> & moonPos, const double & r);

    /* Function to simulate a particle and track where it has been in a spherical
       data grid about the moon. */
    // - Must be in header file because it is a templated function; template necessary
    //   to pass C-style array by reference
    //
    // TODO:
    //  - Check on angles for inclination; Think it will be, e.g., [180-15, 180] for
    // opening angle of 15. Grid data seems to be more designe for [0,15]. Need this
    // in CreateDistributionGrid, GetDistributionIndex, and HoverSim checking if in grid.
    //
    //
#if C_ARRAY
    template <size_t nvphi, size_t nphi, size_t nvr, size_t nr>
    void HoverSim(const double & dt, double *y,
        const double & weight, float (&flux)[nvphi][nphi][nvr][nr])
#else
    void HoverSim(const double & dt, double *y,
        const double & weight, std::unique_ptr<float[]> &flux,
        int nvphi, int nphi, int nvr, int nr)
#endif
    {
        if (CONST_max_rphi < 0) {
            std::cout << "Must call CreateDistributionGrid(...) first!\n";
        }

        float  partPos_x,
               partPos_y,
               partPos_z,
               partVel_x,
               partVel_y,
               partVel_z;
        double dt_next,
               curr_step;
        int    step_success;
        long int previous_ind = -1;

        // Loop To Run Simulation.
        bool within_grid = true;
        while (within_grid)
        {
            // Take one time step of size dt. If desired error has not been reached, reduce step
            // size and repeat until we have completed one step of size dt.
            dt_next   = dt;
            curr_step = 0.;
            while( abs(dt - curr_step) > 1e-10 ) {
                step_success = Solver::BulirschStoerStep(dt_next, y, 1e-12);
                if(step_success > 1e-10) {
                    curr_step += dt_next;
                    dt_next    = dt - curr_step;
                }
                else {
                    dt_next /= 2.;
                }
            }
            // Compute particle position and velocity with respect to Enceladus in 
            // spherical coordinate system centered at the south pole.
            partPos_x = y[0] - y[6];
            partPos_y = y[1] - y[7];
            partPos_z = y[2] - y[8];
            float r_p, phi_r, r_v, phi_v, temp;
            Solver::TransformSpherical(r_p, temp, phi_r, partPos_x,
                partPos_y,partPos_z,y[6],y[7],y[8], true);
            partVel_x = y[3] - y[9];
            partVel_y = y[4] - y[10];
            partVel_z = y[5] - y[11];
            Solver::TransformSpherical(r_v, temp, phi_v, partVel_x,
                partVel_y,partVel_z,y[6],y[7],y[8], false);

            // DEBUG TEST: transform w.r.t. moon position or velocity basis?
            // - Speed and inclination indifferent to this transformation,
            //   which probably means moon-velocity basis functions are the
            //   same as moon-position basis functions
            // partVel_x = y[3] - y[9];
            // partVel_y = y[4] - y[10];
            // partVel_z = y[5] - y[11];
            // float test_v, test_phi;
            // Solver::TransformSpherical(test_v, temp, test_phi, partVel_x,
            //     // partVel_y,partVel_z,y[6],y[7],y[8], false);
            //     partVel_y,partVel_z,y[9],y[10],y[11], false);
            // test_phi = PI - test_phi;
            // DEBUG TEST

            // Treat opening angle as away from negative vertical axis
            phi_r = PI - phi_r;
            phi_v = PI - phi_v;

            // Check if particle outside of data cone in terms of inclination angle
            // or altitude.
            if (phi_r > CONST_max_rphi || r_p > CONST_max_altitude) {
                within_grid = false;
                //if (phi_r > CONST_max_rphi)
                //    std::cout << "Exited grid, phi = " << phi_r*RAD2DEG << ", max angle = " << CONST_max_rphi*RAD2DEG << "\n";
                //else
                //    std::cout << "Exited grid, z = " << partPos_z << ", max altitude = " << CONST_max_altitude << "\n";
                //    std::cout << "vphi = " << phi_v*RAD2DEG << ", rphi = " << phi_r*RAD2DEG << ", "
                //        << "vel = " << r_v << "\n";
                break;
            }
            // Only track density for particles above minimum alttiude of grid
            else if (r_p < CONST_min_altitude) {
                continue;
            }
            else if (phi_v > CONST_max_vphi) {
                break;
            }

            // Add weighted residence time to space-velocity residence time profile
            std::vector<int> ind = Solver::GetDistributionIndex(r_p, phi_r, r_v, phi_v);   
            // DEBUG: make sure indices are not out of range
            if (ind[0] > nvphi) std::cout << "v_phi ind = " << ind[0] << ", max = " << nvphi << "\n";
            else if (ind[1] > nphi) std::cout << "r_phi ind = " << ind[1] << ", max = " << nphi << "\n";
            else if (ind[2] > nvr) std::cout << "|v| ind = " << ind[2] << ", max = " << nvr << "\n";
            else if (ind[3] > nr) std::cout << "|r| ind = " << ind[3] << ", max = " << nr <<
                ", max alt = " << CONST_max_altitude << ", alt = " << r_p << "\n";
            else {
            #if C_ARRAY
                flux[ind[0]][ind[1]][ind[2]][ind[3]] += weight;
            #else
                flux[get4dind(ind[0],ind[1],ind[2],ind[3],nvphi,nphi,nvr,nr)] += weight;
            #endif
            }
        }
        // Set CONST_numVariables and associated variables to inital status if equilibrium
        // charge was reached.
        if(m_change == 1) {
            CONST_numVariables += 1;
            CONST_partRad      *= 1e-6;
            CONST_partMass      = 4./3.*PI * pow(CONST_partRad,3) * 1000.;
            Evaluate            = &Solver::EvaluateDeriv;
            m_change            = 0;
        }       
    }

private:

    vector<vector<double> > m_currents;

    // Equation of motion constants
    double CONST_poleX,
           CONST_poleY,
           CONST_poleZ,
           CONST_planIntegrate;
    int    CONST_numVariables,
           m_change;

    // Conical (2d spherical) grid parameters for for hover simulations.
    // Stop trajectories when outside of cone defined by these parameters.
    float CONST_max_altitude;
    float CONST_min_altitude;
    int CONST_num_altitude;
    float m_d_altitude;
    float CONST_max_velocity;
    int CONST_num_velocity;
    float m_d_velocity;

    float CONST_max_rphi;
    int CONST_rphi;
    float m_rdphi;
    float CONST_max_vphi;
    int CONST_vphi;
    float m_vdphi;

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
           CONST_B0,            // Assume constant strength background field, CONST_B0
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
    int    CONST_min_x,
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
    void   Transform(float & Px, float & Py, float & Pz, const double & Mx,
            const double & My, const double & Mz);
    void TransformSpherical(float & Pd, float & Ptheta, float & Pphi,
            float & Px, float & Py, float & Pz, const double & Mx, const double & My,
            const double & Mz, bool shift_z);
    void   Cartesian2Geog(vector<double> & collisionLocation);
    long int GetDensityIndex(const double & x, const double & y, const double & z);
    std::vector<int> GetDistributionIndex(const float &alt, const float &rphi,
        const float &vel, const float &vphi);
    double X2L(const double & x, const double & y, const double & z);
    double Declination(const double & x, const double & y, const double & z, const double & vx,
            const double & vy, const double & vz);

    // vector<double> Bfield_Connerney(const double & x, const double & y, const double & z, 
    //      const vector<double> & moonPos, const double & r);
    // vector<double> Bfield_Simon(const double & x, const double & y, const double & z, 
    //      const vector<double> & moonPos, const double & r);

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
    int    BulirschStoerStep(const double & dt, double *y_n, const double & tol);
};

#endif
