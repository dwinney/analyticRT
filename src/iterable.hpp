// Abstract class for a trajectory model that can be iterated
// This means being able to save an interpolation of the previous real part
// with wihch to calculate the new imaginary part
//
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef ITERABLE_HPP
#define ITERABLE_HPP

#include "trajectory.hpp"
#include "Math/Interpolator.h"

namespace analyticRT
{
    class iterable : public raw_trajectory
    {
        public:

        iterable(double R_th, std::string id)
        : raw_trajectory(R_th, id)
        {};

        iterable(double R_th, int N, std::string id)
        : raw_trajectory(R_th, N, id)
        {};

        // This is the primary new function, which is to calculate the next iteration
        // We want to control exactly when we iterate in case we need to change parameters
        // in between iterations
        virtual void iterate();

        // Here is the initial guess of the real part
        inline void set_initial_pars(std::array<double,2> alpha, double mu2 = 20) 
        {
            _azInitial = alpha[0]; _apInitial = alpha[1]; _mu2 = mu2;
            this->initialize();
        };
        inline double initial_guess(double s){ return (_azInitial + _apInitial * s) / sqrt(1. + s / _mu2); };

        // Adjust the matching points and parameters involved in the interpolation
        inline void set_interp_pars(int N, std::array<double,2> pars)
        {
            _Ninterp = N; _s1 = pars[0]; _sAsym = pars[1];
        };

        // ---------------------------------------------------------------------------
    
        protected: 

        // Default parameters for the initial guess        
        double _azInitial = 0.5, _apInitial = 0.9, _mu2 = 20;

        // Save an interpolation of the real part evaluated from DR
        int _Ninterp = 500; // Total interpolation will have 2*_Ninterp points
        ROOT::Math::Interpolator _ReAlphaInterp = ROOT::Math::Interpolator(2*_Ninterp, ROOT::Math::Interpolation::kCSPLINE);

        // Or at asymptotic arguments we match to a simple square-root
        double _s1 = 100,_sAsym = 5000, _ReAlphaAsym;

        // Function to run all tasks needed to evaluate the first iteration
        // i.e. populate the interpolations with the inital_guess
        virtual void initialize();
    };
};

#endif