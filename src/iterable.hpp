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
    class raw_iterable : public raw_trajectory
    {
        public:

        raw_iterable(double R_th, std::function<double(double)> F, std::string id)
        : raw_trajectory(R_th, id), _initial_guess(F)
        { initialize(); };

        raw_iterable(double R_th, int N, std::function<double(double)> F, std::string id)
        : raw_trajectory(R_th, N, id), _initial_guess(F)
        { initialize(); };

        // This is the primary new function, which is to calculate the next iteration
        // We want to control exactly when we iterate in case we need to change parameters
        // in between iterations
        void iterate();

        // Adjust the matching points and parameters involved in the interpolation
        inline void set_interp_pars(int N, std::array<double,2> pars)
        {
            _Ninterp = N; _s1 = pars[0]; _sAsym = pars[1];
            this->initialize();
        };
        
        // Evaluate the last saved iteration of the real part
        inline double previousRePart(double s){ return (s < _sAsym) ? _ReAlphaInterp.Eval(s) : _ReAlphaAsym * sqrt(s / _sAsym); };
        
        // ---------------------------------------------------------------------------
    
        protected: 
        
        std::function<double(double)> _initial_guess;

        // Default parameters for the initial guess        
        double _azInitial = 0.5, _apInitial = 0.9, _mu2 = 20;

        // Save an interpolation of the real part evaluated from DR
        int _Ninterp = 100; // Total interpolation will have 2*_Ninterp points
        ROOT::Math::Interpolator _ReAlphaInterp = ROOT::Math::Interpolator(2*_Ninterp, ROOT::Math::Interpolation::kCSPLINE);

        // Or at asymptotic arguments we match to a simple square-root
        double _s1 = 50,_sAsym = 200, _ReAlphaAsym;

        // Function to run all tasks needed to evaluate the first iteration
        // i.e. populate the interpolations with the inital_guess
        void initialize();
    };

    inline std::shared_ptr<raw_iterable> iterable(trajectory alpha)
    {
        std::shared_ptr<raw_iterable> ptr;
        ptr = std::dynamic_pointer_cast<raw_iterable>(alpha);
        return ptr;
    };
};

#endif