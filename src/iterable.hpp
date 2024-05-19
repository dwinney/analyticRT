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

        raw_iterable(double R_th, std::array<double,3> initpars, std::string id)
        : raw_trajectory(R_th, id), _azInitial(initpars[0]), _apInitial(initpars[1]), _mu2(initpars[2])
        { initialize(); };

        raw_iterable(double R_th, int N, std::array<double,3> initpars, std::string id)
        : raw_trajectory(R_th, N, id), _azInitial(initpars[0]), _apInitial(initpars[1]), _mu2(initpars[2])
        { initialize(); };

        // This is the primary new function, which is to calculate the next iteration
        // We want to control exactly when we iterate in case we need to change parameters
        // in between iterations
        void iterate();

        // Iterate multiple times changing the parameters in between each iteration
        template<int N>
        inline void iterate(std::array<std::vector<double>,N> vpars, int n = N)
        {
            for (int i = 0; i < n; i++)
            {
                set_parameters(vpars[i]);
                if (i != n-1) iterate();
            };
        };

        // Adjust the matching points and parameters involved in the interpolation
        inline void set_interp_pars(int N, std::array<double,2> pars)
        {
            _Ninterp = N; _s1 = pars[0]; _sAsym = pars[1];
            this->initialize();
        };
        
        // Evaluate the last saved iteration of the real part
        inline double previous_real(double s){ return (s <= _sAsym) ? _ReAlphaInterp.Eval(s) : _ReAlphaAsym * sqrt(s / _sAsym); };
        inline double previous_imag(double s){ return (s <= _sAsym) ? _ImAlphaInterp.Eval(s) : _ImAlphaAsym * sqrt(s / _sAsym)* log(s)/log(_sAsym); };
        
        // 
        inline complex previous_evaluate(double s){ return previous_real(s) + I*previous_imag(s); };

        // ---------------------------------------------------------------------------
    
        protected: 
        
        // Default parameters for the initial guess        
        double _azInitial = 0.5, _apInitial = 0.9, _mu2 = 20;
        inline double initial_guess(double s){ return (_azInitial + _apInitial*s)/sqrt(1+ s/_mu2); };

        // Save an interpolation of the real part evaluated from DR
        int _Ninterp = 100; // Total interpolation will have 2*_Ninterp points
        ROOT::Math::Interpolator _ReAlphaInterp = ROOT::Math::Interpolator(2*_Ninterp, ROOT::Math::Interpolation::kCSPLINE);
        ROOT::Math::Interpolator _ImAlphaInterp = ROOT::Math::Interpolator(2*_Ninterp, ROOT::Math::Interpolation::kCSPLINE);

        // Or at asymptotic arguments we match to a simple square-root
        double _s1 = 50,_sAsym = 200, _ReAlphaAsym, _ImAlphaAsym;

        // Function to run all tasks needed to evaluate the first iteration
        // i.e. populate the interpolations with the inital_guess
        virtual void initialize();
    };

    inline std::shared_ptr<raw_iterable> iterable(trajectory alpha)
    {
        std::shared_ptr<raw_iterable> ptr;
        ptr = std::dynamic_pointer_cast<raw_iterable>(alpha);
        return ptr;
    };
};

#endif