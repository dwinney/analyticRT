// Trajectory with a sqrt-log asymptotic but which is exactly
// unitarizable at lowest PW
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef UNITARY_HPP
#define UNITARY_HPP

#include "iterative.hpp"

namespace analyticRT
{
    class unitary : public raw_trajectory
    {
        public:

        unitary(double R, int jmin, std::string id)
        : raw_trajectory(R, 4, id), _jmin(jmin)
        { initalize(); };

        protected:
        
        // Parameters are the scale and beta coefficients
        inline void allocate_parameters(std::vector<double> pars)
        {
            // Overall coupling
            set_subtraction(0., pars[0]);
            _Lam2  = pars[1];
            _g     = pars[2]; 
            _gamma = pars[3]; 

            for (int i = 0; i < _Niters; i++) iterate();
        };

        inline std::vector<std::string> parameter_labels(){ return {"alpha(0)", "Lambda^2", "g", "gamma"}; };

        protected:

        // RHC given by a simple square root 
        inline double RHC(double s)
        {
            if (s < _sRHC) return 0.;
            _s = s; // save so we can call it from anywhere 

            // For numerical stability, at asymptotic argumenets, simplify the equation
            if (s > 100) return gamma()*(RePart()*log(q2hat()) + log(beta()/gamma()));

            return gamma()* log(1. + rho()* beta()/gamma() * (1. + pow(q2hat(), RePart())));
        };

        // Save s at each evaluation step to not have to pass it around to each subfunction
        double _s;

        // Members related to the model for the imaginary part along the RHC
        int    _jmin  = 1;           // Lowest physical partial wave 
        double _Lam2  = 1.;          // Scale of elastic unitarity
        double _g = 1., _gamma = 1.; // Overall nearthreshold and high energy couplings
        double _delta = 0.1;

        // xi is the ratio of momenta squared over momenta evaluated at some scale s = Lambda^2
        inline double q2hat(){ return (_s - _sRHC) / 4. / _Lam2; };

        // Phase space factor
        inline double rho(){ return sqrt(1. - _sRHC / _s) / (16.*PI); }

        // Beta is the residue of the lowest physical partial wave
        inline double beta() { return _g / (2.*_jmin + 1.) * pow(q2hat(), _jmin); };
        inline double gamma(){ return _gamma / PI; };

        // Methods related to the interpolation of the real part 
        
        // For relatively small arguments, save an interpolation of the real part evaluated from DR
        int _Ninterp = 100; // Total interpolation will have 2*_Ninterp points
        ROOT::Math::Interpolator _ReAlphaInterp = ROOT::Math::Interpolator(2*_Ninterp, ROOT::Math::Interpolation::kCSPLINE);

        // Or at asymptotic arguments we match to a simple square-root
        double _sAsym = 200, _ReAlphaAsym;

        // Evaluate the last saved iteration of the real part
        inline double RePart(){ return (_s < _sAsym) ? _ReAlphaInterp.Eval(_s) : _ReAlphaAsym * sqrt(_s / _sAsym); };

        // Replace prevously saved interpolation by evaluating the DR
        inline void iterate()
        {
            std::vector<double> s, realpha;

            // Near threshold, (s < s1) we use a lot of points
            double s1 = 50;
            for (int i = 0; i < _Ninterp; i++)
            {
                double si  = _sRHC + (s1 - _sRHC) * double(i) / double(_Ninterp-1); 
                double rei = real_part(si);
                s.push_back(si); realpha.push_back(rei);
            };

            // Repeat for an other N points from s1 to the asymptotic matchpoint
            double s2 = s1 + 1;
            for (int i = 0; i < _Ninterp; i++)
            {
                double si  = s2 + (_sAsym - s2) * double(i) / double(_Ninterp-1); 
                double rei = real_part(si);
                s.push_back(si); realpha.push_back(rei);
            };
        
            _ReAlphaAsym = real_part(_sAsym);
            _ReAlphaInterp.SetData(s, realpha);
        };

        inline void initalize()
        {
             // The "zeroth" iteration, using the linear rho trajectory as an ansatz
            auto initial_guess = [] (double s)
            {
                return (0.5 + 0.9 * s) / sqrt(1. + s / 20); 
            };
            
            std::vector<double> s, realpha;

            for (int i = 0; i < 2*_Ninterp; i++)
            {
                double si  = _sRHC + (_sAsym - _sRHC) * double(i) / double(_Ninterp-1); 
                double rei = initial_guess(si);
                s.push_back(si); realpha.push_back(rei);
            };
            
            _ReAlphaAsym = initial_guess(_sAsym);
            _ReAlphaInterp.SetData(s, realpha);
        };
    };
};

#endif