// Unitarizable trajectory with a sqrt-log asymptotic
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "trajectory.hpp"
#include <Math/Interpolator.h>

namespace analyticRT
{
    class sqrtlog : public raw_trajectory
    {
        public: 

        // Explicitly only allow a RHC
        sqrtlog(double R, std::string id)
        : raw_trajectory(R, id)
        {
            set_Npars(4);
            initalize();
        };

        // Parameters are the scale and beta coefficients
        inline void allocate_parameters(std::vector<double> pars)
        {
            // Overall coupling
            _gamma = pars[0]; 

            // Individual couplings
            _coeffs.clear();
            for (int i = 1; i < Npars(); i++)
            {
                _coeffs.push_back(pars[i]);
            }

            for (int i = 0; i < _Niters; i++)
            {
                iterate();
            };
        };

        inline std::vector<std::string> parameter_labels()
        {
            std::vector<std::string> labels = {"gamma"};

            for (int i = 1; i < Npars(); i++)
            {
                labels.push_back( "c[" + std::to_string(i) + "]");
            };

            return labels;
        };

        // The option is how many terms to consider in the polynomial
        inline void set_option(int n)
        {
            _coeffs.clear();
            set_Npars(n + 1);
        };

        protected:

        // RHC given by a simple square root 
        inline double RHC(double s)
        {
            if (s < _sRHC) return 0.;
            _s = s; // save so we can call it from anywhere 

            // For numerical stability, at asymptotic argumenets, simplify the equation
            if (s > 100)
            {
                return (_gamma / PI) *  (previousRePart() * log(xi()) + log(beta()));
            }

            // else use the Full expression
            return (_gamma / PI) * rho() * log(1. + beta() * pow(xi(), previousRePart()) );
        };

        private:

        // Save s at each evaluation step to not have to pass it around to each subfunction
        double _s;

        // Members related to the model for the imaginary part along the RHC

        double _Lam2 = 10.;          // Scale of elastic unitarity
        double _gamma = 1.;          // Overall coupling
        std::vector<double> _coeffs; // (Real) coefficients of abitrary beta function inside the logarithm

        // xi is the ratio of momenta squared over momenta evaluated at some scale s = Lambda^2
        inline double xi(){ return (_s - _sRHC) / (_Lam2 - _sRHC); };

        // Phase space factor
        inline double rho(){ return sqrt(1. - _sRHC / _s); }

        // Arbitrary rational function entering the coupling at low energies    
        inline double beta()
        {
            double beta = 0;
            for (int i = 0; i < _coeffs.size(); i++)
            {
                beta += _coeffs[i] * pow(_s - _sRHC, i);
            };
            return beta;
        };

        // Methods related to the interpolation of the real part 
        
        // For relatively small arguments, save an interpolation of the real part evaluated from DR
        int _Ninterp = 100; // Total interpolation will have 2*_Ninterp points
        ROOT::Math::Interpolator _ReAlphaInterp = ROOT::Math::Interpolator(2*_Ninterp, ROOT::Math::Interpolation::kCSPLINE);

        // Or at asymptotic arguments we match to a simple square-root
        double _sAsym = 200, _ReAlphaAsym;

        // Evaluate the last saved iteration of the real part
        inline double previousRePart()
        {
            return (_s < _sAsym) ? _ReAlphaInterp.Eval(_s) : _ReAlphaAsym * sqrt(_s / _sAsym);
        };

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
                return (0.5 + 0.9 * s) / sqrt(1. + s / 40); 
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