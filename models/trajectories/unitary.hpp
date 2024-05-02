// Iterative and 
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef UNITARY_HPP
#define UNITARY_HPP

#include "amplitude.hpp"
#include "iterable.hpp"
#include "trajectory.hpp"

namespace analyticRT
{
    class unitary : public raw_iterable
    {
        public:

        unitary(int jmin, std::function<double(double)> F, std::string id)
        : raw_iterable(4*M2_PION, 5, F, id), _jmin(jmin)
        {};

        // RHC given by the logarithmic form 
        inline double RHC(double s)
        {
            if (s < _sRHC) return 0.;

            double  q2hat        = (s - _sRHC) / 4. / _lam2;
            double  rho          = sqrt(1. - _sRHC / s);
            double  gamma        = _gamma / PI;
            // double  delta_alpha  = previous_real(s) -     previous_real(_sRHC);
            double  delta_alpha  = _initial_guess(s) - _initial_guess(_sRHC);

            complex alphas       = previous_real(s) + I * previous_imag(s);

            double beta  = _g / (2.*_jmin + 1.);
            if (_pomeron) beta *= std::norm(1. + _gP/_g*alphas);

            if (s > 100)
            {
                if (delta_alpha < 0) { fatal("unitary::RHC","delta_alpha is negative!"); };
                return gamma*((1 + _jmin + delta_alpha)*log(q2hat) + log(_c*rho*beta/gamma));
            }; 
            return gamma*log(1. + rho/gamma*pow(q2hat, _jmin)*beta*(1. + _c*pow(q2hat, 1 + delta_alpha)));              
        };
        
        static const int kDefault    = 0;
        static const int kAddPomeron = 1;
        inline void set_option(int opt)
        {
            switch (opt)
            {
                case (kAddPomeron) : { _pomeron = true; set_Npars(6); return;}; 
                default: return;
            }
        };

        private:

        // Parameters are the scale and beta coefficients
        inline void allocate_parameters(std::vector<double> pars)
        {
            _lam2  = pars[0];              // Lambda^2 scale
            set_subtraction(pars[1]);      // alpha(0)
            _g     = pars[2];              // Residue 
            _gamma = pars[3];              // High-energy constant
            _c     = pars[4];
            if (_pomeron) _gP  = pars[5];
        };

        // Members related to the model for the imaginary part along the RHC
        int    _jmin  = 1;               // Lowest physical partial wave 
        double _lam2  = 3.;              // Scale of elastic unitarity

        // Free parameters
        double _g     = 1.; // Pole residue
        double _gamma = 1.; // Slope parameter
        double  _c    = 1.; // Inelastic parameter

        // Parameters related to including the constant contribution from a Pomeron
        bool _pomeron = false;
        double _gP = 0;
    };
};

#endif