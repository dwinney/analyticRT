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

        unitary(int jmin, std::array<double,3> F, std::string id)
        : raw_iterable(4*M2_PION, 5, F, id), _jmin(jmin)
        {};

        // RHC given by the logarithmic form 
        inline double RHC(double s)
        {
            if (s < _sRHC) return 0.;

            double  q2hat        = (s - _sRHC) / 4. / _lam2;
            double  rho          = sqrt(1. - _sRHC / s);
            double  gamma        = _gamma / PI;
            double  delta        = previous_real(s) - previous_real(_sRHC);

            double sH = 4*_lam2 + _sRHC;
            complex alpha = (s <= sH)  ? previous_evaluate(s) : previous_evaluate(sH); 
            double beta  = _g / (2.*_jmin + 1.);
            double f = (_constant) ? std::norm(1. + _gp/_g*alpha) : 1.;

            double exponent = 1 + delta;

            if (s >= 200 && exponent > 0 && !is_zero(_c)) return gamma*((_jmin + exponent)*log(q2hat) + log(_c*rho*beta/gamma));
            return gamma*log(1. + rho/gamma*pow(q2hat, _jmin)*beta*(f + _c*pow(q2hat, exponent)));              
        };
        
        static const int kAddConstant    = 1;
        static const int kRemoveConstant = 2;
        inline void set_option(int opt)
        {
            switch (opt)
            {
                case (kAddConstant)    : { _constant = true;  set_Npars(6); return;}; 
                case (kRemoveConstant) : { _constant = false; set_Npars(5); return;}; 
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
            if (_constant) _gp  = pars[5];
        };

        // Members related to the model for the imaginary part along the RHC
        int    _jmin  = 1;               // Lowest physical partial wave 
        double _lam2  = 3.;              // Scale of elastic unitarity

        // Free parameters
        double _g     = 1.; // Pole residue
        double _gamma = 1.; // Slope parameter
        double  _c    = 1.; // Inelastic parameter

        // Parameters related to including the constant contribution from a higher trajectory
        bool _constant = false;
        double _gp = 0;
    };
};

#endif