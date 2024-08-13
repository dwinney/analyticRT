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

        unitary(int jmin, std::string id)
        : raw_iterable(4*M2_PION, 5, [](double s){ return 0.;}, id), _jmin(jmin)
        {};

        // RHC given by the logarithmic form 
        inline double RHC(double s)
        {
            if (s < _sRHC) return 0.;

            double  q2hat   = (s - _sRHC) / 4. / _lam2;
            double  rho     = sqrt(1. - _sRHC / s);
            double  gamma   = _gamma / PI;
            
            double beta = pow(q2hat, _jmin)*_g/(2.*_jmin+1.);
            if (option() == kAddConstant) beta *= std::norm(1. + _gp/_g*previous_evaluate(s));
            if (option() == kExpandAlpha) beta *= std::norm(1. + _gp/_g*(_alphaSUB + _c*pow(q2hat, 4./3)));

            // Whether or not to have 
            bool   second_term  = option() != kExpandAlpha;
            double exponent     = (second_term) ? 1. + previous_real(s) : 0.;
            bool   to_simplify  = (second_term) ? false                 // We have driving term
                                : (s >= 200)                            // q2hat is large
                               && (exponent > 0)                        // exponent is positive
                               && (_c*pow(q2hat, exponent) >= 10*beta); // _c is not too small

            if (to_simplify) return gamma*(exponent*log(q2hat) + log(_c/gamma));
            return gamma*log(1. + rho/gamma*(beta + second_term*_c*pow(q2hat, exponent)));         
        };
        
        static const int kDefault        = 0;
        static const int kAddConstant    = 1;
        static const int kExpandAlpha    = 2;
        inline void set_option(int opt)
        {
            switch (opt)
            {
                case (kDefault)        : { _gp = 0; set_Npars(5); break; };
                case (kAddConstant)    : { set_Npars(6); break; }; 
                case (kExpandAlpha)    : { set_Npars(6); break; };
                default: return;
            }
            _option = opt;
        };

        private:

        // Parameters are the scale and beta coefficients
        inline void allocate_parameters(std::vector<double> pars)
        {
            double sub_point = (option() == kExpandAlpha) ? _sRHC : 0.;
            _lam2  = pars[0];                    // Lambda^2 scale
            set_subtraction(sub_point, pars[1]); // alpha(s_sub)
            _g     = pars[2];                    // Coupling 
            _gamma = pars[3];                    // Slope parameter
            _c     = pars[4];                    // Extra coupling in polynomial
            if (option() != kDefault) _gp  = pars[5];
        };

        // Members related to the model for the imaginary part along the RHC
        int    _jmin  = 1;               // Lowest physical partial wave 
        double _lam2  = 3.;              // Scale of elastic unitarity

        // Free parameters
        double _g     = 1.; // Pole residue
        double _gamma = 1.; // Slope parameter
        double _c     = 0.; // Polynomial coefficient

        // Parameters related to including the constant contribution from a higher trajectory
        double _gp = 0;
    };
};

#endif