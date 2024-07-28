// Iteratable trajectory with a sqrt-log asymptotic.
// We allow an arbitrary n-th order polynomial in the imaginary part
// Near a pole this trajectory recovers the K-matrix formula
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include "trajectory.hpp"
#include "iterable.hpp"
#include <Math/Interpolator.h>

namespace analyticRT
{
    class polynomial : public raw_iterable
    {
        public: 

        polynomial(double sth, int jmin, std::function<double(double)> F, std::string id)
        : raw_iterable(sth, 4, F, id)
        {};

        // RHC given by the log with generic polynomial residue 
        inline double RHC(double s)
        {
            if (s < _sRHC) return 0.;

            double  q2hat   = (s - _sRHC) / 4. / _lam2;
            double  rho     = sqrt(1. - _sRHC / s);
            double gamma    = _gamma/PI;

            double r = 0; // Sum over polynomial of q^2
            for (int i = 0; i < _coeffs.size(); i++) r += _coeffs[i] * pow(q2hat, i*_jmin);
            r += _coeff_last*pow(q2hat, previous_real(s));

            if (s > 200)
            {
                 return (previous_real(s) > 0) ? gamma*(previous_real(s)*log(q2hat) + log(rho/gamma*_coeff_last))
                                               : gamma*(_jmin+_coeffs.size())*log(q2hat);
            };

            // else use the Full expression
            return gamma*log(1. + rho/gamma*r);
        };

        // Parameters are the scale, subtraction, slope parameter, and arbitrarily many beta coefficients
        inline void allocate_parameters(std::vector<double> pars)
        {
            // Overall coupling
            _lam2 = pars[0];
            set_subtraction(0., pars[1]);
            _gamma = pars[2]; 

            // Individual couplings
            _coeffs.clear();
            for (int i = 3; i < Npars()-1; i++) _coeffs.push_back(pars[i]);
            _coeff_last = pars[Npars()-1];
        };

        inline std::vector<std::string> parameter_labels()
        {
            std::vector<std::string> labels = {"Lambda^2", "alpha(0)", "gamma"};
            for (int i = 3; i < Npars()-1; i++) labels.push_back( "c[" + std::to_string(i-2) + "]");
            labels.push_back("c[alpha]");
            return labels;
        };

        // The option is how many terms to consider in the polynomial
        inline void set_option(int n){ _coeffs.clear(); set_Npars(n + 4); };

        protected:

        // Members related to the model for the imaginary part along the RHC

        int _jmin = 1;
        double _lam2  = 5.;          // Scale of elastic unitarity
        double _gamma = 1.;          // Overall coupling
        double _coeff_last = 1.;
        std::vector<double> _coeffs; // (Real) coefficients of abitrary beta function inside the logarithm
    };
};

#endif