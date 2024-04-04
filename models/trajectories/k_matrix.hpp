// Iteratable trajectory with a sqrt-log asymptotic.
// We allow an arbitrary n-th order polynomial in the imaginary part
// Near a pole this trajectory recovers the K-matrix formula
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef K_MATRIX_HPP
#define K_MATRIX_HPP

#include "trajectory.hpp"
#include "iterable.hpp"
#include <Math/Interpolator.h>

namespace analyticRT
{
    class k_matrix : public iterable
    {
        public: 

        k_matrix(double RHC, std::string id)
        : iterable(RHC, 4, id)
        {};

        // RHC given by the log with generic polynomial residue 
        inline double RHC(double s)
        {
            if (s < _sRHC) return 0.;

            double q2hat = (s - _sRHC) / 4. / _Lam2;        // Momenta divided by Regge scale
            double rho   = sqrt(1. - _sRHC / s) / (16.*PI); // Phase space

            double beta = 0; // Sum over polynomial of q^2
            for (int i = 0; i < _coeffs.size(); i++) beta += _coeffs[i] * pow(s - _sRHC, i);

            // For numerical stability, at asymptotic argumenets, simplify the equation
            if (s > 100)
            {
                return (_gamma / PI) *  (previousRePart(s) * log(q2hat) + log(beta));
            }

            // else use the Full expression
            return (_gamma / PI) *  log(1. + (PI / _gamma) * rho * beta * pow(q2hat, previousRePart(s)) );
        };

        // Parameters are the scale, subtraction, slope parameter, and arbitrarily many beta coefficients
        inline void allocate_parameters(std::vector<double> pars)
        {
            // Overall coupling
            _Lam2 = pars[0];
            set_subtraction(0., pars[1]);
            _gamma = pars[2]; 

            // Individual couplings
            _coeffs.clear();
            for (int i = 3; i < Npars(); i++) _coeffs.push_back(pars[i]);
        };

        inline std::vector<std::string> parameter_labels()
        {
            std::vector<std::string> labels = {"Lambda^2", "alpha(0)", "gamma"};
            for (int i = 3; i < Npars(); i++) labels.push_back( "c[" + std::to_string(i-2) + "]");
            return labels;
        };

        // The option is how many terms to consider in the polynomial
        inline void set_option(int n){ _coeffs.clear(); set_Npars(n + 3); };

        protected:

        // Members related to the model for the imaginary part along the RHC

        double _Lam2  = 5.;          // Scale of elastic unitarity
        double _gamma = 1.;          // Overall coupling
        std::vector<double> _coeffs; // (Real) coefficients of abitrary beta function inside the logarithm
    };
};

#endif