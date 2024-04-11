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

        unitary(double R, int jmin, std::string id)
        : raw_iterable(R, 5, id), _jmin(jmin)
        {};

        protected:
        
        // Parameters are the scale and beta coefficients
        inline void allocate_parameters(std::vector<double> pars)
        {
            set_subtraction(0., pars[0]); // alpha(0)
            _Lam2  = pars[1];             // Lambda^2 scale
            _g     = pars[2];             // Residue 
            _gamma = pars[3];             // High-energy constant
            _c     = pars[4];
        };

        inline std::vector<std::string> parameter_labels(){ return {"alpha(0)", "Lambda^2", "g", "gamma", "c"}; };

        // RHC given by the logarithmic form 
        inline double RHC(double s)
        {
            if (s < _sRHC) return 0.;

            double q2hat = (s - _sRHC) / 4. / _Lam2;
            double beta  = _g / (2.*_jmin + 1.);
            double rho   = sqrt(1. - _sRHC / s);
            double gamma = _gamma / PI;
            double delRe = previousRePart(s) - previousRePart(_sRHC);

            // For numerical stability, at asymptotic argumenets, simplify the equation
            if (s > 100) return gamma*((1 + _jmin + delRe)*log(q2hat) + log(_c*rho*beta/gamma));
            return gamma*log(1. + rho*beta/gamma*pow(q2hat, _jmin)*(1. + _c*pow(q2hat, 1 + delRe)));              
        };

        protected:

        // Members related to the model for the imaginary part along the RHC
        int    _jmin  = 1;               // Lowest physical partial wave 
        double _Lam2  = 3.;              // Scale of elastic unitarity

        // Free parameters
        double _g     = 1.; // Pole residue
        double _gamma = 1.; // Slope parameter
        double  _c    = 1.; // Inelastic parameter
    };
};

#endif