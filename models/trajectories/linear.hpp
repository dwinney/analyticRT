// Non-dispersive trajectory that assumes a simple linear form and sqrt imag part
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef LINEAR_HPP
#define LINEAR_HPP

#include "trajectory.hpp"
#include "kinematics.hpp"

namespace analyticRT
{
    class linear : public raw_trajectory
    {
        public: 

        // Explicitly only allow a RHC
        linear(double Rth, std::string id)
        : raw_trajectory(Rth, id)
        {
            set_Npars(3);
        };

        // Instead of defining the discontinuity and letting the dispersion relation
        // calculate everything, override real and imaginary parts directly

        inline complex evaluate(double s)
        {
            return _a0 + _aP*s - _gam*csqrt(_sRHC - s);
        };

        // Parameters are the intercept and coupling
        inline void allocate_parameters(std::vector<double> pars)
        {
            _a0 = pars[0]; _aP = pars[1], _gam = pars[2];
        };

        inline std::vector<std::string> parameter_labels()
        {
            return {"alpha(0)", "alpha'", "gamma"};
        };

        protected:

        double _a0 = 1., _aP = 0., _gam = 0.;

    };
};

#endif