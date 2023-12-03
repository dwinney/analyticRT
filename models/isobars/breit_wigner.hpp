// Simple breit-wigner form
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef BREIT_WIGNER_HPP
#define BREIT_WIGNER_HPP

#include "isobar.hpp"
#include "legendre_P.hpp"

namespace analyticRT
{
    class breit_wigner : public raw_isobar
    {
        public: 

        // Explicitly only allow a RHC
        breit_wigner(key x, unsigned int isospin, int j, std::string id)
        : raw_isobar(x, isospin, id), _j(j)
        {
            set_Nfree(3);
        };

        // Evaluate the full term with angular dependence
        complex evaluate(double s, double zs)
        {
            return _g * legendre_P(_j, zs) / (s - _mass*_mass + I*csqrt(s - STH)*_width);
        };

        // Allocate free parameters
        void allocate_parameters(std::vector<double> pars)
        {
            _g     = pars[0];
            _mass  = pars[1];
            _width = pars[2];
        };

        protected:

        int _j = -1;

        double _g = 1., _mass = 0., _width = 0.;
    };
};

#endif