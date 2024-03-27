// Most implementations of Legendre functions are restricted to their physical domain
// here we by hand code up small cases of Pl for in general complex angles
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef LEGENDRE_P_HPP
#define LEGENDRE_P_HPP

#include "debug.hpp"

namespace analyticRT
{
    template<typename T> 
    T legendre_P(unsigned int j, T z)
    {
        switch (j) 
        {
            case 0: return 1.;
            case 1: return z;
            case 2: return 0.5*(3.*z*z - 1.);
            case 3: return 0.5*z*(5.*z*z - 3.);
            case 4: return (35.*z*z*z*z - 30.*z*z + 3.)/8.;
            case 5: return z*(63.*z*z*z*z - 70.*z*z + 15.)/8.;
            default: {};
        };

        return std::nan("");
    };
};

#endif