//  cgamma.cpp -- Complex gamma function.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns gamma function for complex argument 'z'.
//
//  Returns (1e308,0) if the real part of the argument is a negative integer
//  or 0 or exceeds 171.

#ifndef CGAMMA_HPP
#define CGAMMA_HPP

#include <complex>

namespace analyticRT
{
    // Gamma function of complex argument 'z'
    // optional parameter OPT = 1 can be used to return the log of Gamma(z)
    std::complex<double> cgamma(std::complex<double> z, int OPT = 0);
};

#endif 