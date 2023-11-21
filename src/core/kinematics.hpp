// Here we gather all the kinematic quantities for the pi pi scattering
// Since we have one fixed process we dont need to group this into a movable class
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef KINEMATICS_HPP
#define KINEMATICS_HPP

#include "constants.hpp"

namespace analyticRT
{
    const double M   = M_PION;
    const double M2  = M_PION*M_PION;
    const double STH = 4.*M2;

    // CM momentum
    inline complex p(double s){ return csqrt(Kallen(s, M2, M2))/ 2./ csqrt(s); };
    inline double p2(double s){ return Kallen(s, M2, M2)/ 4./ s; };

    // variables from s-channel quantities
    inline  double t(double s, double zs){ return -2.*p2(s)*(1. - zs); };
    inline  double u(double s, double zs){ return -2.*p2(s)*(1. + zs); };


    // Cross channel angles
    inline double z(double s, double t, double u){ return s*(t - u)/Kallen(s, M2, M2); };
    inline double z_t(double s, double zs){ return z(t(s, zs), s, u(s, zs)); };
    inline double z_u(double s, double zs){ return z(u(s, zs), t(s, zs), s); };
};

#endif