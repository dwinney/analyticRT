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

    const double C[3][3] = 
    {
        {1./3.,  1.,     5./3.},
        {1./3.,  1./2., -5./6.},
        {1./3., -1./2.,  1./6.}
    };

    // CM momentum
    inline double q(double s){  return (s >= 4*M2*M2) ? sqrt(Kallen(s, M2, M2))/ 2./ sqrt(s) : NaN<double>(); };
    inline double q2(double s){ return Kallen(s, M2, M2)/ 4./ s; };

    // Chew-Mandelstam Phase space factor 
    inline complex phase_space_CM(double s)
    {
        double m1 = M_PION, m2 = M_PION;
        complex rho, xi;
        complex result;

        rho    = csqrt(Kallen(s, m1*m1, m2*m2)) / s;
        xi     = 1 - (m1+m2)*(m1+m2)/s;
        result = (rho*log((xi + rho) / (xi - rho)) - xi*(m2-m1)/(m2+m1)*log(m2/m1)) / PI;
        return result / (16*PI);
    };

    // Regular phhase space factor above cut
    inline double phase_space(double s){ return (s >= STH) ? sqrt(1. - STH / s) / (16*PI) : 0.; };

    // variables from s-channel quantities
    inline  double t_man(double s, double zs){ return -2.*q2(s)*(1. - zs); };
    inline  double u_man(double s, double zs){ return -2.*q2(s)*(1. + zs); };


    // Cross channel angles
    inline double z(double s, double t, double u){ return (t - u)/ (s - STH); };
    inline double z_t(double s, double zs){ return z(t_man(s, zs), s, u_man(s, zs)); };
    inline double z_u(double s, double zs){ return z(u_man(s, zs), t_man(s, zs), s); };
};

#endif