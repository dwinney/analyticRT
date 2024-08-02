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

#include "cgamma.hpp"

namespace analyticRT 
{
    // ------------------------------------------------
    std::complex<double> cgamma(std::complex<double> z,int OPT)
    // OPT = 0 for Gamma ; OPT = 1 for log(Gamma)
    {
    std::complex<double> I(0,1);
    std::complex<double> g, infini= 1e308+ 0.*I; // z0,z1
    double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
    double na,t,x1,y1,sr,si;
    int j,k;
    x1=9e9;
    na=9e9;

    static double a[] = {
        8.333333333333333e-02,
        -2.777777777777778e-03,
        7.936507936507937e-04,
        -5.952380952380952e-04,
        8.417508417508418e-04,
        -1.917526917526918e-03,
        6.410256410256410e-03,
        -2.955065359477124e-02,
        1.796443723688307e-01,
        -1.39243221690590};

    x = real(z); x1 = x;
    y = imag(z); y1 = y;
    if (x > 171) return infini;
    if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
        return infini;
    else if (x < 0.0) {
        x = -x;
        y = -y;
    }
    x0 = x;
    if (x <= 7.0) {
        na = (int)(7.0-x);
        x0 = x+na;
    }
    q1 = sqrt(x0*x0+y*y);
    th = atan(y/x0);
    gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(2.0*M_PI);
    gi = th*(x0-0.5)+y*log(q1)-y;
    for (k=0;k<10;k++){
        t = pow(q1,-1.0-2.0*k);
        gr += (a[k]*t*cos((2.0*k+1.0)*th));
        gi -= (a[k]*t*sin((2.0*k+1.0)*th));
    }
    if (x <= 7.0) {
        gr1 = 0.0;
        gi1 = 0.0;
        for (j=0;j<na;j++) {
        gr1 += (0.5*log((x+j)*(x+j)+y*y));
        gi1 += atan(y/(x+j));
        }
        gr -= gr1;
        gi -= gi1;
    }
    if (x1 <= 0.0) {
        q1 = sqrt(x*x+y*y);
        th1 = atan(y/x);
        sr = -sin(M_PI*x)*cosh(M_PI*y);
        si = -cos(M_PI*x)*sinh(M_PI*y);
        q2 = sqrt(sr*sr+si*si);
        th2 = atan(si/sr);
        if (sr < 0.0) th2 += M_PI;
        gr = log(M_PI/(q1*q2))-gr;
        gi = -th1-th2-gi;
        x = x1;
        y = y1;
    }
    if (OPT == 0) {
        g0 = exp(gr);
        gr = g0*cos(gi);
        gi = g0*sin(gi);
    }
    g = gr + I*gi;
    return g;
    }
};