// Simple S-wave K-matrix form
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef KMATRIX_ISO_HPP
#define KMATRIX_ISO_HPP

#include "isobar.hpp"

namespace analyticRT
{
    class k_matrix : public raw_isobar
    {
        public: 

        // Explicitly only allow a RHC
        k_matrix(key x, unsigned int isospin, std::string id)
        : raw_isobar(x, isospin, id)
        {
            set_Npars(3);
        };

        // Evaluate the full term with angular dependence
        inline complex evaluate(double s, double zs)
        {
            double Kinv = _a/(s-_sA)*s/_sA + _b + _c*s;
            return 1./(-Kinv - G(s));
        };

        // Allocate free parameters
        inline void allocate_parameters(std::vector<double> pars)
        {
            _a  = pars[0];
            _b  = pars[1];
            _c  = pars[2];
        };

        inline std::vector<std::string> parameter_labels()
        {
            return {"a", "b", "c"};
        };

        // Chew-Mandelstam phase-space , subtracted at zero
        inline complex G(double s)
        {
            double m1 = M_PION, m2 = M_PION;
            complex rho, xi;
            complex result;

            rho    = csqrt(Kallen(s, m1*m1, m2*m2)) / s;
            xi     = 1 - (m1+m2)*(m1+m2)/s;
            result = (rho*log((xi + rho)/(xi - rho)) - xi*(m2-m1)/(m2+m1)*log(m2/m1))/PI;
            return -result;
        };

        protected:

        double _sA = M2_PION/2, _b = 1, _c = 0, _a = 0.;
    };
};

#endif