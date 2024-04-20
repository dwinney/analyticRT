// Simple S-wave K-matrix form
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef KMATRIX_ISO_HPP
#define KMATRIX_ISO_HPP

#include "isobar.hpp"
#include "trajectory.hpp"

namespace analyticRT
{
    class k_matrix_trajectory; 

    class k_matrix : public raw_isobar
    {
        public: 

        // Explicitly only allow a RHC
        k_matrix(key x, unsigned int isospin, std::string id)
        : raw_isobar(x, isospin, id)
        { set_Npars(3);};

        // Evaluate the full term with angular dependence
        inline complex evaluate(double s, double zs)
        {
            double Kinv = adler_term(s) - (_b + _c*s);
            return 1./(Kinv - G(s));
        };

        inline trajectory get_trajectory()
        {
            // Find s at which Re(1/A) = 0
            auto alpha_sig = [this](double s){ return std::abs(std::real(1/evaluate(s,0))); };
            double s_sig = find_minimum(alpha_sig, {4*M2_PION, 1});
            trajectory alpha = new_trajectory<k_matrix_trajectory>(id());
            alpha->set_parameters({_b - adler_term(s_sig), _c});
            return alpha;
        }

        // Allocate free parameters
        inline void allocate_parameters(std::vector<double> pars){ _a = pars[0]; _b = pars[1]; _c = pars[2]; };
        inline std::vector<std::string> parameter_labels(){ return {"a", "b", "c"}; };

        // Chew-Mandelstam phase-space , subtracted at zero
        static inline complex G(double s)
        {
            double m1 = M_PION, m2 = M_PION;
            complex rho, xi;
            complex result;

            rho    = csqrt(Kallen(s, m1*m1, m2*m2)) / s;
            xi     = 1 - (m1+m2)*(m1+m2)/s;
            result = (rho*log((xi + rho)/(xi - rho)) - xi*(m2-m1)/(m2+m1)*log(m2/m1))/PI;
            return -result;
        };

        inline double adler_term(double s){ return _a/(s-_sA)*(s-STH)/(_sA-STH); };

        protected:

        double _sA = M2_PION/2, _a = 0, _b = 0, _c = 0.;
    };

    class k_matrix_trajectory : public raw_trajectory
    {
        public:
        
        k_matrix_trajectory(double sth, std::string id)
        : raw_trajectory(sth, 2, id)
        {};

        inline void allocate_parameters(std::vector<double> pars){ _b = pars[0]; _c = pars[1]; };
        inline complex evaluate(double s)
        {
            return (_b + _c*s) + k_matrix::G(s);
        };

        protected:

        double _b = 0, _c = 0.;
    };
};

#endif