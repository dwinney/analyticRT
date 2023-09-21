// Unitarizable trajectory with a sqrt-log asymptotic
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "trajectory.hpp"

namespace analyticRT
{
    class sqrtlog : public raw_trajectory
    {
        public: 

        // Explicitly only allow a RHC
        sqrtlog(double R, std::string id)
        : raw_trajectory(R, id)
        {
            set_Npars(4);
        };

        // Parameters are the scale and beta coefficients
        inline void allocate_parameters(std::vector<double> pars)
        {
            // Overall coupling
            _gamma = pars[0]; 

            // Individual couplings
            _coeffs.clear();
            for (int i = 1; i < Npars(); i++)
            {
                _coeffs.push_back(pars[i]);
            }
        };

        inline std::vector<std::string> parameter_labels()
        {
            std::vector<std::string> labels = {"gamma"};

            for (int i = 1; i < Npars(); i++)
            {
                labels.push_back( "c[" + std::to_string(i) + "]");
            };

            return labels;
        };

        // The option is how many terms to consider in the polynomial
        inline void set_option(int n)
        {
            _coeffs.clear();
            set_Npars(n + 1);
        };

        protected:

        // RHC given by a simple square root 
        inline double RHC(double s)
        {
            if (s < _sRHC) return 0.;

            // Zeta is the ratio of momenta squared over some scale
            double zeta = (s - _sRHC) / (_Lam2 - _sRHC);

            // Phase space factors
            double rho = sqrt(1. - _sRHC / s);

            // Arbitrary rational function entering the coupling at low energies    
            double beta = 0;
            for (int i = 0; i < _coeffs.size(); i++)
            {
                beta += _coeffs[i] * pow(s - _sRHC, i);
            };

            double realpha = (0.459075 + 0.9*s) / sqrt(1. + s/70.);

            if (s > 300)
            {
                return (_gamma / PI) * realpha * log(zeta);
            }

            return (_gamma / PI) * log(1. + rho * beta * pow(zeta, realpha) );
        };

        // (Real) coefficients of abitrary beta function inside the logarithm
        std::vector<double> _coeffs; 

        // Scale entering the zeta factors 
        double _Lam2 = 10.;

        // Overall coupling
        double _gamma = 1.;
    };
};