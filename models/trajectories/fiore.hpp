// Dispersive trajectory resulting from the model in [1].
// Its technically iterative but not in the same way since each iteration just
// fits a constants
//
// So despite being iterative we dont inherit from iterable class
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:hep-ph/0011035
// ---------------------------------------------------------------------------

#ifndef FIORE_HPP
#define FIORE_HPP

#include "iterable.hpp"

namespace analyticRT
{
    class fiore : public raw_iterable
    {
        public: 

        // Explicitly only allow a RHC
        fiore(double R, std::string id)
        : raw_iterable(R, id)
        {
            set_Npars(4);
            initialize(); 
        };

        // Parameters are the intercept and coupling
        inline void allocate_parameters(std::vector<double> pars)
        {
            set_subtraction(0., pars[0]);
            for (int i = 1; i < Npars(); i++) _c[i-1] = pars[i];
        };

        inline std::vector<std::string> parameter_labels()
        {
            return {"alpha(0)", "c1", "c2", "cx"};
        };

        // Iteration procedure for lambdas
        inline void iterate()
        {
            for (int i = 0; i < 3; i++) _lam[i] = real_part(_s[i]);
        };
        
        protected:

        // RHC given by a simple square root 
        inline double RHC(double s)
        {
            double result = 0;
            for (int i = 0; i < _s.size(); i++)
            {
                if (s < _s[i]) continue;
                result += _c[i] * sqrt(s - _s[i]) * pow((s-_s[i])/s, _lam[i]);
            };

            return result;
        };

        // Additional thresholds above the lowest one
        std::array<double,3> _s = {_sRHC, 2.12, 30.}, _c, _lam;

        inline void initialize()
        {
            print("wawdfadsf");
            for (int i = 0; i < 3; i++) _lam[i] = 0.491 + 0.874*_s[i];
        };
    };
};

#endif