
// Trajectory resulting from the model in [1]
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:hep-ph/0011035
// ---------------------------------------------------------------------------

#include "trajectory.hpp"

namespace analyticRT
{
    class fiore : public raw_trajectory
    {
        public: 

        // Explicitly only allow a RHC
        fiore(double R, std::string id)
        : raw_trajectory(R, id)
        {
            set_Npars(4);
        };

        // Parameters are the intercept and coupling
        inline void allocate_parameters(std::vector<double> pars)
        {
            set_subtraction(0., pars[0]);
            _c[0] = pars[1];
            _c[1] = pars[2];
            _c[2] = pars[3];

            iterate();
        };

        inline std::vector<std::string> parameter_labels()
        {
            return {"alpha(0)", "c1", "c2", "cx"};
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

        // Iteration procedure for lambdas
        int _Niter = 3;
        inline void iterate()
        {
            for (int i = 0; i < 3; i++)
            {
                _lam[i] = 0.491 + 0.874*_s[i];
            };

            for (int n = 0; n < _Niter; n++)
            {
                for (int i = 0; i < _Niter; i++)
                {
                    _lam[i] = real_part(_s[i]);
                };
            };
        };
    };
};