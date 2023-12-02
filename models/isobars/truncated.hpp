// A truncation of the power series of the hypergeometric function
// To be used only for qhat z < 1 and with a appropriately selected n higher than 
// the number of resonances in energy region of interest
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobar.hpp"
#include "trajectory.hpp"
#include "hyp_2F1.hpp"

namespace analyticRT
{
    class truncated : public raw_isobar
    {
        public: 

        // Explicitly only allow a RHC
        truncated(key x, unsigned int isospin, int nmax, trajectory alpha, std::string id)
        : raw_isobar(x, isospin, id), _alpha(alpha), _nmax(nmax)
        {
            // NO free parameters
            set_Nfree(2 + alpha->Nfree());
        };

        // Evaluate the full term with angular dependence
        complex evaluate(double s, double zs)
        {
            complex as = _alpha->evaluate(s);
            double q2hat = q2(s) / _lambda2;

            complex sum = 0;
            for (int n = 0; n <= _nmax; n++)
            {
                if (( (n + isospin()) % 2) != 0) continue;
                sum += pow(q2hat*zs, n) / (n - as);
            };

            return _g * sum;
        };

        // Allocate free parameters
        void allocate_parameters(std::vector<double> pars)
        {
            _g       = pars[0];
            _lambda2 = pars[1];

            // Rest of the parameters get fed into alpha
            std::vector<double> alpha_pars(pars.begin() + 2, pars.end());  
            _alpha->set_parameters(alpha_pars);
        };

        protected:

        // Each isobar is given by a single trajectory only
        trajectory _alpha = nullptr;

        int _nmax       = 1;  // Max spin to include in sum
        double _g       = 1.; // Coupling
        double _lambda2 = 1.; // Scale Lambda^2

    };
};