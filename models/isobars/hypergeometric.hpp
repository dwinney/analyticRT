// Isobar given by the hypergeometric function
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef HYPERGEOMETRIC_HPP
#define HYPERGEOMETRIC_HPP

#include "isobar.hpp"
#include "trajectory.hpp"
#include "hyp_2F1.hpp"

namespace analyticRT
{
    class hypergeometric : public raw_isobar
    {
        public: 

        hypergeometric(key x, unsigned int isospin, trajectory alpha, std::string id)
        : raw_isobar(x, isospin, id), _alpha(alpha)
        {
            set_Npars(2);
        };

        // Evaluate the full term with angular dependence
        complex evaluate(double s, double zs)
        {
            complex as = _alpha->evaluate(s);
            double q2hat = q2(s) / _lambda2;

            complex fp = hyp_2F1(1., -as, 1.-as, + q2hat*zs);
            complex fm = hyp_2F1(1., -as, 1.-as, - q2hat*zs);

            return - _g * (fp + pow(-1, isospin()) * fm) / as / 2.;
        };

        // Allocate free parameters
        void allocate_parameters(std::vector<double> pars)
        {
            _g       = pars[0];
            _lambda2 = pars[1];

            if (_option == kFixAlpha) return;
            
            // Rest of the parameters get fed into alpha
            std::vector<double> alpha_pars(pars.begin() + 2, pars.end());  
            _alpha->set_parameters(alpha_pars);
        };

        // Option whether or not to fit alpha parameters
        const static int kFixAlpha   = 0; // Default
        const static int kFloatAlpha = 1;
        void set_option(int x)
        {
            _option = x;
            switch (x)
            {
                case kFixAlpha:   { set_Npars(2);                   return; }
                case kFloatAlpha: { set_Npars(2 + _alpha->Npars()); return; }
                default: option_error();
            };
        };

        protected:

        // Each isobar is given by a single trajectory only
        trajectory _alpha = nullptr;

        double _g       = 1.; // Coupling
        double _lambda2 = 1.; // Scale Lambda^2

    };
};

#endif