// A truncation of the power series of the hypergeometric function
// To be used only for qhat z < 1 and with a appropriately selected n higher than 
// the number of resonances in energy region of interest
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef TRUNCATED_HPP
#define TRUNCATED_HPP

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
        { set_Npars(2); };  

        // Evaluate the full term with angular dependence
        complex evaluate(double s, double zs)
        {
            complex as   = _alpha->evaluate(s);
            double q2hat = q2(s) / _lam2;

            complex sum = 0;
            for (int n = 0; n <= _nmax; n++) if ((n + isospin()) % 2 == 0) sum += pow(q2hat*zs, n) / (n - as);
            
            if (_removeZero) sum *= as;
            return (_constant) ? _g*sum - _gC : _g*sum;
        };

        // Allocate free parameters
        void allocate_parameters(std::vector<double> pars)
        {
            _lam2 = pars[0];
            _g    = pars[1];
            if (_constant){ _gC = pars[2]; }  
        };

        trajectory get_trajectory(){ return _alpha; };

        static const int kRemoveZeroPole  = 0;
        static const int kAddConstant     = 1;
        static const int kRemoveConstant  = 2;
        void set_option(int x)
        {
            switch (x)
            {
                case kRemoveZeroPole: { _removeZero = true; return; };
                case kAddConstant:    { _constant   = true; set_Npars(2+1); return; };
                case kRemoveConstant: 
                { 
                    if (!_constant) return; 
                    _constant = false; set_Npars(2); return;
                }
            }
            return;
        };

        protected:

        // Each isobar is given by a single trajectory only
        trajectory _alpha = nullptr;

        int _nmax       = 1;   // Max spin to include in sum
        double _g       = 1.;  // Coupling
        double _lam2    = 1.;  // Scale Lambda^2 indside \hat{q}_s^2
        
        bool _removeZero = false;

        // If to add a constant term in addition to the pole
        bool _constant   = false;
        double _gC = 1.;
    };
};

#endif