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
            for (int n = 0; n <= _nmax; n++) 
            { 
                if ((n + isospin()) % 2 == 0)
                {
                    if (n == 0 && _adler) sum += 1./(adler_pole(s) - as);
                    else sum += pow(q2hat*zs, n) / (n - as);
                };
            };
            return _g * sum;
        };

        // Allocate free parameters
        void allocate_parameters(std::vector<double> pars)
        {
            _lam2 = pars[0];
            _g    = pars[1];
            if (_adler) _gA = pars[2];
        };

        trajectory get_trajectory(){ return _alpha; };

        static const int kAddAdlerZero    = 0;
        static const int kRemoveAdlerZero = 1;
        void set_option(int x)
        {
            switch (x)
            {
                case kAddAdlerZero:    { _adler = true; set_Npars(3); return; };
                case kRemoveAdlerZero: 
                {
                    if (!_adler) return; 
                    _adler = false; set_Npars(2); return;
                };
            }
            return;
        };

        protected:

        // Each isobar is given by a single trajectory only
        trajectory _alpha = nullptr;

        int _nmax       = 1;   // Max spin to include in sum
        double _g       = 1.;  // Coupling
        double _lam2    = 1.;  // Scale Lambda^2 indside \hat{q}_s^2

        // Related to Adler Zero
        inline double adler_pole(double s)
        { 
            double sth = _alpha->sth();
            return _gA / (s - _sA) * (s - sth) / (_sA - sth);
        };
        bool _adler = false;
        double _gA = 1., _sA = M2_PION/2.;
    };
};

#endif