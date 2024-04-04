// GKPY [1] parameterizations for pi pi phase-shifts and inelasticities
// 
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------
// REFERENCES 
// [1] - https://arxiv.org/abs/1102.2183
// ---------------------------------------------------------------------------

#ifndef PIPI_HPP
#define PIPI_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "data_set.hpp"

namespace analyticRT
{
    class pipi
    {
        public: 

        static double phase_shift( int iso, int j, double s);
        static double inelasticity(int iso, int j, double s);
        static inline complex partial_wave(int iso, int j, double s)
        {
            if (s <= 4.*M_PION*M_PION) return 0.;
            
            double delta, eta, k;
            complex amp;

            k     = elastic_mom(s, 4.*M_PION*M_PION);
            delta = phase_shift( iso, j, s);
            eta   = inelasticity(iso, j, s);

            amp  = (eta * exp(2.*I*delta) - 1.) / (2.*I);
            amp /= phase_space(s) * (16.*PI);
            return amp;
        };

        // Produce a data_set with the specified partial wave
        static data_set partial_wave(int iso, int j, int N, std::array<double,2> range);

        private:

        static inline double conformal(double s, double s0)
        {
            return (sqrt(s) - sqrt(s0 - s)) / (sqrt(s) + sqrt(s0 - s));
        };
        static inline double elastic_mom(double s, double sth)
        {
            return sqrt(s - sth) / 2.;
        };
    };
};

#endif