// Abstract class for a trajectory model that can be iterated
// This means being able to save an interpolation of the previous real part
// with wihch to calculate the new imaginary part
//
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "iterable.hpp"

namespace analyticRT
{
    void iterable::initialize()
    {
            // Save an interpolation of the real part using previous iteration
            std::vector<double> s, realpha;

            // Near threshold, (s < s1) we use a lot of points
            for (int i = 0; i < _Ninterp; i++)
            {
                double si  = _sRHC + (_s1 - _sRHC) * double(i) / double(_Ninterp-1); 
                double rei = initial_guess(si);
                s.push_back(si); realpha.push_back(rei);
            };

            // Repeat for an other N points from s1 to the asymptotic matchpoint
            double s2 = _s1 + 1;
            for (int i = 0; i < _Ninterp; i++)
            {
                double si  = s2 + (_sAsym - s2) * double(i) / double(_Ninterp-1); 
                double rei = initial_guess(si);
                s.push_back(si); realpha.push_back(rei);
            };

            // Load everything to the correct interpolators     
            _ReAlphaAsym = real_part(_sAsym);
            _ReAlphaInterp.SetData(s, realpha);
    };

    // Calculate the next iteration
    // Basically the same as initialize() except with the previous real part
    void iterable::iterate()
    {
        // Save an interpolation of the real part using previous iteration
        std::vector<double> s, realpha;
        // Near threshold, (s < s1) we use a lot of points
        double s1 = 100;
        for (int i = 0; i < _Ninterp; i++)
        {
            double si  = _sRHC + (s1 - _sRHC) * double(i) / double(_Ninterp-1); 
            double rei = real_part(si);
            s.push_back(si); realpha.push_back(rei);
        };

        // Repeat for an other N points from s1 to the asymptotic matchpoint
        double s2 = s1 + 1;
        for (int i = 0; i < _Ninterp; i++)
        {
            double si  = s2 + (_sAsym - s2) * double(i) / double(_Ninterp-1); 
            double rei = real_part(si);
            s.push_back(si); realpha.push_back(rei);
        };

        // Load everything to the correct interpolators     
        _ReAlphaAsym = real_part(_sAsym);
        _ReAlphaInterp.SetData(s, realpha);
    };
};