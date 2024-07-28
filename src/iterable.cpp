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
    void raw_iterable::initialize()
    {
        // Save an interpolation of the real part using previous iteration
        std::vector<double> s, realpha, imalpha;

        // Near threshold, (s < s1) we use a lot of points
        for (int i = 0; i < _Ninterp; i++)
        {
            double si  = _sRHC + (_s1 - _sRHC) * double(i) / double(_Ninterp-1); 
            double rei = _initial_guess(si);
            s.push_back(si);
            realpha.push_back(rei);
            imalpha.push_back(0.);
        };

        // Repeat for an other N points from s1 to the asymptotic matchpoint
        double s2 = _s1    + 1;
        double s3 = _sAsym;
        for (int i = 0; i < _Ninterp; i++)
        {
            double si  = s2 + (s3 - s2) * double(i) / double(_Ninterp-1); 
            double rei = _initial_guess(si);
            s.push_back(si); 
            realpha.push_back(rei);
            imalpha.push_back(0.);
        };

        // Load everything to the correct interpolators     
        _ReAlphaInterp.SetData(s, realpha);
        _ReAlphaAsym = _ReAlphaInterp.Eval(_sAsym);

        _ImAlphaInterp.SetData(s, imalpha);
        _ImAlphaAsym = _ImAlphaInterp.Eval(_sAsym);
    };

    // Calculate the next iteration
    // Basically the same as initialize() except with the previous real part
    void raw_iterable::iterate()
    {
        // Save an interpolation of the real part using previous iteration
        std::vector<double> s, realpha, imalpha;
        
        // Near threshold, (s < s1) we use a lot of points
        for (int i = 0; i < _Ninterp; i++)
        {
            double si  = _sRHC + (_s1 - _sRHC) * double(i) / double(_Ninterp-1); 
            complex alphai = evaluate(si);
            s.push_back(si); 
            realpha.push_back(std::real(alphai));
            imalpha.push_back(std::imag(alphai));
        };

        // Repeat for an other N points from s1 to the asymptotic matchpoint
        double s2 = _s1    + 1;
        double s3 = _sAsym;
        for (int i = 0; i < _Ninterp; i++)
        {
            double si  = s2 + (s3 - s2) * double(i) / double(_Ninterp-1); 
            complex alphai = evaluate(si);

            s.push_back(si); 
            realpha.push_back(std::real(alphai));
            imalpha.push_back(std::imag(alphai));
        };

        _ReAlphaInterp.SetData(s, realpha);
        _ReAlphaAsym = _ReAlphaInterp.Eval(_sAsym);

        _ImAlphaInterp.SetData(s, imalpha);
        _ImAlphaAsym = _ImAlphaInterp.Eval(_sAsym);
    };
};