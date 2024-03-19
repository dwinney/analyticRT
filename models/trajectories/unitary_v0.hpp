// Trajectory with a sqrt-log asymptotic but which is exactly
// unitarizable at lowest PW.
// 
// -- Imaginary part goes like q2^jmin leading momentum behavior at any s
// -- hard coded for the rho P-wave
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef UNITARY_V0_HPP
#define UNITARY_V0_HPP

#include "amplitude.hpp"
#include "trajectory.hpp"

namespace analyticRT
{
    class unitary_v0 : public raw_trajectory
    {
        public:

        unitary_v0(double R, int jmin, std::string id)
        : raw_trajectory(R, 4, id), _jmin(jmin)
        { initialize(); };

        unitary_v0(double R, int jmin, amplitude amp, std::string id)
        : raw_trajectory(R, 4, id), _jmin(jmin), _amp(amp)
        { initialize(); };

        protected:
        
        // Parameters are the scale and beta coefficients
        inline void allocate_parameters(std::vector<double> pars)
        {
            // Overall coupling
            set_subtraction(0., pars[0]);
            _Lam2  = pars[1]; _sLam2 = (4*_Lam2 + STH)/2;
            _g     = pars[2]; 
            _gamma = pars[3]; 
            
            for (int i = 0; i < _Niters; i++) iterate();
        };

        inline std::vector<std::string> parameter_labels(){ return {"alpha(0)", "Lambda^2", "g", "gamma"}; };

        protected:

        // For unitarization when including cross-channels we need a pointer to the full
        // amplitude
        amplitude _amp = nullptr;

        // RHC given by a simple square root 
        inline double RHC(double s)
        {
            if (s < _sRHC + EPS) return 0.;
            _s = s; // save so we can call it from anywhere 

            // For numerical stability, at asymptotic argumenets, simplify the equation
            if (s > 100) return gamma()*(RePart()*log(q2hat()) + log(beta()/gamma()));
            
            return gamma()* log(1. + rho()* beta()/gamma() * (1. + pow(q2hat(), RePart())));
        };

        // Save s at each evaluation step to not have to pass it around to each subfunction
        double _s;

        // Members related to the model for the imaginary part along the RHC
        int    _jmin  = 1;               // Lowest physical partial wave 
        double _Lam2  = 3., _sLam2 = 5.; // Scale of elastic unitarity
        double _g     = 1., _gamma = 1.; // Overall nearthreshold and high energy couplings

        // xi is the ratio of momenta squared over momenta evaluated at some scale s = Lambda^2
        inline double q2hat(){ return (_s - _sRHC) / 4. / _Lam2; };

        // Phase space factor
        inline double rho(){ return sqrt(1. - _sRHC / _s) / (16.*PI); }

        // Beta is the residue of the lowest physical partial wave
        inline double gamma(){ return _gamma / PI; };
        inline double beta() 
        {
            double r = _g / (2.*_jmin + 1.) * pow(q2hat(), _jmin);
            if (_amp == nullptr) return r;

            complex alpha = RePart() - I*ImPart();
            return r * std::norm(1. + (_jmin - alpha)/r*Ftilde()); 
        };

        //--------------------------------------------------------------------------------------------------------
        // Methods related to the interpolation of the real part 
        
        // For relatively small arguments, save an interpolation of the real part evaluated from DR
        int _Ninterp = 100; // Total interpolation will have 2*_Ninterp points
        ROOT::Math::Interpolator _ReAlphaInterp = ROOT::Math::Interpolator(2*_Ninterp, ROOT::Math::Interpolation::kCSPLINE);
        ROOT::Math::Interpolator _ImAlphaInterp = ROOT::Math::Interpolator(  _Ninterp, ROOT::Math::Interpolation::kCSPLINE);

        // Or at asymptotic arguments we match to a simple square-root
        double _sAsym = 200, _ReAlphaAsym, _ImAlphaAsym;

        // Evaluate the last saved iteration of the real part
        inline double RePart(){ return (_s < _sAsym) ? _ReAlphaInterp.Eval(_s) : _ReAlphaAsym * sqrt(_s / _sAsym); };
        inline double ImPart(){ return (_s < _sLam2) ? _ImAlphaInterp.Eval(_s) : _ImAlphaAsym; };

        //--------------------------------------------------------------------------------------------------------
        // Methods related to the interpolation of cross channel contributions
        double _fTildeAsym = 0., _sMatch = 0;
        ROOT::Math::Interpolator _ftildeInterp = ROOT::Math::Interpolator(_Ninterp, ROOT::Math::Interpolation::kCSPLINE);

        inline double Ftilde(){ return (_s < _sMatch) ? _ftildeInterp.Eval(_s) : _fTildeAsym; };

        //--------------------------------------------------------------------------------------------------------
        // Replace prevously saved interpolation by evaluating the DR
        inline void iterate()
        {
            // Save an interpolation of the real part using previous iteration
            std::vector<double> s, realpha;
            // Near threshold, (s < s1) we use a lot of points
            double s1 = 50;
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

            if (_amp == nullptr) return;

            // Also save an interpolation of the cross-channel components
            // from sRHC to the scale lambda2
            std::vector<double> ss, ftilde, imalpha;
            for (int i = 0; i < _Ninterp; i++)
            {
                double si  = _sRHC + (_sLam2 - _sRHC) * double(i) / double(_Ninterp-1); 
                double fti =  (si < _sLam2) ? std::real(_amp->cross_projection(_jmin, _jmin, si))
                                            : std::real(_amp->cross_projection(_jmin, _jmin, _sLam2) );
                double imi = RHC(si);
                ss.push_back(si); ftilde.push_back(fti); imalpha.push_back(imi);
            }

            // Load everything to the correct interpolators     
            _ReAlphaAsym = real_part(_sAsym);
            _ReAlphaInterp.SetData(s, realpha);
            _ImAlphaAsym = imalpha.back(); _fTildeAsym  = ftilde.back();

            _sMatch = ss.back() + 0.1;
            ss.push_back(_sMatch); 
            imalpha.push_back(_ImAlphaAsym);     ftilde.push_back(_fTildeAsym);
            _ImAlphaInterp.SetData(ss, imalpha); _ftildeInterp.SetData(ss, ftilde);
        };

        inline void initialize()
        {
            // Load initial interpolation of the Re alpha(s)
            
            std::vector<double> s, realpha, zeros;

            for (int i = 0; i < 2*_Ninterp; i++)
            {
                double si  = _sRHC + (_sAsym - _sRHC) * double(i) / double(2.*_Ninterp-1); 
                double rei = initial_guess(si);
                s.push_back(si); realpha.push_back(rei);
            };
            
            _ReAlphaAsym = initial_guess(_sAsym);
            _ReAlphaInterp.SetData(s, realpha);

            if (_amp == nullptr) return;

            // The "zero-th" iteration assumes zero cross channel contribution
            s.clear();
            for (int i = 0; i < _Ninterp; i++)
            {
                double si  = _sRHC + (_sLam2 - _sRHC) * double(i) / double(_Ninterp-1); 
                s.push_back(si); zeros.push_back(0.);
            };
            
            _ImAlphaInterp.SetData(s, zeros); _ftildeInterp.SetData( s, zeros);
        };
    };
};

#endif