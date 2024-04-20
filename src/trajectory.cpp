// Abstract class for a general Regge trajectory defined by an once-subtracted
// dispersion relation across RHC and optionally also a LHC
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "trajectory.hpp"

namespace analyticRT
{
    // -----------------------------------------------------------------------
    // Options and stuff

    // Intermediate function to check the input vector is right size
    void raw_trajectory::set_parameters(std::vector<double> pars)
    {
        if (pars.size() != _Npars) return error(id() + "::set_parameters() :"
                                                     + " Input vector do not match expected number of free parameters!");
        this->allocate_parameters(pars);
    };

    // Default parameter labels are just par[0], par[1], ...
    std::vector<std::string> raw_trajectory::parameter_labels()
    {
        std::vector<std::string> labels;
        for (int i = 0; i < _Npars; i++)  labels.push_back("par[" + std::to_string(i) + "]");
        return labels;
    };

    // -----------------------------------------------------------------------
    // Evaluation of trajectory

    complex raw_trajectory::evaluate(double s)
    {
        return _alphaSUB + DR_RHC(s) - DR_RHC(_sSUB);
    };

    // Output the imaginary part on the real line
    double raw_trajectory::imaginary_part(double s){  return std::imag(evaluate(s)); };
    double raw_trajectory::real_part(double s){ return std::real(evaluate(s)); };

    double raw_trajectory::width(double s)
    {
        if (s <= _sRHC) return 0.;
        auto f = [&](double x){ return real_part(x); };
        double alphaPrime = boost::math::differentiation::finite_difference_derivative(f, s);
        return imaginary_part(s) / alphaPrime / sqrt(s);
    };
    
    // -----------------------------------------------------------------------
    // Internal functions for evaluation

    // Evaluate dispersion relation on real axis with ieps perscriptions
    complex raw_trajectory::DR_RHC(double s)
    {
        // We split the integration in two parts
        // to properly handle the Principle Value and ieps perscription

        double RHCs = (s > _sRHC) ? RHC(s) : 0.;
        auto fdx = [this, s, RHCs](double x)
        {
            complex integrand;
            integrand  = RHC(x) - RHCs;
            integrand *= (s / x); // One subtraction
            integrand /= (x - s - IEPS); 
            return integrand;
        };

        // bounds of integration
        double low = _sRHC + EPS;
        double high = std::numeric_limits<double>::infinity();

        complex logarithm, integral, result;
        integral  = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(fdx, low, high, 20, 1.E-9, NULL);
        logarithm = are_equal(s, _sRHC, 1E-5) ? 0 : RHCs * log(1. - (s + IEPS) / low);
        result = (integral - logarithm) / M_PI;
    
        return result;
    };

};