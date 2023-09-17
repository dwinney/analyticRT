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
        allocate_parameters(pars);
    };

    // Default parameter labels are just par[0], par[1], ...
    std::vector<std::string> raw_trajectory::parameter_labels()
    {
        std::vector<std::string> labels;
        for (int i = 0; i < _Npars; i++)
        {
            labels.push_back("par[" + std::to_string(i) + "]");
        }
        return labels;
    };

    // -----------------------------------------------------------------------
    // Evaluation of trajectory

    // Output the imaginary part on the real line
    double raw_trajectory::imaginary_part(double s)
    {
        if (s >= _sRHC) return RHC(s);
        if (s <= _sLHC) return LHC(s);
        else return 0.;
    };
};