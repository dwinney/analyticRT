
// Trajectory resulting from a square root RHC
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "trajectory.hpp"

namespace analyticRT
{
    class square_root : public raw_trajectory
    {
        public: 

        // Explicitly only allow a RHC
        square_root(double R, std::string id)
        : raw_trajectory(R, id)
        {
            set_Npars(2);
        };

        // Parameters are the intercept and coupling
        inline void allocate_parameters(std::vector<double> pars)
        {
            set_subtraction(0., pars[0]);
            _coupling = pars[1];
        };

        inline std::vector<std::string> parameter_labels()
        {
            return {"alpha(0)", "coupling"};
        };

        protected:

        // RHC given by a simple square root 
        inline double RHC(double s)
        {
            return _coupling * sqrt(s - _sRHC);
        };

        // Coupling constant to the threshold
        double _coupling = 0;
    };
};