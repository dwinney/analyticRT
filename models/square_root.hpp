
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
            set_Npars(3);
        };

        protected:

        // RHC given by a simple square root 
        inline double RHC(double s)
        {
            return sqrt(s - _sRHC);
        };
    };
};