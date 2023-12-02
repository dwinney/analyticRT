// Dispersive isobar given by omnes function and pi pi scattering phase shift
// 
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobar.hpp"

namespace analyticRT
{
    class omnes : public raw_isobar
    {
        public: 

        // Explicitly only allow a RHC
        omnes(key x, unsigned int isospin, int j, std::string id)
        : raw_isobar(x, isospin, id), _j(j)
        {
            // NO free parameters
            set_Nfree(0);
        };

        // Evaluate the full term with angular dependence
        complex evaluate(double s, double zs)
        {
            ////////////////////////////////
        };

        protected:

        int _j = -1; 
    };
};