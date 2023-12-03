// Assemble a bunch of isobars together into a total amplitude with the correct
// crossing structure
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef AMPLITUDE_HPP
#define AMPLITUDE_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "isobar.hpp"

namespace analyticRT
{
    class amplitude
    {
        public: 
        
        amplitude(){};

        amplitude(std::vector<isobar> isobars)
        : _isobars(isobars)
        {};

        // Evaluate full isopin amplitude by iterating all the isobars
        complex evaluate(unsigned int i, double s, double zs);

        // Evaluate the full partial wave projection
        complex partial_wave(unsigned int i, unsigned int j, double s);

        // Ability to turn off the cross channel contributions
        inline void ignore_cross(bool x){ _ignore_cross = x; };

        inline void add_isobar(isobar x) { _isobars.push_back(x); };
        inline void clear_isobars(){ _isobars.clear(); };
        
        private:

        // Constituent isobars
        std::vector<isobar> _isobars;

        // Whether or not we ignore cross channe contributions
        bool _ignore_cross = false;
    };
};

#endif