// Define an amplitude these are meant to describe the full crossing symmetric
// pi pi scattering. 
// As such it can be a single isobar or sum of isobars in different channels
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef AMPLITUDE_HPP
#define AMPLITUDE_HPP

#include "isobar.hpp"
#include "constants.hpp"

namespace analyticRT
{
    class raw_amplitude;
    using amplitude = std::shared_ptr<raw_amplitude>;

    amplitude operator=(isobar x)
    {
        amplitude a;
        a->_isosbars.push_back(x);
        return a;
    };

    class raw_amplitude : public raw_isobar
    {
    };
};

#endif