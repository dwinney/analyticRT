// Assemble a bunch of isobars together into a total amplitude with the correct
// crossing structure
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "amplitude.hpp"

namespace analyticRT
{
    complex raw_amplitude::direct_channel(unsigned int i, double s, double zs)
    {
        if (i > 2) return std::nan("");

        // Sum over direct channel pieces
        complex result = 0.;
        for (auto isobar : _isobars) if (isobar->isospin() == i) result += isobar->evaluate(s, zs);

        return result;
    };


    complex raw_amplitude::cross_channels(unsigned int i, double s, double zs)
    {
        if (i > 2) return std::nan("");
        if (_ignore_cross) return 0.;

        double t = t_man(s, zs), zt = z_t(s, zs);
        double u = u_man(s, zs), zu = z_u(s, zs);

        // Sum over cross channel pieces
        complex result = 0.;
        for (auto isobar : _isobars)
        {
            int ip = isobar->isospin();
            result += C[i][ip] * isobar->evaluate(t, zt);
            result += C[i][ip] * isobar->evaluate(u, zu) * pow(-1., i + ip);
        };

        return result;
    };

    complex raw_amplitude::direct_projection(unsigned int i, unsigned int j, double s)
    {
        if (i > 2) return std::nan("");

        // Sum over direct channel pieces
        complex result = 0.;
        for (auto isobar : _isobars) if (isobar->isospin() == i) result += isobar->direct_projection(j, s);

        return result;
    };

    complex raw_amplitude::cross_projection(unsigned int i, unsigned int j, double s)
    {
        if (i > 2) return std::nan("");
        if (_ignore_cross) return 0;

        // Sum over cross channel pieces
        complex result = 0.;
        for (auto isobar : _isobars) result += C[i][isobar->isospin()] * isobar->cross_projection(j, s);

        return result;
    };
};
