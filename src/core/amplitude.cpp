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
    complex amplitude::evaluate(unsigned int i, double s, double zs)
    {
        if (i > 2) return std::nan("");

        // Sum over direct channel pieces
        complex result = 0.;
        for (auto isobar : _isobars)
        {
            if (isobar->isospin() == i) result += isobar->evaluate(s, zs);
        };

        if (_ignore_cross) return result;

        double t = t_man(s, zs), zt = z_t(s, zs);
        double u = u_man(s, zs), zu = z_u(s, zs);

        // Sum over cross channel pieces
        for (auto isobar : _isobars)
        {
            int ip = isobar->isospin();
            result += C[i][ip] * isobar->evaluate(t, zt);
            result += C[i][ip] * isobar->evaluate(u, zu) * pow(-1., i + ip);
        };

        return result;
    };

    complex amplitude::partial_wave(unsigned int i, unsigned int j, double s)
    {
        if (i > 2) return std::nan("");

        // Sum over direct channel pieces
        complex result = 0.;
        for (auto isobar : _isobars)
        {
            if (isobar->isospin() == i) result += isobar->direct_projection(j, s);
        };

        if (_ignore_cross) return result;


        // Sum over cross channel pieces
        for (auto isobar : _isobars)
        {
            int ip = isobar->isospin();
            result += C[i][ip] * isobar->cross_projection(j, s);
            result += C[i][ip] * isobar->cross_projection(j, s) * pow(-1., i + ip);
        };

        return result;
    };
};
