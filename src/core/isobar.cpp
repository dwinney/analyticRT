// Define an isobar, these work like the full amplitude but specify
// only a single channel and therefore fixed isospin
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobar.hpp"

namespace analyticRT
{
    // Make sure the parameters have the correct size and then feed it to allocate_parameters
    void raw_isobar::set_parameters(std::vector<double> pars)
    {
        if (pars.size() != N_free()) 
        {
            warning(id(), "Incorrect number of parameters passed! Expected " + std::to_string(N_free()) + 
                          " but received " + std::to_string(pars.size()) + "!");
            return;
        };

        allocate_parameters(pars);
    };

    // By default we calculate the partial wave projection numerically
    complex raw_isobar::direct_projection(unsigned int j, double s)
    {
        auto dF = [j, s, this](double z)
        {
            return legendre_P(j, z) * this->evaluate(s, z) / 2.;
        };

        return boost::math::quadrature::gauss_kronrod<double, 15>::integrate(dF, -1, 1, 5, 1.E-9, NULL);
    };

    complex raw_isobar::cross_projection(unsigned int j, double s)
    {
        // No factor of 2 because t and u channels are summed
        auto dF = [j, s, this](double z)
        {
            return legendre_P(j, z) * this->evaluate(t(s, z), z_t(s, z));
        };

        return boost::math::quadrature::gauss_kronrod<double, 15>::integrate(dF, -1, 1, 5, 1.E-9, NULL);
    };
};