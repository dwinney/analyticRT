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
    // ---------------------------------------------------------------------------
    // Parameter handling 

    // Make sure the parameters have the correct size and then feed it to allocate_parameters
    void raw_isobar::set_parameters(std::vector<double> pars)
    {
        if (pars.size() != N_pars()) 
        {
            warning(id(), "Incorrect number of parameters passed! Expected " + std::to_string(N_pars()) + 
                          " but received " + std::to_string(pars.size()) + "!");
            return;
        };

        allocate_parameters(pars);
    };

    void raw_isobar::allocate_parameters(std::vector<double> x)
    {
        // For each subamplitude grab a subvector of its parameters
        auto running_iter = x.begin();
        for (isobar amp : _isobars)
        {
            auto sub_pars = std::vector<double>(running_iter, running_iter + amp->N_pars());
            amp->set_parameters(sub_pars);

            running_iter += amp->N_pars();
        };
        return;
    };

    // ---------------------------------------------------------------------------
    // Default virtual methods 

    // The default evaluation is to iterate over sub-isobars
    complex raw_isobar::evaluate(double s, double zs)
    {
        complex result = 0.;
        for (auto x : _isobars)
        {
            result += x->evaluate(s, zs);
        };
        return result;
    };

    // By default we calculate the partial wave projection numerically
    complex raw_isobar::direct_projection(unsigned int j, double s)
    {
        auto dF = [j, s, this](double z)
        {
            return legendre_P(j, z) * this->evaluate(s, z) / 2.;
        };

        if (_adaptive) return boost::math::quadrature::gauss_kronrod<double, 15>::integrate(dF, -1, 1, 5, 1.E-9, NULL);
        return boost::math::quadrature::gauss<double, 15>::integrate(dF, -1, 1);
    };

    complex raw_isobar::cross_projection(unsigned int j, double s)
    {
        // No factor of 2 because t and u channels are summed
        auto dF = [j, s, this](double z)
        {
            return legendre_P(j, z) * this->evaluate(t_man(s, z), z_t(s, z));
        };

        if (_adaptive) return boost::math::quadrature::gauss_kronrod<double, 15>::integrate(dF, -1, 1, 5, 1.E-9, NULL);
        return boost::math::quadrature::gauss<double, 15>::integrate(dF, -1, 1);
    };

    // ------------------------------------------------------------------------------
    // External methods for sums of isobars
    isobar operator+(isobar a, isobar b)
    {
        bool are_compatible = (a->isospin() == b->isospin());
        if (!are_compatible)
        {
            warning("Tried adding two non-compatible isobars (" + a->id() + " and " + b->id() + ")!");
            return nullptr;
        }

        std::string id = a->id() + " + " + b->id();

        // This allows the saved pointers to all ways be "base" level isobars in stead of other sums
        std::vector<isobar> from_a = get_isobars(a);
        std::vector<isobar> from_b = get_isobars(b);

        from_a.insert(from_a.end(), from_b.begin(), from_b.end());

        // Initialize and return a new isobar
        return new_isobar<raw_isobar>(a->isospin(), from_a, id);
    };
};