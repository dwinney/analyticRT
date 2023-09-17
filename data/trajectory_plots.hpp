// Methods to shortcut plotting curves on top of data
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef TRAJECTORY_PLOTS_HPP
#define TRAJECTORY_PLOTS_HPP

#include "trajectory_data.hpp"
#include "elementwise.hpp"
#include "plotter.hpp"
#include "plot.hpp"

namespace analyticRT
{

    inline std::vector<plot> isovector_plots(plotter & plotter)
    {
        data_set rhos = rho_spectrum();
        data_set as   = a_spectrum();

        plot re = plotter.new_plot();

        std::vector<double> r_s  = square_elementwise(rhos._x);
        std::vector<double> r_ds = 2.*multiply_elementwise(rhos._x, rhos._dx);
        std::vector<double> r_J  = rhos._y;
        std::vector<double> r_dJ = rhos._dy;
        re.add_data({r_s, r_ds}, {r_J, r_dJ}, jpacColor::Orange);

        std::vector<double> a_s  = square_elementwise(as._x);
        std::vector<double> a_ds = 2.*multiply_elementwise(as._x, as._dx);
        std::vector<double> a_J  = as._y;
        std::vector<double> a_dJ = as._dy;
        re.add_data({a_s, a_ds}, {a_J, a_dJ}, jpacColor::Orange);

        re.set_ranges({-2, 7.5}, {-1.5, 8.});
        re.set_labels("#it{s} [GeV^{2}]", "#alpha(#it{s})");
        re.add_vertical(  0, {kBlack, kSolid});
        re.add_horizontal(0, {kBlack, kSolid});

        plot gam = plotter.new_plot();

        std::vector<double> r_gamma  = rhos._z;
        std::vector<double> r_dgamma = rhos._dz;
        gam.add_data({r_s, r_ds}, {r_gamma, r_dgamma}, jpacColor::Orange);

        std::vector<double> a_gamma  = as._z;
        std::vector<double> a_dgamma = as._dz;
        gam.add_data({a_s, a_ds}, {a_gamma, a_dgamma}, jpacColor::Orange);

        gam.set_ranges({0., 7.5}, {0., 0.8});
        gam.set_labels("#it{s} [GeV^{2}]", "#Gamma(#it{s})  [GeV]");

        return {re, gam};
    };
};

#endif