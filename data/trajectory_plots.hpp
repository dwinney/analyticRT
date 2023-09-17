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
        data_set isovectors = isovector_spectrum();

        plot re = plotter.new_plot();

        std::vector<double>  s = square_elementwise(isovectors._x);
        std::vector<double> ds = 2.*multiply_elementwise(isovectors._x, isovectors._dx);
        std::vector<double>  J = isovectors._y;
        std::vector<double> dJ = isovectors._dy;

        re.add_data({s, ds}, {J, dJ}, jpacColor::Orange);
        re.set_ranges({-2, 6.5}, {-1.5, 6.5});
        re.set_labels("#it{s} [GeV^{2}]", "Re #alpha(#it{s})");
        re.add_vertical(  0, {kBlack, kSolid});
        re.add_horizontal(0, {kBlack, kSolid});

        plot gam = plotter.new_plot();

        std::vector<double>  gamma = isovectors._z;
        std::vector<double> dgamma = isovectors._dz;

        gam.add_data({s, ds}, {gamma, dgamma}, jpacColor::Orange);
        gam.set_ranges({-2, 6.5}, {-0.05, 0.7});
        gam.set_labels("#it{s} [GeV^{2}]", "#Gamma(#it{s})");
        gam.add_vertical(  0, {kBlack, kSolid});
        gam.add_horizontal(0, {kBlack, kSolid});

        return {re, gam};
    };
};

#endif