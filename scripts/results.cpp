// Reproduction of the rho-a2 trajectory from [1] compared to fit results using sqrtlog
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:hep-ph/0011035
// ---------------------------------------------------------------------------

#include "trajectory_plots.hpp"
#include "trajectory.hpp"
#include "fiore.hpp"
#include "sqrtlog.hpp"

void results()
{   
    using namespace analyticRT;

    // Set up Fiore-type trajectory with parameters from [1]
    trajectory alpha_fiore = new_trajectory<fiore>(4.*M2_PION);
    alpha_fiore->set_parameters({0.491, 0.140, 0.902, 28.031});

    trajectory alpha = new_trajectory<sqrtlog>(4.*M2_PION);
    alpha->set_subtraction(0., 0.477);
    alpha->set_option(2);
    alpha->max_iterations(3);
    alpha->set_parameters({1.06038145, 0.72811094, 0.700916037});

    auto re = [](trajectory alpha)
    {
        return [alpha](double s){return alpha->real_part(s);};
    };
    auto im = [](trajectory alpha)
    {
        return [alpha](double s){return alpha->imaginary_part(s);};
    };
    auto gamma = [](trajectory alpha)
    {
        return [alpha](double s){return alpha->width(s);};
    }; 

    // Set up plotter
    plotter plotter;
    
    plot J_plot = isovector_spins(plotter);
    J_plot.set_curve_points(200);
    J_plot.add_curve( {-2, 7.5}, re(alpha), "Real");
    J_plot.add_dashed({-2, 7.5}, re(alpha_fiore));
    J_plot.add_curve( {-2, 7.5}, im(alpha), "Imaginary");
    J_plot.add_dashed({-2, 7.5}, im(alpha_fiore));
    J_plot.set_legend(0.4, 0.7);
    J_plot.set_legend_spacing(0.02);

    plot width_plot = isovector_widths(plotter);
    width_plot.color_offset(2);
    width_plot.set_curve_points(200);
    width_plot.add_curve( {0., 7.5}, gamma(alpha));
    width_plot.add_dashed({0., 7.5}, gamma(alpha_fiore));

    // Output to file
    plotter.combine({2,1}, {J_plot, width_plot}, "results.pdf");
};