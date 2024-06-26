// Reproduction of the rho-a2 trajectory from [1] compared to fit results using sqrtlog
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:hep-ph/0011035
// ---------------------------------------------------------------------------

#include "spectrum_plots.hpp"
#include "trajectory.hpp"
#include "trajectories/fiore.hpp"
#include "trajectories/k_matrix.hpp"

void results()
{   
    using namespace analyticRT;

    // Set up Fiore-type trajectory with parameters from [1]
    trajectory alpha_fiore = new_trajectory<fiore>(4.*M2_PION, "Fiore");
    alpha_fiore->set_parameters({0.491, 0.140, 0.902, 28.031});

    trajectory alpha = new_trajectory<k_matrix>(4.*M2_PION, "Us");
    alpha->set_option(1);
    alpha->set_parameters({1.0, 0.532629947, 1.04260251, 43.4117231});

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