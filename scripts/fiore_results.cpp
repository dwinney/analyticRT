// Reproduction of the rho-a2 trajectory from [1]
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

void fiore_results()
{   
    using namespace analyticRT;

    // Set up Fiore-type trajectory with parameters from [1]
    trajectory alpha = new_trajectory<fiore>(4.*M2_PION);
    alpha->set_parameters({0.491, 0.140, 0.902, 28.031});

    // Curves to be plotted below
    auto realpha = [alpha](double s)
    {
        return alpha->real_part(s);
    };
    auto imalpha = [alpha](double s)
    {
        return alpha->imaginary_part(s);
    };  
    auto width = [alpha](double s)
    {
        return alpha->width(s);
    };  

    // Set up plotter
    plotter plotter;
    auto ps = isovector_plots(plotter);
    
    plot J_plot = ps[0];
    J_plot.set_curve_points(200);
    J_plot.add_curve({-2, 7.5}, realpha, "Real");
    J_plot.add_curve({-2, 7.5}, imalpha, "Imaginary");
    J_plot.set_legend(0.4, 0.7);
    J_plot.set_legend_spacing(0.02);

    plot width_plot = ps[1];
    width_plot.color_offset(2);
    width_plot.set_curve_points(200);
    width_plot.add_curve({0., 7.5}, width);

    // Output to file
    plotter.combine({2,1}, {J_plot, width_plot}, "fiore.pdf");
};