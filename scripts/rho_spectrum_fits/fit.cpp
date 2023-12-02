// conducts fit to isovector resonance masses and widths and outputs results
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:hep-ph/0011035
// ---------------------------------------------------------------------------

#include "trajectory_data.hpp"
#include "trajectory_plots.hpp"
#include "fitter.hpp"
#include "trajectories/iterative.hpp"
#include "trajectories/fiore.hpp"

void fit()
{   
    using namespace analyticRT;

    trajectory alpha = new_trajectory<iterative>(4.*M2_PION, "Iterated trajectory");
    alpha->set_option(1);

    fitter fitter(alpha);
    fitter.add_data( isovector_spectrum() );
    
    fitter.set_parameter_limits("alpha(0)", {0.,1});
    fitter.set_parameter_limits("gamma", {0.1, 2.});
    fitter.set_parameter_limits("c[1]",  {0., 10.});
    // fitter.set_parameter_limits("c[2]",  {0., 10.});
    
    fitter.iterative_fit(3);

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
    
    plot J_plot = isovector_spins(plotter);
    J_plot.set_curve_points(200);
    J_plot.add_curve({-2, 7.5}, realpha, "Real");
    J_plot.add_curve({-2, 7.5}, imalpha, "Imaginary");
    J_plot.set_legend(0.4, 0.7);
    J_plot.set_legend_spacing(0.02);

    plot width_plot = isovector_widths(plotter);
    width_plot.color_offset(2);
    width_plot.set_curve_points(200);
    width_plot.add_curve({0., 7.5}, width);

    // Output to file
    plotter.combine({2,1}, {J_plot, width_plot}, "results.pdf");
};