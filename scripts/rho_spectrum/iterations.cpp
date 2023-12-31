// Test of iteration convergence from fit results. 
// Plots trajectory at each iterations up to n = 5
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

void iterations()
{   
    using namespace analyticRT;

    trajectory alpha = new_trajectory<iterative>(4.*M2_PION);
    alpha->set_option(1);
    
    std::vector<double> pars = {1.0, 0.532629947, 1.04260251, 43.4117231};

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
    
    plot replot = isovector_spins(plotter);
    replot.set_curve_points(200);

    auto plot_iter = [&](int n)
    {
        alpha->max_iterations(n);
        alpha->set_parameters(pars);
        replot.add_curve( {-2, 7.5}, realpha, "n = " + std::to_string(n));
        replot.add_dashed({-2, 7.5}, imalpha);
    };

    plot_iter(0);
    plot_iter(1);
    plot_iter(2);
    plot_iter(3);

    replot.set_legend(0.4, 0.7);
    replot.set_legend_spacing(0.02);
    replot.save("iterations.pdf");
};