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

#include "spectrum_data.hpp"
#include "spectrum_plots.hpp"
#include "fitter.hpp"
#include "trajectories/k_matrix.hpp"

void iterations()
{   
    using namespace analyticRT;

    trajectory alpha = new_trajectory<k_matrix>(4.*M2_PION, "Iterative");
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

    int n = 0;
    auto plot_iter = [&]()
    {
        alpha->set_parameters(pars);
        replot.add_curve( {-2, 7.5}, realpha, "n = " + std::to_string(n));
        replot.add_dashed({-2, 7.5}, imalpha);
        alpha->iterate(); n++;
    };

    for (int i = 0; i <= 5; i++) plot_iter();

    replot.set_legend(0.4, 0.7);
    replot.set_legend_spacing(0.02);
    replot.save("iterations.pdf");
};