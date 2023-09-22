#include "trajectory_data.hpp"
#include "trajectory_plots.hpp"
#include "fitter.hpp"
#include "trajectory.hpp"
#include "sqrtlog.hpp"

void iterations()
{   
    using namespace analyticRT;

    trajectory alpha = new_trajectory<sqrtlog>(4.*M2_PION);
    alpha->set_subtraction(0., 0.477);

    std::vector<double> pars = {1.06044486, 2.15706648,  2.07668753};

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