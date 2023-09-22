#include "trajectory_data.hpp"
#include "trajectory_plots.hpp"
#include "fitter.hpp"
#include "trajectory.hpp"
#include "sqrtlog.hpp"
#include "fiore.hpp"

void fit()
{   
    using namespace analyticRT;

    trajectory alpha = new_trajectory<sqrtlog>(4.*M2_PION, "Sqrt-Log trajectory");
    alpha->set_subtraction(0., 0.477);

    fitter fitter(alpha);
    fitter.add_data( isovector_spectrum() );
    
    fitter.set_parameter_limits("gamma", {0., 5.});
    fitter.set_parameter_limits("c[1]",  {0., 10.});
    fitter.set_parameter_limits("c[2]",  {0., 10.});

    fitter.iterative_fit(10);

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