// Comparison of the resulting dispersive trajectory with a linear form
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "trajectory_plots.hpp"
#include "trajectories/iterative.hpp"

void compare()
{   
    using namespace analyticRT;

    trajectory alpha = new_trajectory<iterative>(4.*M2_PION);
    alpha->set_option(2);
    alpha->max_iterations(3);
    alpha->set_parameters({0.477, 1.06038145, 0.72811094, 0.700916037});

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

    auto reLine = [](double s){ return 0.5 + 0.9 * s; };
    
    // Set up plotter
    plotter plotter;
    
    plot low = isovector_spins(plotter);
    low.set_curve_points(200);
    low.add_curve( {-2, 7.5}, re(alpha), "Real");
    low.add_dashed({-2, 7.5}, [](double s){ return 0.5 + 0.9*s; });
    low.add_curve( {-2, 7.5}, im(alpha), "Imaginary");
    low.add_dashed({-2, 7.5}, [](double s){ return (s > 4.*M2_PION) ? 0.1*sqrt(s - 4.*M2_PION) : 0.; });
    low.set_legend(0.4, 0.7);
    low.set_legend_spacing(0.02);

    plot high = plotter.new_plot();
    high.set_logscale(true, true);
    high.set_curve_points(10000);
    high.set_labels("#it{s} [GeV^{2}]", "#alpha(#it{s})");
    high.add_curve( {0.1, 10000}, re(alpha), "Real");
    high.add_dashed({0.1, 10000}, [](double s){ return 0.5 + 0.9*s; });
    high.add_curve( {0.1, 10000}, im(alpha), "Imaginary");
    high.add_dashed({0.1, 10000}, [](double s){ return (s > 4.*M2_PION) ? 0.1*sqrt(s - 4.*M2_PION) : 0.; });

    plotter.combine({2,1}, {low, high}, "compare.pdf");
};