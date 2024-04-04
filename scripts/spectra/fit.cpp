// conducts fit to isovector resonance masses and widths and outputs results
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
#include "trajectories/fiore.hpp"

using namespace analyticRT;

struct spectrum_fit
{
    static std::string data_type(int i){ return "Masses and Widths"; };

    // We add a small fictious error to the spins from the uncertainty in the masses
    // Therefore we minimize chi2 instead of difference of squares
    static double fcn(const std::vector<data_set> &data_sets, isobar iso, trajectory traj)
    { 
        double chi2 = 0; 
        for (auto &data : data_sets)
        {
            for (int i = 0; i < data._N; i++)
            {
                double s = pow(data._x[i], 2);

                double spin_th  = traj->real_part(s);
                double spin_ex  = data._y[i];
                double spin_err = data._dy[i];
                chi2 += pow((spin_th - spin_ex)/spin_err, 2);
                
                double width_th  = traj->width(s);
                double width_ex  = data._z[i];
                double width_err = data._dz[i];
                chi2 += pow((width_th - width_ex)/width_err, 2);
            };
        };
        return chi2; 
    };
};

void fit()
{   
    using namespace analyticRT;

    trajectory alpha = new_trajectory<k_matrix>(4.*M2_PION, "Iterated trajectory");
    alpha->set_option(2);

    fitter<spectrum_fit> fitter(nullptr, alpha);
    fitter.add_data( rho_spectrum() );
    fitter.add_data( a_spectrum() );
    
    fitter.fix_parameter("Lambda^2", 1.5);

    fitter.set_parameter_limits("alpha(0)", {0.,  1.});

    fitter.set_parameter_posdef("gamma");
    fitter.set_parameter_posdef("c[1]");
    fitter.set_parameter_posdef("c[2]");
    
    fitter.do_iterative_fit(5);

    // Set up plotter
    plotter plotter;
    
    plot J_plot = isovector_spins(plotter);
    J_plot.set_curve_points(200);
    J_plot.add_curve({-2, 7.5}, [alpha](double s){ return alpha->real_part(s); },      "Real");
    J_plot.add_curve({-2, 7.5}, [alpha](double s){ return alpha->imaginary_part(s); }, "Imaginary");
    J_plot.set_legend(0.4, 0.7);
    J_plot.set_legend_spacing(0.02);

    plot width_plot = isovector_widths(plotter);
    width_plot.color_offset(2);
    width_plot.set_curve_points(200);
    width_plot.add_curve({0., 7.5},  [alpha](double s){ return alpha->width(s); });

    // Output to file
    plotter.combine({2,1}, {J_plot, width_plot}, "results.pdf");
};