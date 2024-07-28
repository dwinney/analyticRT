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
#include "trajectories/polynomial.hpp"
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
                chi2 += pow(spin_th - spin_ex, 2);
                
                double width_th  = width(traj, s);
                double width_ex  = data._z[i];
                chi2 += pow(width_th - width_ex, 2);
            };
        };
        return chi2; 
    };
    static double width(trajectory alpha, double s)
    {
        return alpha->imaginary_part(s) / 0.9 / sqrt(s);
    };  
};

void fit()
{   
    using namespace analyticRT;

    auto guess = [](double s){ return (0.5 + s)/sqrt(1. + s/10.); };

    trajectory alpha = new_trajectory<polynomial>(4*M2_PION, 1, guess, "Iterated trajectory");
    alpha->set_option(0);

    fitter<spectrum_fit> fitter(nullptr, alpha);

    data_set rhos = rho_spectrum();
    data_set as   = a_spectrum();
    fitter.add_data( rhos );
    fitter.add_data( as );
    
    fitter.fix_parameter("Lambda^2", 2.0);
    fitter.fix_parameter("alpha(0)", 0.5);

    fitter.set_parameter_posdef("gamma");
    fitter.set_parameter_posdef("c[alpha]");
    
    fitter.do_iterative_fit({1.3, 1.4}, 20);


    trajectory alpha_fiore = new_trajectory<fiore>(4.*M2_PION, "Fiore");
    alpha_fiore->set_parameters({0.491, 0.140, 0.902, 28.031});
    alpha_fiore->set_integrator_depth(25);

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
    
    plot J_plot = plotter.new_plot();
    J_plot.set_curve_points(200);
    J_plot.add_curve( {-2, 7.5}, re(alpha), "Real");
    J_plot.add_dashed({-2, 7.5}, re(alpha_fiore));
    J_plot.add_curve( {-2, 7.5}, im(alpha), "Imaginary");
    J_plot.add_dashed({-2, 7.5}, im(alpha_fiore));
    J_plot.add_curve({-2, 7.5}, [](double s){ return 0.5 + 0.9*s;}, "0.5 + 0.9 #it{s}" );
    J_plot.set_legend(0.4, .7);
    J_plot.set_legend_spacing(0.02);

    J_plot.add_data({square_elementwise(rhos._x), {}}, {rhos._y, {}}, jpacColor::DarkGrey);
    J_plot.add_data({square_elementwise(as._x), {}}, {as._y, {}}, jpacColor::DarkGrey);

    J_plot.set_ranges({-2, 7.5}, {-1.5, 8.});
    J_plot.set_labels("#it{s} [GeV^{2}]", "#alpha(#it{s})");
    J_plot.add_vertical(  0, {kBlack, kSolid});
    J_plot.add_horizontal(0, {kBlack, kSolid}); 
    J_plot.save("jplot.pdf");

    // plot width_plot = isovector_widths(plotter);
    // width_plot.color_offset(2);
    // width_plot.set_curve_points(200);
    // width_plot.add_curve( {0., 7.5}, gamma(alpha));
    // width_plot.add_dashed({0., 7.5}, gamma(alpha_fiore));

    plot HE_plot = plotter.new_plot();
    HE_plot.set_curve_points(300);
    HE_plot.set_logscale(true, true);
    HE_plot.add_curve( {1,  1E5},  re(alpha), "Real");
    HE_plot.add_dashed( {1, 1E5}, re(alpha_fiore));
    HE_plot.add_curve( {1,  1E5}, im(alpha), "Imaginary");
    HE_plot.add_dashed( {1, 1E5}, im(alpha_fiore));
    HE_plot.add_curve({1,   1E5}, [](double s){ return 0.5 + 0.9*s;}, "0.5 + 0.9 #it{s}" );
    HE_plot.save("heplot.pdf");

    // Output to file
    plotter.combine({2,1}, {J_plot, HE_plot}, "results.pdf");
};