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
#include "hyp_2F1.hpp"
#include "cgamma.hpp"


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

    std::array<double,3> lambdas, sths, cs;
    lambdas = {0.531836806,    2.33587715,     24.6842368};
    sths    = {0.0779191396,   2.12,           30};
    cs      = {0.143824451,    0.502572062,    37.1979956};

    auto ImFiore = [lambdas,sths,cs] (double s)
    {
        double result = 0;
        for (int i = 0; i < 3; i++)
        {
            if (s > sths[i]) result+= cs[i]*sqrt(s-sths[i])*pow(1-sths[i]/s, lambdas[i]);
        }
        return result;
    };

    auto ReFiore = [lambdas,sths,cs] (double s)
    {
        double result = 0.491;
        for (int i = 0; i < 3; i++)
        {
            if (s > sths[i]) result += std::real(2./sqrt(PI)*cs[i]*cgamma(lambdas[i] + 3./2)/cgamma(lambdas[i]+1.)*sqrt(sths[i])*hyp_2F1(-lambdas[i], 1., 3./2, sths[i]/s));
            else result += std::real(s/sqrt(PI)*cs[i]*cgamma(lambdas[i] + 3./2)/cgamma(lambdas[i]+2.)/sqrt(sths[i])*hyp_2F1(1., 1/2., lambdas[i]+2., s/sths[i]));
        }
        return result;
    };

    // Set up plotter
    plotter plotter;
    
    plot J_plot = plotter.new_plot();
    J_plot.set_curve_points(200);
    J_plot.add_curve( {-2, 7.5}, [alpha](double s){ return alpha->real_part(s); },      "Real");
    J_plot.add_dashed({-2, 7.5}, ReFiore);
    J_plot.add_curve( {-2, 7.5}, [alpha](double s){ return alpha->imaginary_part(s); }, "Imaginary");
    J_plot.add_dashed({-2, 7.5}, ImFiore);
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

    plot HE_plot = plotter.new_plot();
    HE_plot.set_curve_points(500);
    HE_plot.set_logscale(true, true);
    HE_plot.add_curve( {1,  1E5}, [alpha](double s){ return alpha->real_part(s); },      "Real");
    HE_plot.add_dashed({1,  1E5}, ReFiore);
    HE_plot.add_curve( {1,  1E5}, [alpha](double s){ return alpha->imaginary_part(s); }, "Imaginary");
    HE_plot.add_dashed({1,  1E5}, ImFiore);
    HE_plot.add_curve({1,   1E5}, [](double s){ return 0.5 + 0.9*s;}, "0.5 + 0.9 #it{s}" );
    HE_plot.set_labels("#it{s} [GeV^{2}]", "#alpha(#it{s})");
    HE_plot.save("heplot.pdf");
};