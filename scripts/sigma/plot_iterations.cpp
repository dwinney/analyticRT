// Comparison of the S-wave projection of hypergeometric isobar and pipi data
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobars/truncated.hpp"
#include "trajectories/unitary.hpp"

#include <sstream>
#include "kinematics.hpp"
#include "pipi.hpp"
#include "trajectory.hpp"
#include "isobar.hpp"
#include "plotter.hpp"
#include "fitter.hpp"
#include "data_set.hpp"

using namespace analyticRT;

void plot_iterations()
{
    using namespace analyticRT;
    using complex = std::complex<double>;

    // Global constants
    int iso = 0, J = 0;

    // -------------------------------- ------------------------------------------
    // Set up the unitary dispersive trajectory
    auto guess = [](double s){ return (-0.2 + 0.1*s)/sqrt(1.+s/20.); };

    trajectory alpha = new_trajectory<unitary>(iso, guess, "sigma");
    alpha->set_option(unitary::kAddConstant);
    // alpha->set_integrator_depth(10);
    iterable(alpha)->set_interp_pars(600, {50, 1500});

    // The trajectory defines an isobar
    isobar f0 = new_isobar<truncated>(iso, 0, alpha, "I = 0 only");
    f0->set_option(truncated::kAddConstant);

    // -------------------------------- ------------------------------------------
    // Instead of fitting to the partial waves we fit to the trajectory

    // Iterate once to begin fitting
    std::vector<std::vector<double>> pars; 

    pars.push_back({2, -0.181640798, 0.227942035, 1.16959135, 0.522641701, 1.05323747});
    pars.push_back({2, -0.15564638,  0.361631951, 2.24314035, 0.111740855, 0.00377618901});
    pars.push_back({2, -0.145902785, 0.532873878, 3.58606037, 0.082107926, 0.0000211061485});

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#it{f}^{0}_{0}(#it{s})");
    p1.set_legend(0.25, 0.7);
    p1.set_ranges({0, 0.9}, {-0.3, 1.6});

    double smax = 1.5;
    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#alpha_{#sigma}(#it{s})");
    p2.set_ranges({-0.2, smax}, {-0.25, 0.23});
    p2.set_curve_points(200);
    p2.set_legend(0.3, 0.7);
    p2.shade_region({STH, 0.6});
    p2.add_vertical(  0, {kBlack, kSolid});
    p2.add_horizontal(0, {kBlack, kSolid}); 

    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::real(pipi::partial_wave(0, 0, s));}, solid(jpacColor::DarkGrey, "GKPY"));
    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::imag(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey        ));

    p2.add_curve( {-0.2, smax}, guess, solid(jpacColor::DarkGrey, "Initial guess"));
    p2.add_curve( {-0.2, smax}, [&](double s){return 0.;}, dashed(jpacColor::DarkGrey));

    std::vector<std::string> labels = {"0^{th} iter.", "1^{st} iter.", "2^{nd} iter."};
    for (int i = 0; i < pars.size(); i++) 
    { 
        alpha->set_parameters(pars[i]); 
        f0->set_parameters({pars[i][0], pars[i][2], pars[i][3]});

        p1.add_curve(  {0, 0.9},   [&](double s){ return std::real(f0->direct_projection(0, s)); }, labels[i]);
        p1.add_dashed( {0, 0.9},   [&](double s){ return std::imag(f0->direct_projection(0, s)); });

        p2.add_curve(  {-0.2,  smax}, [&](double s){ return alpha->real_part(s); },                labels[i]);
        p2.add_dashed( {-0.2,  smax}, [&](double s){ return alpha->imaginary_part(s); });
        alpha->iterate(); 
    };

    p1.save("f00_iters.pdf");
    p2.save("alpha_sigma_iters.pdf");
};