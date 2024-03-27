// Comparison of the P-wave projection of hypergeometric isobar and pipi data
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobars/hypergeometric.hpp"
#include "isobars/truncated.hpp"
#include "trajectories/unitary.hpp"
#include "trajectories/linear.hpp"

#include "trajectory_plots.hpp"
#include "kinematics.hpp"
#include "trajectory.hpp"
#include "amplitude.hpp"
#include "pipi.hpp"
#include "print.hpp"
#include "plotter.hpp"

void pwave()
{
    using namespace analyticRT;

    // ---------------------------------------------------------------------------
    // The overall pi pi amplitude is accessed through amplitude class

    amplitude amp = new_amplitude("Full");

    // ---------------------------------------------------------------------------
    // For the I = 1 we can assume for now that there is a single trajectory

    // Given by the dispersive form
    trajectory alpha = new_trajectory<unitary>(4.*M2_PION, 1, "#rho");

    // The trajectory defines an isobar
    isobar rho = new_isobar<truncated>(1, 5, alpha, "I = 1");

    // Register our isobar into the full amplitude
    amp->add_isobar(rho);
    amp->ignore_cross(true); // We'll ignore the cross-channel projections in the partial wave

    // Unitarity forces the parameters of the trajectory and isobar
    // To be related
    double lambda2 = 1.5;   // Regge scale

    // Couplings if amp ignores cross channels
    double alpha0  = 0.557; // Trajectory intercept
    double g       = 145;   // Coupling constant
    double h       = 1.05;  // Asymptotic constant

    rho->set_parameters(  {g, lambda2});
    alpha->max_iterations(0);
    alpha->set_initial_guess({0.557, 0.9}, 20);
    alpha->set_parameters({alpha0, lambda2, g, h});

    print(alpha->evaluate(1)); exit(1);

    // ---------------------------------------------------------------------------
    // Make plot

    std::array<double,2> plot_bounds = {EPS, 2.0};
    plotter plotter;

    // Trajectory vs rho masses
    plot p1 = plotter.new_plot();
    p1.set_legend_spacing(0.03);
    p1.set_legend(0.25, 0.7);
    p1.set_curve_points(200);
    p1.add_curve({0,  1}, [alpha](double s){ return alpha->real_part(s); },   "Real");
    p1.add_curve({0,  1}, [alpha](double s){ return alpha->imaginary_part(s); },             "Imaginary");
    p1.save("traj.pdf");

    // // P-wave real part vs data
    // entry_style data_style;
    // data_style._color = jpacColor::DarkGrey;
    // data_style._style = kDashed;
    // data_style._label = "Data";

    // plot p3 = plotter.new_plot();
    // p3.set_labels("#it{s}  [GeV^{2}]", "Re _{1}(s)");
    // p3.set_legend(0.64, 0.20);
    // p3.add_curve(  plot_bounds, [amp](double s){ return std::real(amp->partial_wave     (1, 1, s)); }, amp->id());
    // p3.add_curve(  plot_bounds, [amp](double s){ return std::real(amp->direct_projection(1, 1, s)); }, "Direct Channel");
    // p3.add_curve(  plot_bounds, [amp](double s){ return std::real(amp->cross_projection (1, 1, s)); }, "Cross Channel");
    // p3.add_curve( {STH + EPS, plot_bounds[1]}, [](double s){ return 16*PI*std::real(pipi::partial_wave(1, 1, s)); }, data_style);

    // plot p4 = plotter.new_plot();
    // p4.color_offset(1);
    // p4.set_labels("#it{s}  [GeV^{2}]", "Im f_{1}(s)");
    // p4.set_legend(0.64, 0.65);
    // p4.add_curve(  plot_bounds, [amp](double s){ return std::imag(amp->partial_wave(1, 1, s)); }, amp->id());
    // p4.add_dashed({STH + EPS, plot_bounds[1]}, [amp](double s){ return phase_space(s)*std::norm(amp->partial_wave(1, 1, s)); } );
    // p4.add_curve( {STH + EPS, plot_bounds[1]}, [](double s){ return 16*PI*std::imag(pipi::partial_wave(1, 1, s)); }, data_style);
    // plotter.combine({3,1}, {p1, /* p2, */ p3, p4}, "pwave.pdf");
};