// Comparison of the P-wave projection of hypergeometric isobar and pipi data
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobars/truncated.hpp"
#include "trajectories/unitary.hpp"
#include "spectrum_data.hpp"

#include "kinematics.hpp"
#include "pipi.hpp"
#include "trajectory.hpp"
#include "isobar.hpp"
#include "plotter.hpp"
#include "data_set.hpp"

using namespace analyticRT;

void replot()
{
    using namespace analyticRT;

    int iso = 1, J = 1;

    // --------------------------------------------------------------------------
    // For the I = 1 we can assume for now that there is a single trajectory

    auto guess = [](double s){ return (0.5 + 0.9*s)/sqrt(1. + s/20.); };

    // This defines the dispersive form
    trajectory alpha = new_trajectory<unitary>(J, guess, "#rho");
    alpha->set_integrator_depth(20);

    // The trajectory defines an isobar
    isobar rho = new_isobar<truncated>(iso, 5, alpha, "truncated, n = 5");

    // ---------------------------------------------------------------------------
    // Import parameters from file

    std::string dir = "scripts/rho/";
    std::string file_prefix = "rho_fit_2GeV_";

    auto iso_pars  = import_transposed<21>(dir + file_prefix + "iso_pars.txt");
    auto traj_pars = import_transposed<21>(dir + file_prefix + "traj_pars.txt");

    iterable(alpha)->iterate<21>(traj_pars, 21);
    rho->set_parameters(iso_pars.back());

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    data_set rhos = rho_spectrum();
    data_set as   = a_spectrum();

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#alpha_{#rho}(#it{s})");
    p1.set_ranges({-0.1, 1.5}, {0., 2.0});
    p1.set_legend(0.4, 0.7);
    p1.add_curve(  {-0.1, 1.5},      [alpha](double s){ return alpha->real_part(s);} ,      "Real");
    p1.add_curve(  {-0.1, 1.5},      [alpha](double s){ return alpha->imaginary_part(s); }, "Imaginary");
    p1.add_curve(  {-0.1, 1.5},      [alpha](double s){ return 0.5+0.9*s; }, dashed(jpacColor::DarkGrey, "0.5 + 0.9 #it{s}"));
    p1.add_data({square_elementwise(rhos._x), {}}, {rhos._y, {}}, jpacColor::DarkGrey);
    p1.add_vertical(  0, {kBlack, kSolid});
    p1.shade_region({STH,1});
    p1.save("timelike.pdf");

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#it{f}_{#kern[-0.7]{1}}^{#kern[-0.7]{1}}(#it{s})");
    p2.set_ranges({0, 1}, {-0.7, 1.3});
    p2.set_legend(0.25, 0.7);
    p2.color_offset(2);
    p2.add_curve(  {0,1},      [rho](double s){ return std::real(rho->direct_projection(1, s));} , "Real");
    p2.add_curve(  {0,1},      [rho](double s){ return std::imag(rho->direct_projection(1, s)); }, "Imaginary");
    p2.add_curve( {0, 1},    [rho](double s){ return (s > STH) ? sqrt(1.- STH/s) * std::norm(rho->direct_projection(1,s)) : 0; }, dashed(jpacColor::Orange, "Exact Unitarity"));
    p2.add_curve( {STH + EPS, 1}, [](double s){ return std::real(pipi::partial_wave(1, 1, s));}, dashed(jpacColor::DarkGrey, "GKPY"));
    p2.add_curve( {STH + EPS, 1}, [](double s){ return std::imag(pipi::partial_wave(1, 1, s));}, dashed(jpacColor::DarkGrey));
    p2.save("pw_rho.pdf");

    auto dat = import_data<4>("data/charge_exchange.dat");
    plot p3 = plotter.new_plot();
    p3.set_labels("#it{s}  [GeV^{2}]", "Re #alpha_{#rho}(#it{s})");
    p3.set_ranges({-1.5, 0}, {-0.7, 0.6});
    p3.add_curve(  {-1.5, 0},      [alpha](double s){ return alpha->real_part(s);});
    p3.add_curve(  {-1.5, 0},      [alpha](double s){ return 0.5+0.9*s; }, dashed(jpacColor::DarkGrey));
    p3.add_data({-dat[0], dat[1]/2}, {dat[2], dat[3]});
    p3.save("spacelike.pdf");
};