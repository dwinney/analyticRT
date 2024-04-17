// Comparison of the P-wave projection of hypergeometric isobar and pipi data
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobars/truncated.hpp"
#include "isobars/k_matrix.hpp"
#include "trajectories/unitary.hpp"

#include "kinematics.hpp"
#include "pipi.hpp"
#include "trajectory.hpp"
#include "isobar.hpp"
#include "plotter.hpp"
#include "fitter.hpp"

using namespace analyticRT;

void swave()
{
    using namespace analyticRT;
    using complex = std::complex<double>;

    // Global constants
    int iso = 0, J = 0;

    // --------------------------------------------------------------------------
    // We will compare to the inverse of the K-matrix amplitude

    isobar kmatrix = new_isobar<k_matrix>(iso, "K-matrix");

    // Fit parameters
    double a, b, c; 
    a = -0.33855517;
    b =  35.481475;
    c =  1.2369145;
    kmatrix->set_parameters({a, b, c});
    
    // --------------------------------------------------------------------------
    // Set up the unitary dispersive trajectory

    trajectory alpha = new_trajectory<unitary>(STH, iso, "Dispersive");
    iterable(alpha)->set_initial_pars({-8, 0.2}, 1);
    alpha->set_subtraction(STH, 0); // Move subtraction to threshold instead of 0

    // The trajectory defines an isobar
    isobar sigma = new_isobar<truncated>(iso, 4, alpha, "#it{I} = 0");

    double lam2 = 1.5, g_A, gamma, g_sig, c_sig, alpha0;

    g_sig  =  0.2;
    g_A    = -0.08;
    gamma  = 4.;
    c_sig  = 0.00001;
    alpha0 = -4.8*g_sig;

    sigma->set_parameters({g_sig, lam2});
    alpha->set_parameters({alpha0, lam2, g_A, g_sig, gamma, c_sig});

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    // P-wave real part vs data
    entry_style dashed;
    dashed._color = jpacColor::DarkGrey;
    dashed._style = kDashed;
    dashed._add_to_legend = true;
    dashed._label = "K-matrix";

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#alpha^{#sigma}_{#it{s}} / #it{g}_{#sigma}");
    p1.set_legend(0.70, 0.2);
    p1.set_ranges({0, 1}, {-5, 2});
    p1.add_curve(  {STH + EPS, 1},  [alpha, g_sig](double s){ return alpha->real_part(s) / g_sig;},      "Real");
    p1.add_curve(  {STH + EPS, 1},  [alpha, g_sig](double s){ return alpha->imaginary_part(s)/ g_sig;}, "Imaginary");
    p1.add_curve(  {STH + EPS, 1},  [kmatrix, J](double s){ return std::real(-1/kmatrix->direct_projection(J, s));}, dashed);
    dashed._add_to_legend = false;
    p1.add_curve(  {STH + EPS, 1},  [kmatrix, J](double s){ return std::imag(-1/kmatrix->direct_projection(J, s));}, dashed);
    p1.add_vertical(STH);

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#it{A}^{(0)}_{0}(s)");
    p2.set_legend(0.60, 0.2);
    p2.add_curve(  {STH, 0.8},  [sigma](double s){ return std::real(sigma->direct_projection(0, s)); }, "Real");
    p2.add_curve(  {STH, 0.8},  [sigma](double s){ return std::imag(sigma->direct_projection(0, s)); }, "Imaginary");
    dashed._add_to_legend = true;
    dashed._label = "GKPY";
    p2.add_curve( {STH + EPS, 0.8}, [](double s){ return std::real(pipi::partial_wave(0, 0, s));}, dashed);
    dashed._add_to_legend = false;
    p2.add_curve( {STH + EPS, 0.8}, [](double s){ return std::imag(pipi::partial_wave(0, 0, s));}, dashed);

    plotter.combine({2,1}, {p1, p2}, "a00.pdf");
};