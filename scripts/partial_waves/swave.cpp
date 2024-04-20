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

struct pipi_fit
{
    static std::string data_type(int i){ return "GKPY PW"; };

    // Sum difference of squares for both real and imaginary part
    static double fcn(const std::vector<data_set> &data_sets, isobar iso, trajectory traj)
    { 
        double dos = 0; // difference of squares 
        for (auto &data : data_sets)
        {
            for (int i = 0; i < data._N; i++)
            {
                double s     = data._x[i];
                std::complex<double> ex(data._y[i], data._z[i]); 
                std::complex<double> th = iso->direct_projection(0, s);

                dos += std::norm( std::imag(th) - std::imag(ex) );
                dos += std::norm( std::real(th) - std::real(ex) );
            };
        };
        return dos; 
    };
};


void swave()
{
    using namespace analyticRT;
    using complex = std::complex<double>;

    // Global constants
    int iso = 0, J = 0;
    
    // --------------------------------------------------------------------------
    // Compare with the "trajectory" from the K-matrix

    isobar kmatrix = new_isobar<k_matrix>(0, "K-matrix");
    kmatrix->set_parameters({0.33855466, -4.2438229, 1.2369133});

    // -------------------------------- ------------------------------------------
    // Set up the unitary dispersive trajectory

    auto initial_sigma = [](double s){ return (-0.7 + 1.0*s) / sqrt(1 + s/10); };

    trajectory alpha = new_trajectory<unitary>(iso, initial_sigma, "Dispersive");
    alpha->set_subtraction(STH, 0); // Move subtraction to threshold instead of 0

    // The trajectory defines an isobar
    isobar sigma = new_isobar<truncated>(iso, 4, alpha, "#it{I} = 0");
    sigma->set_option(truncated::kAddAdlerZero);

    data_set pipi_pwave = pipi::partial_wave(iso, J, 10, {0.1, 0.7});

    fitter<pipi_fit> fitter(sigma, alpha);
    fitter.set_parameter_labels({"lam2 (iso)", "g (iso)", "g_A", "lam2", "alpha(sth)", "g", "gamma", "c"});
    fitter.add_data( pipi_pwave );

    // Sync isobar's parameters to the trajectory as required by unitarity
    fitter.sync_parameter("g (iso)",    "g");
    fitter.sync_parameter("lam2 (iso)", "lam2");
    
    fitter.fix_parameter("lam2", 4.);
    fitter.fix_parameter("g",    1.0);
    fitter.set_parameter_posdef("gamma");
    
    fitter.do_fit({0.36526, -4.4484817, 3.3043884, 109.54624});

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#tilde{#alpha}^{#sigma}_{#it{s}} / #it{g}_{#sigma}");
    p1.set_legend(0.70, 0.2);
    p1.set_ranges({0, 3}, {-5, 4});
    p1.add_curve(  {STH + EPS, 3},  [sigma](double s){   return std::real(-1/sigma->direct_projection(0, s));},      "Real");
    p1.add_curve(  {STH + EPS, 3},  [sigma](double s){   return std::imag(-1/sigma->direct_projection(0, s));}, "Imaginary");
    p1.add_curve(  {STH + EPS, 3},  [kmatrix](double s){ return std::real(-1/kmatrix->direct_projection(0, s));}, dashed(jpacColor::DarkGrey, "K-matrix"));
    p1.add_curve(  {STH + EPS, 3},  [kmatrix](double s){ return std::imag(-1/kmatrix->direct_projection(0, s));}, dashed(jpacColor::DarkGrey));

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#it{A}^{(0)}_{0}(s)");
    p2.set_legend(0.450, 0.2);
    p2.set_ranges({STH, 0.9}, {-0.5, 1.5});
    p2.add_curve( {STH + EPS, 0.9},  [sigma](double s){ return std::real(sigma->direct_projection(0, s)); }, "Real");
    p2.add_curve( {STH + EPS, 0.9},  [sigma](double s){ return std::imag(sigma->direct_projection(0, s)); }, "Imaginary");
    p2.add_curve( {STH + EPS, 0.9}, [](double s){ return std::real(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey, "GKPY"));
    p2.add_curve( {STH + EPS, 0.9}, [](double s){ return std::imag(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey));

    plotter.combine({2,1}, {p2, p1}, "a00.pdf");
};