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

    trajectory alpha_K = kmatrix->get_trajectory();

    // -------------------------------- ------------------------------------------
    // Set up the unitary dispersive trajectory

    auto initial_sigma = [](double s){ return (-0.7 + 1.0*s) / sqrt(1 + s/20); };

    trajectory alpha = new_trajectory<unitary>(iso, initial_sigma, "Dispersive");

    // The trajectory defines an isobar
    isobar sigma = new_isobar<truncated>(iso, 4, alpha, "#it{I} = 0");
    sigma->set_option(truncated::kAddAdlerZero);

    // "data" to fit against
    data_set pipi_pwave = pipi::partial_wave(iso, J, 10, {0.1, 0.7});

    fitter<pipi_fit> fitter(sigma, alpha);
    fitter.set_parameter_labels({"lam2 (iso)", "g (iso)", "g_A", "m_sigma", "lam2", "alpha(0)", "g", "gamma", "c"});
    fitter.add_data( pipi_pwave );

    // Sync isobar's parameters to the trajectory as required by unitarity
    fitter.sync_parameter("g (iso)",    "g");
    fitter.fix_parameter( "g",    1.0);
    fitter.sync_parameter("lam2 (iso)", "lam2");
    fitter.fix_parameter( "lam2", 1.5);
    fitter.fix_parameter("m_sigma", 0.820);
    
    // Fit the remaining parameters g_A, alpha(sth), gamma, and c
    fitter.set_parameter_posdef("gamma");
    fitter.set_parameter_posdef("c");
    fitter.do_fit({0.36539468, -0.33808994, 1.2604883, 127.01431});

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#it{A}^{(0)}_{0}(s)");
    p1.set_legend(0.25, 0.2);
    p1.set_ranges({STH, 0.9}, {-0.3, 1.3});
    p1.add_curve( {STH + EPS, 0.9},  [sigma](double s){ return std::real(sigma->direct_projection(0, s)); }, "Real");
    p1.add_curve( {STH + EPS, 0.9},  [sigma](double s){ return std::imag(sigma->direct_projection(0, s)); }, "Imaginary");
    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::real(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey, "GKPY"));
    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::imag(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey));
    p1.add_curve( {STH + EPS, 0.9}, [kmatrix](double s){ return std::real(kmatrix->direct_projection(0, s));}, dashed(jpacColor::Green, "K-matrix"));
    p1.add_curve( {STH + EPS, 0.9}, [kmatrix](double s){ return std::imag(kmatrix->direct_projection(0, s));}, dashed(jpacColor::Green));

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#alpha^{#sigma}_{#it{s}}");
    p2.set_legend(0.70, 0.2);
    p2.set_ranges({0, 1}, {-0.75, 1.5});
    p2.add_curve( {0, 1},  [alpha](double s){ return alpha->real_part(s); }, "Real");
    p2.add_curve( {0, 1},  [alpha](double s){ return alpha->imaginary_part(s); }, "Imaginary");
    p2.add_curve( {EPS, 1}, [alpha_K](double s){ return alpha_K->real_part(s);}, dashed(jpacColor::Green, "K-matrix"));
    p2.add_curve( {EPS, 1}, [alpha_K](double s){ return alpha_K->imaginary_part(s);}, dashed(jpacColor::Green));

    plotter.combine({2,1}, {p1, p2}, "a00_dispersive.pdf");
};