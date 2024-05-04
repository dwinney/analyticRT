// Comparison of the S-wave projection of hypergeometric isobar and pipi data
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobars/truncated.hpp"
#include "trajectories/unitary.hpp"

#include "kinematics.hpp"
#include "pipi.hpp"
#include "trajectory.hpp"
#include "isobar.hpp"
#include "plotter.hpp"
#include "fitter.hpp"
#include "data_set.hpp"

using namespace analyticRT;

struct swave_fit
{
    static std::string data_type(int i){ return (i == 0) ? "GKPY PW" : "Trajectory"; };


    // Sum difference of squares for both real and imaginary part
    static double fcn(const std::vector<data_set> &data_sets, isobar f, trajectory alpha)
    { 
        double dos = 0; // difference of squares 
        for (auto &data : data_sets)
        {
            for (int i = 0; i < data._N; i++)
            {
                double s     = data._x[i];
                std::complex<double> ex(data._y[i], data._z[i]); 
                std::complex<double> th = (data._type == 0) ? f->direct_projection(0, s) : alpha->evaluate(s);

                dos += std::norm( th - ex );
            };
        };

        // In addition to the data above, we add two "rubber band" points
        dos += std::norm(iterable(alpha)->previous_evaluate(25)  - alpha->evaluate(25) );
        dos += std::norm(iterable(alpha)->previous_evaluate(100) - alpha->evaluate(100));

        return dos; 
    };
};


void swave()
{
    using namespace analyticRT;
    using complex = std::complex<double>;

    // Global constants
    int iso = 0, J = 0;

    // -------------------------------- ------------------------------------------
    // Set up the unitary dispersive trajectory

    auto initial_sigma = [](double s){ return (-1. + 0.8*s) / sqrt(1 + s/20); };

    trajectory alpha = new_trajectory<unitary>(iso, initial_sigma, "Dispersive");
    alpha->set_option(unitary::kAddConstant);

    // The trajectory defines an isobar
    isobar sigma = new_isobar<truncated>(iso, 0, alpha, "truncated, n = 4");
    sigma->set_option(truncated::kAddConstant);

    // -------------------------------- ------------------------------------------
    // GKPY partial waves

    // "data" to fit against
    data_set pipi_swave = pipi::partial_wave(iso, J, 10, {0.1, 0.7});

    // -------------------------------- ------------------------------------------
    // Instead of fitting to the partial waves we fit to the trajectory

    // Iterate once to begin fitting
    std::vector<double> zeroth_pars  = {1.5, -0.83683862, 0.89805832, 0.64101313,  21.942810,  0.93737249};
    alpha->set_parameters(zeroth_pars);
    alpha->iterate();
  
    // --------------------------------------------------------------------------
    // If fitting doing a fit uncomment this

    fitter<swave_fit> fitter(sigma, alpha);
    fitter.set_parameter_labels({"lam2 (iso)", "g (iso)", "gp (iso)", "lam2", "alpha(0)", "g", "gamma", "c", "gp"});
    fitter.add_data( pipi_swave );

    // Sync isobar's parameters to the trajectory as required by unitarity
    fitter.sync_parameter("g (iso)",    "g");
    fitter.sync_parameter("lam2 (iso)", "lam2");
    fitter.sync_parameter("gp (iso)",  "gp");

    fitter.fix_parameter( "lam2",      1.5);

    fitter.set_parameter_limits("alpha(0)", {-1, 1});
    fitter.set_parameter_posdef("gamma");
    fitter.set_parameter_posdef("g");
    fitter.set_parameter_posdef("gp");
    fitter.set_parameter_posdef("c");

    fitter.do_iterative_fit({-0.75930719, 0.78648168, 1.06071710,  3.7111063,  0.89640185}, 5, "a00");

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#it{A}^{(0)}_{0}(#it{s})");
    p1.color_offset(2);
    p1.set_legend(0.25, 0.7);
    p1.set_ranges( {0, 0.9}, {-0.3, 1.3});
    p1.add_curve(  {0, 0.9},  [sigma]( double s){ return std::real(sigma->direct_projection(0, s)); }, "Real");
    p1.add_curve(  {0, 0.9},  [sigma]( double s){ return std::imag(sigma->direct_projection(0, s)); }, "Imaginary");
    p1.add_dashed( {0, 0.9},  [sigma]( double s){ return (s > STH) ? sqrt(1.- STH/s) * std::norm(sigma->direct_projection(0,s)) : 0; });
    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::real(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey, "GKPY"));
    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::imag(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey        ));

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#alpha(#it{s})");
    p2.set_legend(0.70, 0.2);
    p2.add_curve( {-0.5, 250}, [alpha](double s){ return alpha->real_part(s); },              "Real");
    p2.add_curve( {-0.5, 250}, [alpha](double s){ return alpha->imaginary_part(s); },         "Imaginary");
    p2.add_curve( {STH, 250},  [alpha](double s){ return iterable(alpha)->previous_real(s); }, dashed(jpacColor::DarkGrey));
    p2.add_curve( {STH, 250},  [alpha](double s){ return iterable(alpha)->previous_imag(s); }, dashed(jpacColor::DarkGrey));

    plotter.combine({2,1}, {p2, p1}, "a00.pdf");
};