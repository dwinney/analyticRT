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

struct swave_fit
{
    using complex = std::complex<double>;

    static std::string data_type(int i){ return  "GKPY PW"; };

    // Sum difference of squares for both real and imaginary part
    static double fcn(const std::vector<data_set> &data_sets, isobar f, trajectory alpha)
    { 
        double dos = 0; // difference of squares 
        for (auto &data : data_sets)
        {
            for (int i = 0; i < data._N; i++) dos += std::norm( f->direct_projection(0, data._x[i]) - complex(data._y[i], data._z[i]) );
        };
        return dos; 
    };
};

void fit()
{
    using namespace analyticRT;
    using complex = std::complex<double>;

    // Global constants
    int iso = 0, J = 0;

    // -------------------------------- ------------------------------------------
    // Set up the unitary dispersive trajectory

    trajectory alpha = new_trajectory<unitary>(iso, "sigma");
    alpha->set_option(unitary::kExpandAlpha);

    // The trajectory defines an isobar
    isobar f0 = new_isobar<truncated>(iso, 10, alpha, "I = 0 only");
    f0->set_option(truncated::kAddConstant);

    // --------------------------------------------------------------------------
    // If fitting doing a fit uncomment this

    fitter<swave_fit> fitter(f0, alpha);
    fitter.set_parameter_labels({"lam2 (iso)", "g (iso)", "gp (iso)", "lam2", "alpha(sth)", "g", "gp", "gamma", "c"});
    fitter.add_data( pipi::partial_wave(iso, J,  10, {0.2, 0.60}) );
    fitter.add_data( pipi::partial_wave(iso, J,  10, {STH, 0.1}) );

    // Sync isobar's parameters to the trajectory as required by unitarity
    fitter.sync_parameter("g (iso)",    "g");
    fitter.sync_parameter("lam2 (iso)", "lam2");
    fitter.sync_parameter("gp (iso)",  "gp");

    fitter.fix_parameter("lam2",  2.0);

    fitter.set_parameter_posdef("gamma");
    fitter.set_parameter_posdef("g");
    fitter.set_parameter_posdef("c");

    // k = 2 and up to 0.40
    fitter.do_fit({ -0.056164337, 0.21948343, 3.6842419, 0.0069999613, 41.747148});

    print("alpha(0) = ", alpha->real_part(0));

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#it{f}^{0}_{0}(#it{s})");
    p1.color_offset(2);
    p1.set_legend(0.25, 0.7);
    p1.set_ranges({0, 0.9}, {-0.3, 1.7});
    p1.add_curve( {0, 0.9},  [f0]( double s){ return std::real(f0->direct_projection(0, s)); }, "Real");
    p1.add_curve( {0, 0.9},  [f0]( double s){ return std::imag(f0->direct_projection(0, s)); }, "Imaginary");
    p1.add_curve( {0, 0.9},  [f0]( double s){ return (s > STH) ? sqrt(1.- STH/s) * std::norm(f0->direct_projection(0,s)) : 0; }, dashed(jpacColor::Orange, "Exact Unitarity"));
    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::real(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey, "GKPY"));
    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::imag(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey        ));
    p1.save("a00_PW.pdf");

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#alpha_{#sigma}(#it{s})");
    p2.set_ranges({-0.2, 1.5}, {-0.08, 0.04});
    p2.set_legend(0.65, 0.4);
    p2.shade_region({STH, 0.6});
    p2.add_vertical(  0, {kBlack, kSolid});
    p2.add_horizontal(0, {kBlack, kSolid}); 
    p2.add_curve( {-0.2,  1.5}, [alpha](double s){ return alpha->real_part(s); },              "Real");
    p2.add_curve( {-0.2,  1.5}, [alpha](double s){ return alpha->imaginary_part(s); },         "Imaginary");
    p2.add_data({std::vector<double>({0.5*0.5}),{}}, {std::vector<double>({0.}), {}}, jpacColor::DarkGrey);
    p2.save("alpha_sigma.pdf");

    plotter.combine({2,1}, {p2,p1}, "a00_results.pdf");
};