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

    trajectory alpha = new_trajectory<unitary>(iso, std::array<double,3>({-0.3, -0.2, 1/20}), "sigma");
    alpha->set_option(unitary::kAddConstant);
    iterable(alpha)->set_interp_pars(200, {50, 200});

    // The trajectory defines an isobar
    isobar f0 = new_isobar<truncated>(iso, 0, alpha, "I = 0 only");
    f0->set_option(truncated::kAddConstant);

    // -------------------------------- ------------------------------------------
    // GKPY partial waves

    // "data" to fit against
    data_set pipi_swave    = pipi::partial_wave(iso, J,  10, {STH, 0.45});

    // -------------------------------- ------------------------------------------
    // Instead of fitting to the partial waves we fit to the trajectory

    // Iterate once to begin fitting
    std::vector<std::vector<double>> pars; 

    // {-0.15, -0.2, 1/20}
    pars.push_back({1.5, -0.167538494 ,  0.916527549,  0.0273027551,  0,  5.4430882   });  // 0 

    // --------------------------------------------------------------------------
    // If fitting doing a fit uncomment this

    fitter<swave_fit> fitter(f0, alpha);
    fitter.set_parameter_labels({"lam2 (iso)", "g (iso)", "gp (iso)", "lam2", "alpha(0)", "g", "gamma", "c", "gp"});
    fitter.add_data( pipi_swave );

    // Sync isobar's parameters to the trajectory as required by unitarity
    fitter.sync_parameter("g (iso)",    "g");
    fitter.sync_parameter("lam2 (iso)", "lam2");
    fitter.sync_parameter("gp (iso)",  "gp");
    fitter.fix_parameter("lam2",  1.5);
    fitter.fix_parameter("c",       0);

    fitter.set_parameter_posdef("gamma");
    fitter.set_parameter_posdef("g");
    fitter.set_parameter_posdef("c");
    fitter.set_parameter_posdef("gp");

    std::vector<double> initpars = pars.back();
    // for (auto par : pars){ alpha->set_parameters(par); alpha->iterate(); };
    
    // fitter.do_fit({initpars[1], initpars[2], initpars[3], initpars[5]});
    // alpha->iterate();
    // fitter.do_fit({initpars[1], initpars[2], initpars[3], initpars[4], initpars[5]});    

    // alpha->iterate();
    // fitter.do_fit({initpars[1], initpars[2], initpars[3],initpars[4], initpars[5]});

    fitter.do_iterative_fit({initpars[1], initpars[2], initpars[3], initpars[5]}, 1, "a00");
    

    // --------------------------------------------------------------------------

    // // IF JUST PLOTTING
    // for (int i = 0; i < pars.size(); i++) 
    // { 
    //     alpha->set_parameters(pars[i]); 
    //     if (i == pars.size() - 1) f0->set_parameters({pars[i][0], pars[i][2], pars[i][5]});
    //     else alpha->iterate(); 
    // };

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#it{A}^{(0)}_{0}(#it{s})");
    p1.color_offset(2);
    p1.set_legend(0.25, 0.7);
    p1.set_ranges({0, 0.9}, {-0.3, 1.3});
    // p1.print_to_terminal(true);
    p1.add_curve( {0, 0.9},  [f0]( double s){ return std::real(f0->direct_projection(0, s)); }, "Real");
    p1.add_curve( {0, 0.9},  [f0]( double s){ return std::imag(f0->direct_projection(0, s)); }, "Imaginary");
    p1.add_curve( {0, 0.9},  [f0]( double s){ return (s > STH) ? sqrt(1.- STH/s) * std::norm(f0->direct_projection(0,s)) : 0; }, dashed(jpacColor::Orange, "Exact Unitarity"));
    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::real(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey, "GKPY"));
    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::imag(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey        ));

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#alpha(#it{s})");
    p2.set_legend(0.65, 0.2);
    // p2.print_to_terminal(true);
    p2.add_curve( {-0.5,  2}, [alpha](double s){ return alpha->real_part(s); },              "Real");
    p2.add_curve( {-0.5,  2}, [alpha](double s){ return alpha->imaginary_part(s); },         "Imaginary");
    p2.add_curve( {STH+EPS,  2},  [alpha](double s){ return iterable(alpha)->previous_real(s); }, dashed(jpacColor::DarkGrey, "Previous iteration"));
    p2.add_curve( {STH+EPS,  2},  [alpha](double s){ return iterable(alpha)->previous_imag(s); }, dashed(jpacColor::DarkGrey));

    plotter.combine({2,1}, {p2,p1}, "a00_results.pdf");

    plot p3 = plotter.new_plot();
    p3.set_labels("#it{s}  [GeV^{2}]", "#alpha(#it{s})");
    p3.set_legend(0.65, 0.2);
    p3.set_logscale(true, false);
    p3.add_curve( {1, 10000}, [alpha](double s){ return alpha->real_part(s);      },  "Real");
    p3.add_curve( {1, 10000}, [alpha](double s){ return alpha->imaginary_part(s); },  "Imaginary");

    plot p4 = plotter.new_plot();
    p4.set_labels("#it{s}  [GeV^{2}]", "#alpha(#it{s})");
    p4.set_legend(0.65, 0.2);
    p4.add_curve(  {50, 500}, [alpha](double s){ return alpha->real_part(s);      },   "Real");
    p4.add_curve(  {50, 500}, [alpha](double s){ return alpha->imaginary_part(s); },   "Imaginary");

    plotter.combine({2,2}, {p2, p1, p3, p4}, "a00_checks.pdf");
};