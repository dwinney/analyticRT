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
    auto guess = [](double s){ return (-0.4 + 0.01*s)/sqrt(1+1.2*s); };

    // auto guess = [](double s){ return (-0.43 + 0.01*s)/sqrt(1+2*s); };

    trajectory alpha = new_trajectory<unitary>(iso, guess, "sigma");
    alpha->set_option(unitary::kAddConstant);
    alpha->set_integrator_depth(10);
    iterable(alpha)->set_interp_pars(400, {50, 200});

    // The trajectory defines an isobar
    isobar f0 = new_isobar<truncated>(iso, 0, alpha, "I = 0 only");
    f0->set_option(truncated::kAddConstant);

    // -------------------------------- ------------------------------------------
    // Instead of fitting to the partial waves we fit to the trajectory

    // Iterate once to begin fitting
    std::vector<std::vector<double>> pars; 
 

    pars.push_back({1.5, -0.467211,  9.30962144, 0.0218357688, 0, 19.5350503});

    // {-0.46, -0.2} (s/sh)
    // pars.push_back({1.5, -0.4767211,  9.30962144, 0.0218357688, 0, 19.5350503});
    // pars.push_back({1.5, -0.47383736, 13.9306818, 0.0338991508, 0, 29.4      });

    // --------------------------------------------------------------------------
    // If fitting doing a fit uncomment this

    fitter<swave_fit> fitter(f0, alpha);
    fitter.set_parameter_labels({"lam2 (iso)", "g (iso)", "gp (iso)", "lam2", "alpha(0)", "g", "gamma", "c", "gp"});
    // fitter.add_data( pipi::partial_wave(iso, J,  5, {0.3, 0.5}) );
    fitter.add_data( pipi::partial_wave(iso, J,  10, {STH, 0.50}) );
    fitter.add_data( pipi::partial_wave(iso, J,  2, {STH, 0.15}) );

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
    std::vector<std::vector<double>> fitpars; std::vector<double> ds;

    fitter.do_fit({initpars[1], initpars[2], initpars[3], initpars[5]});
    // fitpars.push_back(fitter.pars()); ds.push_back(fitter.fcn_dof());

    // alpha->iterate();
    // fitter.do_fit({initpars[1], initpars[2], initpars[3], initpars[5]});
    // fitpars.push_back(fitter.pars()); ds.push_back(fitter.fcn_dof());


    // alpha->iterate();
    // fitter.fix_parameter("gp", 29.4);
    // fitter.do_fit({initpars[1], initpars[2], initpars[3]});
    // fitpars.push_back(fitter.pars()); ds.push_back(fitter.fcn_dof());
    
    // divider(6); centered(6, "FIT RESULTS"); divider(6);
    // print("i", "alpha(0)", "g_sig", "gamma", "g_pom", "d^2");
    // divider(6);
    // for (int i = 0; i < fitpars.size(); i++)
    // {
    //     print(i, fitpars[i][4], fitpars[i][5], fitpars[i][6], fitpars[i][8], ds[i]);
    // };
    
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
    p1.save("a00_PW.pdf");

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#alpha(#it{s})");
    p2.set_legend(0.65, 0.2);
    // p2.print_to_terminal(true);
    p2.add_curve( {-0.5,  2}, [alpha](double s){ return alpha->real_part(s); },              "Real");
    p2.add_curve( {-0.5,  2}, [alpha](double s){ return alpha->imaginary_part(s); },         "Imaginary");
    p2.add_curve( {STH+EPS,  2},  [alpha](double s){ return iterable(alpha)->previous_real(s); }, solid(jpacColor::DarkGrey, "Previous iteration"));
    p2.add_curve( {STH+EPS,  2},  [alpha,guess](double s){ return guess(s); }, dashed(jpacColor::Purple, "Previous iteration"));

    plotter.combine({2,1}, {p2,p1}, "a00_results.pdf");

    // plot p3 = plotter.new_plot();
    // p3.set_labels("#it{s}  [GeV^{2}]", "#alpha(#it{s})");
    // p3.set_legend(0.65, 0.2);
    // p3.set_logscale(true, false);
    // p3.add_curve( {1, 10000}, [alpha](double s){ return alpha->real_part(s);      },  "Real");
    // p3.add_curve( {1, 10000}, [alpha](double s){ return alpha->imaginary_part(s); },  "Imaginary");

    // plot p4 = plotter.new_plot();
    // p4.set_labels("#it{s}  [GeV^{2}]", "#alpha(#it{s})");
    // p4.set_legend(0.65, 0.2);
    // p4.add_curve(  {50, 500}, [alpha](double s){ return alpha->real_part(s);      },   "Real");
    // p4.add_curve(  {50, 500}, [alpha](double s){ return alpha->imaginary_part(s); },   "Imaginary");

    // plotter.combine({2,2}, {p2, p1, p3, p4}, "a00_checks.pdf");
};