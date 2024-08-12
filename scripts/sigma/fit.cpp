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

    // pars.push_back({2, -0.175403, 0.198617, 0.554793, 1.0784, 1.01672});
    // pars.push_back({2, -0.192996, 0.584351, 0.119025, 0.877235, 3});
    // pars.push_back({2, -0.144277, 0.326603, 0.0970915, 0.00222155, 2.1789});
    // pars.push_back({2, -0.132981, 0.518653, 0.0598452, 3.71065e-05, 3.83722});

    for (auto p : pars)
    {
        alpha->set_parameters(p);
        alpha->iterate();
    };

    // --------------------------------------------------------------------------
    // If fitting doing a fit uncomment this

    fitter<swave_fit> fitter(f0, alpha);
    fitter.set_parameter_labels({"lam2 (iso)", "g (iso)", "gp (iso)", "lam2", "alpha(0)", "g", "gamma", "c", "gp"});
    fitter.add_data( pipi::partial_wave(iso, J,  10, {0.1, 0.70}) );
    fitter.add_data( pipi::partial_wave(iso, J,  2,  {STH, 0.15}) );

    // Sync isobar's parameters to the trajectory as required by unitarity
    fitter.sync_parameter("g (iso)",    "g");
    fitter.sync_parameter("lam2 (iso)", "lam2");
    fitter.sync_parameter("gp (iso)",  "gp");

    fitter.fix_parameter("lam2",  2.0);
    fitter.set_parameter_limits("gp", {3,10});

    fitter.set_parameter_posdef("gamma");
    fitter.set_parameter_posdef("g");
    fitter.set_parameter_posdef("c");


    std::vector<double> initpars = pars.back();

    fitter.do_fit({initpars[1], initpars[2], initpars[3], initpars[4], initpars[5]});

    auto fit_pars = fitter.pars();
    std::stringstream ss;
    ss << "pars.push_back({";
    for (int i = 0; i < 6; i++)
    {
        ss << fit_pars[3+i];
        if (3+i < fit_pars.size()-1) ss << ", ";
    };
    ss << "});";
    print(ss.str());

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

    plotter.combine({2,1}, {p2,p1}, "a00_results.pdf");
};