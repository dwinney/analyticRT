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

        // In addition to the data above, we add "rubber band" points to help force convergence
        std::vector<double> rubber_band_points = {10};
        for (auto s : rubber_band_points){ dos += std::norm(iterable(alpha)->previous_real(s) - alpha->real_part(s) ); };

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

    trajectory alpha = new_trajectory<unitary>(iso, std::array<double,3>({-1.0, 0.8, 20}), "Guess = {-1, 0.8, 20}");
    alpha->set_option(unitary::kAddConstant);
    iterable(alpha)->set_interp_pars(50, {50, 200});

    // The trajectory defines an isobar
    isobar f0 = new_isobar<truncated>(iso, 0, alpha, "truncated, n = 4");
    f0->set_option(truncated::kAddConstant);

    // -------------------------------- ------------------------------------------
    // GKPY partial waves

    // "data" to fit against
    data_set pipi_swave = pipi::partial_wave(iso, J,  10, {0.1, 0.6});

    // -------------------------------- ------------------------------------------
    // Instead of fitting to the partial waves we fit to the trajectory

    // Iterate once to begin fitting
    std::vector<std::vector<double>> pars; 
    pars.push_back({1.5, -0.92508322,  1.1485475, 0.48401511, 301.17225, 1.1122341} );  // 0
    pars.push_back({1.5, -0.97842706,  1.0304921, 1.2021969,  36.921587,  0.86899737}); // 1
    // pars.push_back({1.5, -1,  1.0123121,   1.3919758,  26.812212,  0.82535992});  // 2
    // pars.push_back({1.5, -1,  0.96206883,  1.4974325,  26.570337,  0.76244214}); // 3
    // pars.push_back({1.5, -1,  0.94599787,  1.5687863,  24.144071,  0.74362752}); // 4
    // pars.push_back({1.5, -1,  0.95816777,  1.5447815,  23.826851,  0.75874368}); // 5
    // pars.push_back({1.5, -1,  0.96087587,  1.5362693,  23.684327,  0.76224477}); // 6
    // pars.push_back({1.5, -1,  0.96198883,  1.538402,   23.689577,  0.76331045}); // 7
    // pars.push_back({1.5, -1,  0.95898862,  1.5638811,  23.42519,   0.75898323}); // 8
    
    // IF FITTING
    for (auto p : pars) { alpha->set_parameters(p); alpha->iterate(); };

    // // IF JUST PLOTTING
    // for (int i = 0; i < pars.size(); i++) 
    // { 
    //     alpha->set_parameters(pars[i]); 
    //     if (i == pars.size() - 1) f0->set_parameters({pars[i][0], pars[i][2], pars[i][5]});
    //     else alpha->iterate(); 
    // };

    // // IF IMPORTING FROM FILE
    // auto apars = import_transposed<6>("/scripts/partial_waves/a00_traj_pars.txt");
    // auto fpars = import_transposed<6>("/scripts/partial_waves/a00_iso_pars.txt");

    // int i = 5;
    // iterable(alpha)->iterate<6>(apars, i);
    // f0->set_parameters(fpars[i]);

    // --------------------------------------------------------------------------
    // If fitting doing a fit uncomment this

    fitter<swave_fit> fitter(f0, alpha);
    fitter.set_parameter_labels({"lam2 (iso)", "g (iso)", "gp (iso)", "lam2", "alpha(0)", "g", "gamma", "c", "gp"});
    fitter.add_data( pipi_swave );

    // Sync isobar's parameters to the trajectory as required by unitarity
    fitter.sync_parameter("g (iso)",    "g");
    fitter.sync_parameter("lam2 (iso)", "lam2");
    fitter.sync_parameter("gp (iso)",  "gp");

    fitter.fix_parameter( "lam2", 1.5);
    fitter.fix_parameter( "alpha(0)", -1.);

    // fitter.set_parameter_limits("alpha(0)", {-1, 1});
    fitter.set_parameter_posdef("gamma");
    fitter.set_parameter_posdef("g");
    fitter.set_parameter_posdef("gp");
    fitter.set_parameter_posdef("c");

    std::vector<double> initpars = pars.back();
    // fitter.do_fit({initpars[2], initpars[3], initpars[4], initpars[5]});
    fitter.do_iterative_fit({initpars[2], initpars[3], initpars[4], initpars[5]}, 10, "a00");

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#it{A}^{(0)}_{0}(#it{s})");
    p1.color_offset(2);
    p1.set_legend(0.25, 0.7);
    p1.set_ranges( {0, 0.9}, {-0.3, 1.3});
    p1.add_curve(  {0, 0.9},  [f0]( double s){ return std::real(f0->direct_projection(0, s)); }, "Real");
    p1.add_curve(  {0, 0.9},  [f0]( double s){ return std::imag(f0->direct_projection(0, s)); }, "Imaginary");
    p1.add_dashed( {0, 0.9},  [f0]( double s){ return (s > STH) ? sqrt(1.- STH/s) * std::norm(f0->direct_projection(0,s)) : 0; });
    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::real(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey, "GKPY"));
    p1.add_curve( {STH + EPS, 0.9}, [](double s){ return std::imag(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey        ));

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#alpha(#it{s})");
    p2.set_legend(0.65, 0.2);
    p2.add_curve( {-0.5, 2}, [alpha](double s){ return alpha->real_part(s); },              "Real");
    p2.add_curve( {-0.5, 2}, [alpha](double s){ return alpha->imaginary_part(s); },         "Imaginary");
    p2.add_curve( {STH,  2},  [alpha](double s){ return iterable(alpha)->previous_real(s); }, dashed(jpacColor::DarkGrey, "Previous iteration"));
    p2.add_curve( {STH,  2},  [alpha](double s){ return iterable(alpha)->previous_imag(s); }, dashed(jpacColor::DarkGrey));

    plot p3 = plotter.new_plot();
    p3.set_labels("#it{s}  [GeV^{2}]", "#alpha(#it{s})");
    p3.set_legend(0.65, 0.2);
    p3.add_curve( {-0.5, 500}, [alpha](double s){ return alpha->real_part(s); },              "Real");
    p3.add_curve( {-0.5, 500}, [alpha](double s){ return alpha->imaginary_part(s); },         "Imaginary");
    p3.add_curve( {STH,  500},  [alpha](double s){ return iterable(alpha)->previous_real(s); }, dashed(jpacColor::DarkGrey, "Previous iteration"));
    p3.add_curve( {STH,  500},  [alpha](double s){ return iterable(alpha)->previous_imag(s); }, dashed(jpacColor::DarkGrey));

    plotter.combine({3,1}, {p2, p3, p1}, "a00.pdf");
};