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
                std::complex<double> th = iso->direct_projection(1, s);

                dos += std::norm( std::imag(th) - std::imag(ex) );
                dos += std::norm( std::real(th) - std::real(ex) );
            };
        };
        return dos; 
    };
};

void fit()
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
    isobar f1 = new_isobar<truncated>(iso, 5, alpha, "truncated, n = 5");
    
    data_set pipi_pwave = pipi::partial_wave(iso, J, 10, {0.1, 0.9});

    fitter<pipi_fit> fitter(f1, alpha);
    fitter.set_parameter_labels({"lam2 (iso)", "g (iso)", "lam2", "alpha(0)", "g", "gamma", "c"});
    fitter.add_data( pipi_pwave );

    // Sync isobar's parameters to the trajectory as required by unitarity
    fitter.sync_parameter("g (iso)",    "g");
    fitter.sync_parameter("lam2 (iso)", "lam2");
    
    fitter.fix_parameter("lam2",     2.0);
    fitter.fix_parameter("alpha(0)", 0.491);
    fitter.set_parameter_posdef("g");
    fitter.set_parameter_posdef("gamma");
    fitter.set_parameter_posdef("c");

    fitter.do_iterative_fit({3.16, 1.2, 0.1}, 20, "rho_fit_2GeV");

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    data_set rhos = rho_spectrum();
    data_set as   = a_spectrum();

    plot p1 = plotter.new_plot();
    p1.add_logo(false);
    p1.set_labels("#it{s}  [GeV^{2}]", "#alpha_{#rho}(#it{s})");
    p1.set_ranges({-0.2, 1.5}, {0., 2.0});
    p1.set_legend(0.4, 0.7);
    auto reAlp = p1.add_curve(  {-0.2, 1.5},      [alpha](double s){ return alpha->real_part(s);} ,      "Real");
    auto imAlp = p1.add_curve(  {-0.2, 1.5},      [alpha](double s){ return alpha->imaginary_part(s); }, "Imaginary");
    auto Lin   = p1.add_curve(  {-0.2, 1.5},      [alpha](double s){ return 0.5+0.9*s; }, dashed(jpacColor::DarkGrey, "0.5 + 0.9 #it{s}"));
    p1.add_data({square_elementwise(rhos._x), {}}, {rhos._y, {}}, jpacColor::DarkGrey);
    p1.add_vertical(  0, {kBlack, kSolid});
    p1.shade_region({STH,1});
    p1.save("timelike.pdf");

    plot p2 = plotter.new_plot();
    p2.add_logo(false);
    p2.set_labels("#it{s}  [GeV^{2}]", "#it{f}_{1}^{1}(#it{s})");
    p2.set_ranges({0, 1}, {-0.7, 1.3});
    p2.set_legend(0.25, 0.7);
    p2.color_offset(2);
    auto reF00  = p2.add_curve(  {0,1}, [f1](double s){ return std::real(f1->direct_projection(1, s));} , "Real");
    auto imF00  = p2.add_curve(  {0,1}, [f1](double s){ return std::imag(f1->direct_projection(1, s)); }, "Imaginary");
    auto imUni  = p2.add_curve(  {0,1}, [f1](double s){ return (s > STH) ? sqrt(1.- STH/s) * std::norm(f1->direct_projection(1,s)) : 0; }, dashed(jpacColor::Orange, "Exact Unitarity"));
    auto reGKPY = p2.add_curve(  {0,1}, []  (double s){ return (s > STH) ? std::real(pipi::partial_wave(1, 1, s)) : 0.;}, dashed(jpacColor::DarkGrey, "GKPY"));
    auto imGKPY = p2.add_curve(  {0,1}, []  (double s){ return (s > STH) ? std::imag(pipi::partial_wave(1, 1, s)) : 0.;}, dashed(jpacColor::DarkGrey));
    p2.save("f11.pdf");

    auto dat = import_data<4>("data/charge_exchange.dat");
    plot p3 = plotter.new_plot();
    p3.add_logo(false);
    p3.set_labels("#it{s}  [GeV^{2}]", "Re #alpha_{#rho}(#it{s})");
    p3.set_ranges( {-1.5, 0}, {-0.7, 0.6});
    auto reTL = p3.add_curve(  {-1.5, 0},      [alpha](double s){ return alpha->real_part(s);});
    auto liTL = p3.add_curve(  {-1.5, 0},      [alpha](double s){ return 0.5+0.9*s; }, dashed(jpacColor::DarkGrey));
    p3.add_data({-dat[0], dat[1]/2}, {dat[2], dat[3]});
    p3.save("spacelike.pdf");

    // print_to_file<6>("fig5.txt", {"s", "ReF00", "ImF00", "Ex_Uni", "ReGKPY", "ImGKPY"}, {reF00[0], reF00[1], imF00[1], imUni[1], reGKPY[1], imGKPY[1]});
    // print_to_file<4>("fig6.txt", {"s", "ReAlpha", "ImAlpha", "Linear"}, {reAlp[0], reAlp[1], imAlp[1], Lin[1]});
    // print_to_file<3>("fig7.txt", {"s", "ReAlpha", "Linear"}, {reTL[0], reTL[1], liTL[1]});
};