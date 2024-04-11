// Comparison of the P-wave projection of hypergeometric isobar and pipi data
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

void pwave()
{
    using namespace analyticRT;

    // --------------------------------------------------------------------------
    // For the I = 1 we can assume for now that there is a single trajectory

    // Given by the dispersive form
    trajectory alpha = new_trajectory<unitary>(4.*M2_PION, 1, "#rho");
    iterable(alpha)->set_interp_pars(200, {50., 200.});

    // The trajectory defines an isobar
    isobar rho = new_isobar<truncated>(1, 5, alpha, "I = 1");

    data_set pipi_pwave = pipi::partial_wave(1, 1, 10, {0.1, 1.0});

    fitter<pipi_fit> fitter(rho, alpha);
    fitter.set_parameter_labels({"g (iso)", "lam2 (iso)", "alpha(0)", "lam2", "g", "gamma" , "c"});
    fitter.add_data( pipi_pwave );
    fitter.set_guess_range({0, 200});

    // Sync isobar's parameters to the trajectory as required by unitarity
    fitter.sync_parameter("g (iso)", "g");
    fitter.sync_parameter("lam2 (iso)", "lam2");
    
    fitter.fix_parameter("lam2", 1.5);
    fitter.fix_parameter("alpha(0)", 0.5);
    fitter.set_parameter_posdef("g");
    fitter.set_parameter_posdef("gamma");
    fitter.set_parameter_posdef("c");
    
    fitter.do_iterative_fit({3.16, 1.2, 0.4}, 20);

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    // P-wave real part vs data
    entry_style dashed;
    dashed._color = jpacColor::DarkGrey;
    dashed._style = kDashed;
    dashed._add_to_legend = true;

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#it{#alpha}^{#it{#rho}}_{#it{s}}");
    p1.set_ranges({-1, 4}, {-0.5, 5});
    p1.set_legend(0.25, 0.7);
    p1.add_curve(  {-1,4},      [alpha](double s){ return alpha->real_part(s);} , "Real");
    p1.add_curve(  {-1,4},      [alpha](double s){ return alpha->imaginary_part(s); }, "Imaginary");

    dashed._label = "0.5 + 0.9 #it{s}";
    p1.add_curve(  {-1,4},      [](double s){ return 0.5 + 0.9*s; }, dashed);
    p1.add_vertical((4*1.5 + STH)/2.);
    p1.save("rho_alpha.pdf");

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#it{A}_{1}^{(1)}(#it{s})");
    p2.set_ranges({0,1}, {-0.7, 1.3});
    p2.set_legend(0.25, 0.7);
    p2.color_offset(2);
    p2.add_curve(  {0,1},      [rho](double s){ return std::real(rho->direct_projection(1, s));} , "Real");
    p2.add_curve(  {0,1},      [rho](double s){ return std::imag(rho->direct_projection(1, s)); }, "Imaginary");
    dashed._label = "GKPY";
    p2.add_curve( {STH + EPS, 1}, [](double s){ return std::real(pipi::partial_wave(1, 1, s));}, dashed);
    dashed._add_to_legend = false;
    p2.add_curve( {STH + EPS, 1}, [](double s){ return std::imag(pipi::partial_wave(1, 1, s));}, dashed);
    p2.save("rho_PW.pdf");
};