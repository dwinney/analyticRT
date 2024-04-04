// Comparison of the P-wave projection of hypergeometric isobar and pipi data
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobars/truncated.hpp"
#include "trajectories/unitary.hpp"

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

                dos += std::norm( std::imag(th) - 16*PI*std::imag(ex) );
            };
        };
        return dos; 
    };
};

void pwave()
{
    using namespace analyticRT;

    // ---------------------------------------------------------------------------
    // For the I = 1 we can assume for now that there is a single trajectory

    // Given by the dispersive form
    trajectory alpha = new_trajectory<unitary>(4.*M2_PION, 1, "#rho");

    // The trajectory defines an isobar
    isobar rho = new_isobar<truncated>(1, 5, alpha, "I = 1");

    data_set pipi_pwave = pipi::partial_wave(1, 1, 5, {0., 0.6});

    fitter<pipi_fit> fitter(rho, alpha);
    fitter.set_parameter_labels({"g (iso)", "lam2 (iso)", "alpha(0)", "lam2", "g", "gamma"});
    fitter.add_data( pipi_pwave );
    fitter.set_guess_range({0, 200});

    // Sync isobar's parameters to the trajectory as required by unitarity
    fitter.sync_parameter("g (iso)", "g");
    fitter.sync_parameter("lam2 (iso)", "lam2");
    
    fitter.fix_parameter("lam2", 1.5);
    fitter.set_parameter_posdef("g");
    fitter.set_parameter_posdef("gamma");
    fitter.set_parameter_limits("alpha(0)", {0.3, 1});

    fitter.do_iterative_fit({0.55, 165, 1.2}, 1, "rhofit");

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    // Trajectory vs rho masses
    plot p1 = plotter.new_plot();
    p1.set_legend_spacing(0.03);
    p1.set_legend(0.25, 0.7);
    p1.set_curve_points(50);
    p1.add_curve({0, 1}, [alpha](double s){ return alpha->real_part(s); },      "Real");
    p1.add_curve({0, 1}, [alpha](double s){ return alpha->imaginary_part(s); }, "Imaginary");

    // P-wave real part vs data
    entry_style data_style;
    data_style._color = jpacColor::DarkGrey;
    data_style._style = kDashed;
    data_style._label = "Data";

    plot p3 = plotter.new_plot();
    p3.set_labels("#it{s}  [GeV^{2}]", "Re _{1}(s)");
    p3.set_legend(0.64, 0.20);
    p3.add_curve(  {0,1}, [rho](double s){ return std::real(rho->direct_projection(1, s)); }, rho->id());
    p3.add_curve( {STH + EPS, 1}, [](double s){ return 16*PI*std::real(pipi::partial_wave(1, 1, s)); }, data_style);

    plot p4 = plotter.new_plot();
    p4.color_offset(1);
    p4.set_labels("#it{s}  [GeV^{2}]", "Im f_{1}(s)");
    p4.set_legend(0.64, 0.65);
    p4.add_curve(  {0,1}, [rho](double s){ return std::imag(rho->direct_projection(1, s)); }, rho->id());
    p4.add_curve( {STH + EPS, 1}, [](double s){ return 16*PI*std::imag(pipi::partial_wave(1, 1, s)); }, data_style);
    plotter.combine({3,1}, {p1, p3, p4}, "pwave.pdf");
};