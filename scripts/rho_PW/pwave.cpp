// Comparison of the P-wave projection of hypergeometric isobar compared to pipi data
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobars/hypergeometric.hpp"
#include "trajectories/iterative.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "pipi.hpp"
#include "print.hpp"
#include "plotter.hpp"

void pwave()
{
    using namespace analyticRT;
    using complex = std::complex<double>;
    
    std::array<double,2> plot_bounds = {EPS, 1.4};

    double lambda2 = 0.8, g = 40.4;
    trajectory alpha = new_trajectory<iterative>("iterated #rho");
    alpha->max_iterations(5);
    alpha->set_parameters({lambda2, 0.53, 1.04260251, g});
    
    isobar rho = new_isobar<hypergeometric>(1, alpha, "rho");
    rho->set_parameters({3.*g, lambda2});

    amplitude amp;
    amp.add_isobar(rho);
    amp.ignore_cross(true);

    plotter plotter;
    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#alpha(s)");
    p1.set_legend(0.3, 0.7);
    p1.add_curve( plot_bounds, [alpha](double s){ return alpha->real_part(s);      }, "Real");
    p1.add_curve( plot_bounds, [alpha](double s){ return alpha->imaginary_part(s); }, "Imag");

    plot p2 = plotter.new_plot();
    p2.color_offset(2);
    p2.set_labels("#it{s}  [GeV^{2}]", "#it{f}_{1}(s)");
    p2.set_legend(0.7, 0.7);
    p2.add_curve(  plot_bounds, [&amp](double s){ return std::real(amp.partial_wave(1, 1, s)); }, "Real");
    p2.add_curve(  plot_bounds, [&amp](double s){ return std::imag(amp.partial_wave(1, 1, s)); }, "Imag");

    p2.add_dashed( {STH, plot_bounds[1]}, [rho](double s){ return phase_space(s)*std::norm(rho->direct_projection(1, s)); } );
    
    p2.add_curve(  {0.1, plot_bounds[1]}, [](double s){ return 16*PI*std::imag(pipi::partial_wave(1, 1, s)); }, "Data" );
    p2.add_dashed( {0.1, plot_bounds[1]}, [](double s){ return 16*PI*std::real(pipi::partial_wave(1, 1, s)); });

    plotter.combine({2,1}, {p1, p2}, "pwave.pdf");
};