// Test of how many terms in the power series of hypergeometric functions
// are needed to approximate the isobar in region of interest
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobars/truncated.hpp"
#include "isobars/hypergeometric.hpp"
#include "trajectories/iterative.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "print.hpp"
#include "plotter.hpp"

void truncation_test()
{
    using namespace analyticRT;
    using complex = std::complex<double>;
    
    std::array<double,2> plot_bounds = {EPS, 1.4};
    
    trajectory alpha = new_trajectory<iterative>("linear");
    alpha->max_iterations(5);
    alpha->set_parameters({0.8, 0.53, 1.04260251, 40.4});

    isobar pwave_1 = new_isobar<truncated>(1, 1, alpha, "rho");
    pwave_1->set_parameters({1., 1.});

    isobar pwave_3 = new_isobar<truncated>(1, 3, alpha, "rho");
    pwave_3->set_parameters({1., 1.});

    isobar pwave_5 = new_isobar<truncated>(1, 5, alpha, "rho");
    pwave_5->set_parameters({1., 1.});

    isobar hyp     = new_isobar<hypergeometric>(1, alpha, "hypergeo");
    hyp->set_parameters({1., 1.});

    plotter plotter;

    int j = 1;
    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#tilde{f}_{1}(s)");
    p1.set_legend(0.4, 0.7);
    p1.add_curve(  plot_bounds, [pwave_1, j](double s){ return std::real(pwave_1->cross_projection(j, s)); }, "n_{max} = 1");
    p1.add_curve(  plot_bounds, [pwave_3, j](double s){ return std::real(pwave_3->cross_projection(j, s)); }, "n_{max} = 3");
    p1.add_curve(  plot_bounds, [pwave_5, j](double s){ return std::real(pwave_5->cross_projection(j, s)); }, "n_{max} = 5");
    p1.add_curve(  plot_bounds, [hyp, j]    (double s){ return std::real(hyp->cross_projection(j, s));     }, "n_{max} = #infty");

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "f_{1}(s)");
    p2.set_legend(0.74, 0.6);
    p2.add_curve(  plot_bounds, [pwave_1, j](double s){ return std::real(pwave_1->direct_projection(j, s)); }, "n_{max} = 1");
    p2.add_dashed( plot_bounds, [pwave_1, j](double s){ return std::imag(pwave_1->direct_projection(j, s)); });
    p2.add_curve(  plot_bounds, [pwave_3, j](double s){ return std::real(pwave_3->direct_projection(j, s)); }, "n_{max} = 3");
    p2.add_dashed( plot_bounds, [pwave_3, j](double s){ return std::imag(pwave_3->direct_projection(j, s)); });
    p2.add_curve(  plot_bounds, [pwave_5, j](double s){ return std::real(pwave_5->direct_projection(j, s)); }, "n_{max} = 5");
    p2.add_dashed( plot_bounds, [pwave_5, j](double s){ return std::imag(pwave_5->direct_projection(j, s)); });
    p2.add_curve(  plot_bounds, [hyp, j]    (double s){ return std::real(hyp->direct_projection(j, s));     }, "n_{max} = #infty");
    p2.add_dashed( plot_bounds, [hyp, j]    (double s){ return std::imag(hyp->direct_projection(j, s));     });

    plotter.combine({2,1}, {p2, p1}, "truncated.pdf");
};