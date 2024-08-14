// Comparison of the S-wave projection of hypergeometric isobar and pipi data
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "trajectories/unitary.hpp"

#include <sstream>
#include "kinematics.hpp"
#include "pipi.hpp"
#include "trajectory.hpp"
#include "isobar.hpp"
#include "plotter.hpp"
#include "data_set.hpp"

using namespace analyticRT;

void yukawa()
{
    using namespace analyticRT;
    using complex = std::complex<double>;

    // -------------------------------- ------------------------------------------
    // Set up the unitary dispersive trajectories

    trajectory sigma = new_trajectory<unitary>(0, "sigma");
    sigma->set_option(unitary::kExpandAlpha);
    sigma->set_parameters({2., -0.056164348, 0.2194834, 3.6842423, 0.0069999547, 41.746907});

    trajectory rho = new_trajectory<unitary>(1, [](double s){ return (0.5 + 0.9*s)/sqrt(1. + s/20.); }, "#rho");
    rho->set_integrator_depth(20);
    
    std::string dir = "scripts/rho/";
    std::string file_prefix = "rho_fit_2GeV_";
    auto traj_pars = import_transposed<21>(dir + file_prefix + "traj_pars.txt");
    iterable(rho)->iterate<21>(traj_pars, 20);

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    int N = 100000;
    double smin = -1000, smax = 5000;
    std::vector<double> s, im_s, re_s;
    for (int i = 0; i < N; i++)
    {
        double si = smin + double(i)*(smax - smin)/double(N-1);

        complex sigmai = sigma->evaluate(si);
        re_s.push_back(std::real(sigmai));
        im_s.push_back(std::imag(sigmai));
    };

    plot p1 = plotter.new_plot();
    p1.color_offset(1);
    p1.set_ranges({-0.25, -0.04}, {-0.01, 0.1});
    p1.set_labels("Re #alpha_{#sigma}(#it{s})", "Im #alpha_{#sigma}(#it{s})");
    p1.add_curve(re_s, im_s);
    p1.save("sigma.pdf");

    N = 100;
    smin = -5, smax = 20;
    std::vector<double> im_r, re_r;
    for (int i = 0; i < N; i++)
    {
        double si = smin + double(i)*(smax - smin)/double(N-1);

        s.push_back(si);
        re_r.push_back(rho->real_part(si));
        im_r.push_back(rho->imaginary_part(si));
    };

    plot p2 = plotter.new_plot();
    p2.set_ranges({-4, 17}, {-1, 7});
    p2.set_labels("Re #alpha_{#rho}(#it{s})", "Im #alpha_{#rho}(#it{s})");
    p2.add_curve(re_r, im_r);
    p2.save("rho.pdf");

};