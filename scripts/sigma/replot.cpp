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

void replot()
{
    using namespace analyticRT;
    using complex = std::complex<double>;

    // Global constants
    int iso = 0, J = 0;

    // -------------------------------- ------------------------------------------
    // Set up the unitary dispersive trajectory

    trajectory alpha = new_trajectory<unitary>(iso, "sigma");
    alpha->set_option(unitary::kExpandAlpha);
    alpha->set_parameters({2.0, -0.089, 0.34765838, 3.6826296, 0.011097201, 66.172156});

    std::string dir = "scripts/sigma/";
    auto jrp_curves = import_data<5>(dir + "alpha-sig-error-plb-central.dat");

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#alpha_{#sigma}(#it{s})");
    p2.set_ranges({-0.2, 1.5}, {-0.15, 0.23});
    p2.set_legend(0.75, 0.6);
    p2.shade_region({STH, 0.6});
    p2.add_vertical(  0, {kBlack, kSolid});
    p2.add_horizontal(0, {kBlack, kSolid}); 
    p2.add_curve( {-0.2,  1.5}, [alpha](double s){ return alpha->real_part(s); },              "Real");
    p2.add_dashed( jrp_curves[0], jrp_curves[1]);
    p2.add_curve( {-0.2,  1.5}, [alpha](double s){ return alpha->imaginary_part(s); },         "Imaginary");
    p2.add_dashed( jrp_curves[0], jrp_curves[2]);

    // p2.add_data({std::vector<double>({0.5*0.5}),{}}, {std::vector<double>({0.}), {}}, jpacColor::DarkGrey);
    p2.save("alpha_sigma.pdf");
};