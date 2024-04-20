// Fit the S-wave pi pi phase shifts to a simple K-matrix
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "isobars/k_matrix.hpp"

#include "kinematics.hpp"
#include "isobar.hpp"
#include "plotter.hpp"
#include "fitter.hpp"
#include "pipi.hpp"

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
                std::complex<double> th = iso->direct_projection(0, s);

                dos += std::norm( std::imag(th) - std::imag(ex) );
                dos += std::norm( std::real(th) - std::real(ex) );
            };
        };
        return dos; 
    };
};

void kmatrix_fit()
{
    using namespace analyticRT;

    // --------------------------------------------------------------------------
    // For the I = 1 we can assume for now that there is a single trajectory

    isobar   A00        = new_isobar<k_matrix>(0, "K-matrix");
    data_set pipi_pwave = pipi::partial_wave(0, 0, 10, {0.1, 0.8});

    fitter<pipi_fit> fitter(A00, nullptr); // fitter takes an iosbar but not a trajectory
    fitter.add_data( pipi_pwave );
    fitter.set_guess_range({-10, 10});
    fitter.do_fit(20);

    // After the fit, extract the trajectory
    trajectory alpha = A00->get_trajectory();

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#it{A}_{0}^{(0)}(#it{s})");
    p1.set_legend(0.25, 0.7);
    p1.set_ranges({0, 0.8}, {-0.3, 1.3});
    p1.add_header("#it{K}-matrix");
    p1.add_curve(  {EPS, 0.8},      [A00](double s){ return std::real(A00->direct_projection(0, s));} , "Real");
    p1.add_curve(  {EPS, 0.8},      [A00](double s){ return std::imag(A00->direct_projection(0, s)); }, "Imaginary");
    p1.add_curve( {STH + EPS, 1}, [](double s){ return std::real(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey, "GKPY"));
    p1.add_curve( {STH + EPS, 1}, [](double s){ return std::imag(pipi::partial_wave(0, 0, s));}, dashed(jpacColor::DarkGrey));

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#tilde{#alpha}^{#sigma}_{#it{s}} / #it{g}_{#sigma}");
    p2.set_ranges({0, 0.8}, {-5, 2.5});
    p2.set_legend(0.7, 0.25);
    p2.add_header("#it{K}-matrix");
    p2.add_curve(  {EPS, 0.8},      [A00](double s){ return std::real(-1./A00->direct_projection(0, s));} , "Real");
    p2.add_curve(  {EPS, 0.8},      [A00](double s){ return std::imag(-1./A00->direct_projection(0, s)); }, "Imaginary");

    plotter.combine({2,1}, {p1, p2}, "a00_kmatrix.pdf");

    plot p3 = plotter.new_plot();
    p3.set_labels("#it{s}  [GeV^{2}]", "#alpha^{#sigma}_{#it{s}} / #it{g}_{#sigma}");
    p3.set_ranges({0, 1.2}, {-0.5, 1.5});
    p3.set_legend(0.25, 0.75);
    p3.add_header("#it{K}-matrix");

    p3.add_curve(  {EPS, 1.2}, [alpha](double s){ return alpha->real_part(s);},      "Real");
    p3.add_curve(  {EPS, 1.2}, [alpha](double s){ return alpha->imaginary_part(s);}, "Imaginary");
    p3.add_horizontal(0, {kBlack, kSolid});
    p3.save("alpha_kmatrix.pdf");
};