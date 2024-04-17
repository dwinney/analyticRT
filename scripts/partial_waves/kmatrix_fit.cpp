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

    isobar A00 = new_isobar<k_matrix>(1, "K-matrix");

    data_set pipi_pwave = pipi::partial_wave(0, 0, 10, {0.1, 0.8});

    fitter<pipi_fit> fitter(A00, nullptr);
    fitter.add_data( pipi_pwave );
    fitter.set_guess_range({-10, 10});
    fitter.do_fit(100);

    std::vector<double> pars = fitter.pars();

    // ---------------------------------------------------------------------------
    // Make plot

    plotter plotter;

    // P-wave real part vs data
    entry_style dashed;
    dashed._color = jpacColor::DarkGrey;
    dashed._style = kDashed;
    dashed._add_to_legend = true;

    plot p1 = plotter.new_plot();
    p1.set_labels("#it{s}  [GeV^{2}]", "#it{A}_{0}^{(0)}(#it{s})");
    p1.set_ranges({0, 0.8}, {-0.3, 1.3});
    p1.set_legend(0.25, 0.7);
    p1.add_header("#it{K}-matrix");
    p1.add_curve(  {EPS, 0.8},      [A00](double s){ return std::real(A00->direct_projection(0, s));} , "Real");
    p1.add_curve(  {EPS, 0.8},      [A00](double s){ return std::imag(A00->direct_projection(0, s)); }, "Imaginary");
    dashed._label = "GKPY";
    p1.add_curve( {STH + EPS, 1}, [](double s){ return std::real(pipi::partial_wave(0, 0, s));}, dashed);
    dashed._add_to_legend = false;
    p1.add_curve( {STH + EPS, 1}, [](double s){ return std::imag(pipi::partial_wave(0, 0, s));}, dashed);

    plot p2 = plotter.new_plot();
    p2.set_labels("#it{s}  [GeV^{2}]", "#tilde{#alpha}^{#sigma}_{#it{s}}");
    p2.set_ranges({0, 0.8}, {-5, 2.5});
    p2.set_legend(0.7, 0.25);
    p2.add_header("#it{K}-matrix");
    p2.add_curve(  {STH + EPS, 0.8},      [A00](double s){ return std::real(-1/A00->direct_projection(0, s));} , "Real");
    p2.add_curve(  {STH + EPS, 0.8},      [A00](double s){ return std::imag(-1/A00->direct_projection(0, s)); }, "Imaginary");
    p2.add_vertical(STH);

    plotter.combine({2,1}, {p1, p2}, "a00_kmatrix.pdf");

    plot p3 = plotter.new_plot();
    p3.set_labels("#it{s}  [GeV^{2}]", "#alpha^{#sigma}_{#it{s}}");
    p3.set_ranges({0, 0.8}, {-0.5, 1.3});
    p3.set_legend(0.7, 0.2);
    p3.add_header("#it{K}-matrix");

    A00->set_parameters({0, pars[1], pars[2]});
    double s_sig = 0.34766819735151394;
    double sub = std::real(-1/A00->direct_projection(0, s_sig));
    
    p3.add_curve(  {EPS, 0.8}, [A00,sub](double s){ return std::real(-1/A00->direct_projection(0, s)) - sub;}, "Real");
    p3.add_curve(  {EPS, 0.8}, [A00](double s){ return std::imag(-1/A00->direct_projection(0, s));}, "Imaginary");
    p3.add_vertical(STH);
    p3.save("alpha_kmatrix.pdf");
};