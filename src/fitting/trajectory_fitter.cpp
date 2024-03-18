// Class for taking in mass/width data and performing a fit a given trajectory model
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "trajectory_fitter.hpp"

namespace analyticRT
{
    // -----------------------------------------------------------------------
    // Set limits, labels, and fix parameters

    void trajectory_fitter::reset_parameters()
    {
        _pars.clear();
        _Nfree = _trajectory->Nfree();

        // populate parameters vector of appropriate size
        for (int i = 0; i < _trajectory->Nfree(); i++)
        {
            _pars.push_back(i);
        };
        set_parameter_labels(_trajectory->parameter_labels());
    };


    // -----------------------------------------------------------------------
    // Methods for managing data

    void trajectory_fitter::add_data(data_set data)
    {
        switch (data._type)
        {
            case regge_trajectory:
            {
                _data.push_back(data); _N += data._N;
                break;
            };
            default:
            {
                warning("fitter::add_data", "data_set " + data._id + " of unsupported type!");
                return;
            }
        };
    };

    void trajectory_fitter::clear_data()
    {
        _N = 0; _data.clear();
    };

    // ---------------------------------------------------------------------------
    // Chi-squared functions

    // Total chi2, this combines all data sets and is the function that gets called 
    // by the minimizer
    double trajectory_fitter::fit_chi2(const double * cpars)
    {
        // First convert the C string to a C++ vector
        std::vector<double> pars = convert(cpars);

        // Pass parameters to the amplitude
        _trajectory->set_parameters(pars);

        // Then sum over all data sets and observables
        double totchi2 = 0;
        for (auto data : _data)
        {
            totchi2 += chi2(data);
        };

        return totchi2;
    };

    // Calculate the chi2 from integrated xsection data
    double trajectory_fitter::chi2(data_set data)
    {
        // Sum over data points
        double chi2 = 0;
        for (int i = 0; i < data._N; i++)
        {
            double s = data._x[i];

            complex alpha_th  = _trajectory->evaluate(s);
            if (_fit_imag) chi2 += std::norm(alpha_th - (data._y[i] + I*data._z[i]));
            else           chi2 += std::norm(std::real(alpha_th) - data._y[i]);
        };

        return chi2;
    };

    // ---------------------------------------------------------------------------
    // Status messages printed to command line

    // Summary of data sets that have been recieved
    void trajectory_fitter::data_info()
    {
        using std::cout; 
        using std::left;
        using std::setw;
        using std::endl;

        cout << left;
        divider();
        cout << "Fitting trajectory (\"" << _trajectory->id() << "\") to " << _N << " data points:" << endl;
        line();
        cout << setw(30) << "DATA SET"         << setw(20) << "POINTS" << endl;
        cout << setw(30) << "----------------" << setw(20) << "-------" << endl;

        for (auto data : _data)
        {   
            cout << setw(30) << data._id  << setw(20) << data._x.size() << endl;  
        };
    };
};