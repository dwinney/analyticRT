// Class for taking in mass/width data and performing a fit a given trajectory model
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "spectrum_fitter.hpp"

namespace analyticRT
{
    // -----------------------------------------------------------------------
    // Set limits, labels, and fix parameters

    void spectrum_fitter::reset_parameters()
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

    void spectrum_fitter::add_data(data_set data)
    {
        switch (data._type)
        {
            case spectrum:
            {
                _spectrum_data.push_back(data); 
                (_ignore_widths) ?  _N += data._N : _N += 2*data._N;
                break;
            };
            case timelike: 
            {
                _timelike_data.push_back(data); 
                _N += data._N;
                break;
            };
            default:
            {
                warning("fitter::add_data", "data_set " + data._id + " of unsupported type!");
                return;
            }
        };
    };

    void spectrum_fitter::clear_data()
    {
        _N = 0;
        _spectrum_data.clear();
        _timelike_data.clear();
    };

    // ---------------------------------------------------------------------------
    // Chi-squared functions

    // Total chi2, this combines all data sets and is the function that gets called 
    // by the minimizer
    double spectrum_fitter::fit_chi2(const double * cpars)
    {
        // First convert the C string to a C++ vector
        std::vector<double> pars = convert(cpars);

        // Pass parameters to the amplitude
        _trajectory->set_parameters(pars);

        // Then sum over all data sets and observables
        double chi2 = 0;
        for (auto data : _spectrum_data)
        {
            chi2 += chi2_spectrum(data);
        };
        for (auto data : _timelike_data)
        {
            chi2 += chi2_timelike(data);
        };

        return chi2;
    };

    // Calculate the chi2 from integrated xsection data
    double spectrum_fitter::chi2_spectrum(data_set data)
    {
        // Sum over data points
        double chi2 = 0;
        for (int i = 0; i < data._N; i++)
        {
            double s = pow(data._x[i], 2);

            double spin_th  = _trajectory->real_part(s);
            double spin_ex  = data._y[i];
            double spin_err = data._dy[i];
            chi2 += pow((spin_th - spin_ex)/spin_err, 2);
            
            if (_ignore_widths) continue; 

            double width_th  = _trajectory->width(s);
            double width_ex  = data._z[i];
            double width_err = data._dz[i];
            chi2 += pow((width_th - width_ex)/width_err, 2);
        };

        return chi2;
    };

    // Calculate chi2 from differential data set
    double spectrum_fitter::chi2_timelike(data_set data)
    {   
        double chi2 = 0;
        for (int i = 0; i < data._N; i++)
        {
            double s = data._x[i];
            double alpha_th =  _trajectory->real_part(s);
            double alpha_ex = data._y[i];
            double error    = data._dy[i];

            chi2 += pow((alpha_th - alpha_ex)/error, 2);
        };

        return chi2;
    };
    
    // ---------------------------------------------------------------------------
    // Status messages printed to command line

    // Summary of data sets that have been recieved
    void spectrum_fitter::data_info()
    {
        using std::cout; 
        using std::left;
        using std::setw;
        using std::endl;

        cout << left;
        divider();
        cout << "Fitting trajectory (\"" << _trajectory->id() << "\") to " << _N << " data points:" << endl;
        line();
        cout << setw(30) << "DATA SET"         << setw(20) << "TYPE"     << setw(10) << "POINTS" << endl;
        cout << setw(30) << "----------------" << setw(20) << "--------------" << setw(10) << "-------" << endl;

        for (auto data : _spectrum_data)
        {
            std::string num_points = (_ignore_widths) ? std::to_string(data._x.size()) : std::to_string(data._x.size()) + " x 2";
            cout << setw(30) << data._id  << setw(20)  << "Spectrum"   << setw(10) << num_points << endl;  
        };
        for (auto data : _timelike_data)
        {   
            cout << setw(30) << data._id  << setw(20) << "Timelike" << setw(10) << data._x.size() << endl;  
        };
    };
};