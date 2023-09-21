// Class which allows an amplitude to be fit to data based on chi2 minimization
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "fitter.hpp"

namespace analyticRT
{
    // -----------------------------------------------------------------------
    // Methods for managing data

    void fitter::add_data(data_set data)
    {
        switch (data._type)
        {
            case spectrum:
            {
                _spectrum_data.push_back(data); 
                _N += 2*data._N;
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

    void fitter::clear_data()
    {
        _N = 0;
        _spectrum_data.clear();
        _timelike_data.clear();
    };

    // -----------------------------------------------------------------------
    // Set limits, labels, and fix parameters

    void fitter::reset_parameters()
    {
        _pars.clear();
        _Nfree = _trajectory->Npars();

        // populate parameters vector of appropriate size
        for (int i = 0; i < _trajectory->Npars(); i++)
        {
            _pars.push_back(i);
        };
        set_parameter_labels(_trajectory->parameter_labels());
    };

    // Give each parameter a label beyond their default par[i] name
    void fitter::set_parameter_labels(std::vector<std::string> labels)
    {
        if (labels.size() != _pars.size())
        {
            warning("fitter::set_parameter_labels", "Labels vector does not match number of parameters!");
            return;
        }

        for (int i = 0; i < _pars.size(); i++)
        {
            _pars[i]._label = labels[i];
        };
    };

    // Given a parameter label, find the corresponding index
    int fitter::find_parameter(std::string label)
    {
        for (auto par : _pars)
        {
            if (par._label == label) return par._i;
        };
        return error("fitter::find_parameter : Cannot find parameter labeled " + label + "!", -1);
    };

    // Set limits and/or a custom stepsize
    void fitter::set_parameter_limits(parameter& par, std::array<double,2> bounds, double step)
    {
        par._custom_limits = true;
        par._lower         = bounds[0];
        par._upper         = bounds[1];
        par._step          = step;
    };

    void fitter::set_parameter_limits(std::string label, std::array<double,2> bounds, double step)
    {
        int index = find_parameter(label);
        if (index == -1) return;
        return set_parameter_limits(_pars[index], bounds, step);
    }

    void fitter::fix_parameter(parameter& par, double val)
    {
        // If parameter is already fixed, just update the fixed val
        // otherwise flip the fixed flag and update the number of free pars
        if (!par._fixed) _Nfree--;

        par._fixed = true;
        par._value = val;
        _fit = false;
    };

    void fitter::fix_parameter(std::string label, double val)
    {
        int index = find_parameter(label);
        if (index == -1) return;
        return fix_parameter(_pars[index], val);
    };

    void fitter::free_parameter(parameter& par)
    {
        // if not fixed, this does nothing
        if (!par._fixed) return;

        par._fixed = false;
        _fit = false;
    };

    void fitter::free_parameter(std::string label)
    {
        int index = find_parameter(label);
        return free_parameter(_pars[index]);
    };

    // Given a C-style array of size _Nfree
    // Convert to a C++ style std::vector and populate
    // fixed value parameters in the expected order
    std::vector<double> fitter::convert(const double * cpars)
    {
        std::vector<double> result;

        // Move along the pars index when a parameter is not fixed
        int i = 0;

        for (auto par : _pars)
        {
            if (par._fixed)
            {
                result.push_back(par._value);
            }
            else 
            {
                result.push_back(cpars[i]);
                i++;
            }
        };

        if (i != _Nfree) warning("fitter::convert", "Something went wrong in converting parameter vector.");

        return result;
    };

    // ---------------------------------------------------------------------------
    // Chi-squared functions

    // Total chi2, this combines all data sets and is the function that gets called 
    // by the minimizer
    double fitter::fit_chi2(const double * cpars)
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
    double fitter::chi2_spectrum(data_set data)
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

            double width_th  = _trajectory->width(s);
            double width_ex  = data._z[i];
            double width_err = data._dz[i];
            chi2 += pow((width_th - width_ex)/width_err, 2);
        };

        return chi2;
    };

    // Calculate chi2 from differential data set
    double fitter::chi2_timelike(data_set data)
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
    // Load parameter information into minuit

    void fitter::set_up(std::vector<double> starting_guess)
    {
        _minuit->Clear();
        _minuit->SetTolerance(_tolerance);
        _minuit->SetPrintLevel(_print_level);
        _minuit->SetMaxFunctionCalls(_max_calls);

        // Iterate over each _par but also keep track of the index in starting_guess
        int i = 0;
        for ( auto par : _pars)
        {   
            if (par._fixed) continue;
            
            _minuit->SetVariable(i, par._label, starting_guess[i], par._step);

            if (par._custom_limits) _minuit->SetVariableLimits(i, par._lower, par._upper);

            // If we've made it this far, this par isnt fixed, and we move to the next index 
            i++;
        };
    
        _fcn = ROOT::Math::Functor(this, &fitter::fit_chi2, _Nfree);
        _minuit->SetFunction(_fcn);
    };

    // Actually do the fit given a vector of size amp->Npars() as starting values
    // Prints results to command line and sets best fit parameters into _amplitude
    void fitter::do_fit(std::vector<double> starting_guess, bool show_data)
    {
        if (starting_guess.size() != _Nfree) 
        {
            warning("fitter::do_fit", "Starting guess not the correct size! Expected " + std::to_string(_Nfree) + " parameters!");
            return;
        };

        set_up(starting_guess);

        if (show_data) { line(); data_info(); };
        parameter_info(starting_guess);


        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Beginning fit..." << std::flush; 

        if (_print_level != 0) line();   
        _minuit->Minimize();
        if (_print_level != 0) line();   

        std::cout << "Done! \n";

        // Timing info
        auto stop     = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast< std::chrono::seconds>(stop - start);
        std::cout << std::left << "Finished in " << duration.count() << " s" << std::endl;

        line();
        print_results();
    };

    // Repeat do_fit N times and find the best fit
    // Parameters are randomly initialized each time on the interval [-5, 5] unless custom limits are set
    void fitter::do_fit(int N)
    {
        if (N == 0) return;

        divider();
        std::cout << std::left << "Commencing N = " + std::to_string(N) + " fit iterations." << std::endl;

        // Initial guess
        std::vector<double> guess;

        for (int i = 1; i <= N; i++)
        {
            guess.clear();

            // Whether this is the first iteration
            bool first_fit =  !(i-1);
            
            // Initialize the guess for each parameter
            for (auto par : _pars)
            {
                if (par._fixed) continue;

                if (par._custom_limits) guess.push_back(_guesser->Uniform(par._lower, par._upper)); 
                else                    guess.push_back(_guesser->Uniform(_guess_range[0], _guess_range[1]));
            };

            // Do out fit with this random guess
            if (!first_fit) std::cout << std::left << "Fit (" + std::to_string(i) + "/" + std::to_string(N) + ")" << std::endl;
            do_fit(guess, first_fit);

            // Compare with previous best and update
            if ( first_fit || chi2dof() < _best_chi2dof)
            {
                _best_chi2dof = _minuit->MinValue() / (_N - _minuit->NFree());
                _best_chi2    = _minuit->MinValue();
                _best_pars    = convert(_minuit->X());
                _best_errs    = convert(_minuit->Errors());
            };
        };

        // After looping, set the best_pars to the amplitude
        _trajectory->set_parameters(_best_pars);

        // And set the global saved pars to the best_fit
        std::cout << std::left << "Best fit found after N = " + std::to_string(N) + " iterations" << std::endl;
        line();
        print_results(false);
    };

    void fitter::iterative_fit(int N)
    {
        // Do the first fit on the "zeroth" solution: 
         if (N == 0) return;

        std::vector< std::vector<double> > values, errors;
        std::vector<std::array<double,2>> chi2s;

        divider();
        std::cout << std::left << "Doing initial fit to uniterated trajectory..." << std::endl;

        // Initial guess
        std::vector<double> guess;
        
        // Initialize the guess for each parameter
        for (auto par : _pars)
        {
            if (par._fixed) continue;

            if (par._custom_limits) guess.push_back(_guesser->Uniform(par._lower, par._upper)); 
            else                    guess.push_back(_guesser->Uniform(_guess_range[0], _guess_range[1]));
        };

        // Do out fit with this random guess
        do_fit(guess, true);

        // And set the global saved pars to the best_fit
        line();
        print_results(true);

        divider();
        std::cout << std::left << "Commencing iterative refitting..." << std::endl;
        line();

        for (int i = 0; i < N; i++)
        {
            std::vector<double> last_pars = convert(_minuit->X());
            values.push_back(last_pars);
            errors.push_back(convert(_minuit->Errors()));
            chi2s.push_back({_minuit->MinValue(), _minuit->MinValue() / (_N - _minuit->NFree())});

            _trajectory->iterate();
            std::cout << std::left << "Iteration (" + std::to_string(i+1) + "/" + std::to_string(N) + "):" << std::endl;

            do_fit(last_pars, false);
        };
        values.push_back(convert(_minuit->X()));
        errors.push_back(convert(_minuit->Errors()));
        chi2s.push_back({_minuit->MinValue(), _minuit->MinValue() / (_N - _minuit->NFree())});

        line();
        divider();
        line();

        divider(3);
        print("ITERATION", "chi2", "chi2/dof");
        divider(3);
        for (int n = 0; n < N+1; n++)
        {
            print(n, chi2s[n][0], chi2s[n][1]);
        };
        
        line();
        divider(3);
        print("ITERATION", "FIT VALUE", "ERROR");
        divider(3);
        std::vector<std::string> labels = _trajectory->parameter_labels();
        for (int i = 0; i < _trajectory->Npars(); i++)
        {
            divider(3);
            centered(3, labels[i]);
            divider(3);
            for (int j = 0; j < N + 1; j++)
            {
                print(j, values[j][i], errors[j][i]);
            };
            divider(3);
            line();
        };
    };

    // ---------------------------------------------------------------------------
    // Status messages printed to command line

    // Summary of data sets that have been recieved
    void fitter::data_info()
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
            cout << setw(30) << data._id  << setw(20)  << "Spectrum"   << setw(10) << std::to_string(data._x.size()) + " x 2"  << endl;  
        };
        for (auto data : _timelike_data)
        {   
            cout << setw(30) << data._id  << setw(20) << "Timelike" << setw(10) << data._x.size() << endl;  
        };
    };

    // Print out a little table of the starting values of all the parameters and other set options
    // Here starting_guess should be of size _Nfree!
    void fitter::parameter_info(std::vector<double> starting_guess)
    {  
        using std::cout; 
        using std::left;
        using std::setw;
        using std::endl;

        cout << std::setprecision(8);
        cout << left;

        line(); divider();
        // Print message at the beginning of the fit
        cout << "Fitting " + std::to_string(_Nfree) << " (of " << std::to_string(_pars.size()) << ") parameters:" << endl;

        line();

        cout << left << setw(10) << "N"     << setw(17) << "PARAMETER"  << setw(20) << "START VALUE"  << endl;
        cout << left << setw(10) << "-----" << setw(17) << "----------" << setw(20) << "------------" << endl;

        // Moving index from the guess vector
        int i = 0;
        for (auto par : _pars)
        {
            // Parse whether a parameter has extra options 
            // such as custom limits
            std::string extra = "";
            if (par._custom_limits)
            {   
                std::stringstream ss;
                ss << std::setprecision(5) << "[" << par._lower << ", " << par._upper << "]";
                extra = ss.str();
            };

            // Or is fixed
            double par_val;
            if (par._fixed)
            {
                par_val = par._value;
                extra   = "[FIXED]";
            }
            else
            {
                par_val = starting_guess[i];
                i++;
            }

            cout << left << setw(10) << par._i << setw(17) << par._label << setw(20) << par_val << setw(20) << extra << endl;
        };
        line(); divider(); line();
    };

    // At the end of a fit, print out a table sumarizing the fit results
    // if last_fit == true, we grab the results from the most recent fit in _minuit
    // else we print out the ones saved in _best_fit
    void fitter::print_results(bool last_fit)
    {
        using std::cout; 
        using std::left;
        using std::setw;
        using std::endl;

        cout << std::setprecision(8);
        cout << left;

        int dof                  = _N - _minuit->NFree();
        double chi2              = (last_fit) ? _minuit->MinValue()               : _best_chi2;
        double chi2dof           = (last_fit) ? _minuit->MinValue() / double(dof) : _best_chi2dof;
        std::vector<double> pars = (last_fit) ? convert(_minuit->X())             : _best_pars;
        std::vector<double> errs = (last_fit) ? convert(_minuit->Errors())        : _best_errs;


        divider();
        std::cout << std::left << std::setw(10) << "chi2 = "      << std::setw(15) << chi2;
        std::cout << std::left << std::setw(10) << "chi2/dof = "  << std::setw(15) << chi2dof << "\n";

        line();

        cout << left << setw(10) << "N"     << setw(16) << "PARAMETER"  << setw(18) << "FIT VALUE"    << setw(18) << "ERROR"        << endl;
        cout << left << setw(10) << "-----" << setw(16) << "----------" << setw(18) << "------------" << setw(18) << "------------" << endl;

        for (auto par : _pars)
        {
            double val;
            std::string err;
            std::stringstream ss;
            ss << std::setprecision(8);

            if (par._fixed)
            {
                val = par._value;
                err = "[FIXED]";
            }
            else
            {
                val = pars[par._i];
                ss << errs[par._i];
                err = ss.str();
            }

            cout << left << setw(10) << par._i << setw(16) << par._label << setw(18) << val << setw(18) << err << endl;
        };
        line(); divider(); line();
        
        // At the end update the amplitude parameters to include the fit results
        _trajectory->set_parameters(pars);

        _chi2     = chi2;
        _chi2dof  = chi2dof;
        _fit_pars = pars;

        // Let rest of the fitter that a fit result has been saved
        _fit = true;
    }
};