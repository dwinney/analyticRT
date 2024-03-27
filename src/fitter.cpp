// Class for taking in data and performing a fit a given trajectory model
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "fitter.hpp"

namespace analyticRT
{
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
        
        // After looping, set the best_pars to the amplitude being fit
        pass_pars(_best_pars);

        // And set the global saved pars to the best_fit
        print( "Best fit found after N = " + std::to_string(N) + " iterations" );
        line();
        print_results(false);
    };

    void fitter::iterative_fit(int N)
    {
        // Do the first fit on the "zeroth" solution: 
         if (N == 0) return;

        std::vector<std::array<double,2>> chi2s;

        divider(); print("Doing initial fit to uniterated trajectory...");

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

        divider(); print("Commencing iterative refitting..."); line();

        // For each subsequent fit we grab the previous best fit parameters
        // and use them as the seed for the new fit after calling trajectory::iterate()
        std::vector<double> last_pars, last_errors; 
        for (int i = 0; i < N; i++)
        {
            last_pars   = convert(_minuit->X());
            last_errors = convert(_minuit->Errors());
            std::vector<double> next_pars;

            // We have to filter fixed parameters
            for (auto& par : _pars)
            {
                if (!par._fixed) next_pars.push_back( last_pars[par._i] );

                // save each iteration inside the paramter for calling later in summar
                par._itval.push_back( last_pars[   par._i ] );
                par._iterr.push_back( last_errors[ par._i ] );
            };
            chi2s.push_back({_minuit->MinValue(), _minuit->MinValue() / (_N - _minuit->NFree())});

            iterate();
            std::cout << std::left << "Iteration (" + std::to_string(i+1) + "/" + std::to_string(N) + "):" << std::endl;

            do_fit(next_pars, false);
        };

        last_pars   = convert(_minuit->X());
        last_errors = convert(_minuit->Errors());
        for (auto& par : _pars)
        {
            par._itval.push_back( last_pars[   par._i ] );
            par._iterr.push_back( last_errors[ par._i ] );
        };
        chi2s.push_back({_minuit->MinValue(), _minuit->MinValue() / (_N - _minuit->NFree())});

        line(); divider(); line();
        divider(3); print("ITERATION", "chi2", "chi2/dof"); divider(3);
        for (int n = 0; n < N+1; n++)
        {
            print(n, chi2s[n][0], chi2s[n][1]);
        };
        
        line(); divider(3); print("ITERATION", "FIT VALUE", "ERROR"); divider(3);
        for (auto par : _pars)
        {
            if (par._fixed) continue;

            divider(3); centered(3, par._label); divider(3);

            for (int j = 0; j < N + 1; j++)
            {
                print(j, par._itval[j], par._iterr[j]);
            };
            divider(3); line();
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
        
        pass_pars(pars);

        _chi2     = chi2;
        _chi2dof  = chi2dof;
        _fit_pars = pars;

        // Let rest of the fitter that a fit result has been saved
        _fit = true;
    }
};