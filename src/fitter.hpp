// Class which allows an isobar and a trajectory to be fit with Minuit
// Currently this is a single channel and single trajectory fit but can be generalized
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef FITTER_HPP
#define FITTER_HPP

#include "kinematics.hpp"
#include "isobar.hpp"
#include "trajectory.hpp"
#include "iterable.hpp"
#include "data_set.hpp"
#include "print.hpp"

#include <stdio.h>
#include <chrono>
#include <string>
#include <vector>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom.h"

namespace analyticRT
{
    // ---------------------------------------------------------------------------
    // Structs for storing relevant info inside the fitter

    // Each free parameter of a model has associated with is a bunch of options
    class parameter
    {
        public:

        // Default initialized
        // This initializes a fixed parameter for normalizations
        parameter()
        : _fixed(true), _value(1.), 
          _label(default_label(0)), _i(-1),
          _lower(0), _upper(5)
        {};
        
        // Constructor for amplitude variables
        parameter(int i)
        : _i(i), _label(default_label(i))
        {};

        int         _i;
        std::string _label;
        std::string _message;
        bool   _fixed         = false;
        double _value         = 0;

        bool   _custom_limits = false;
        double _upper         = 0;
        double _lower         = 0;
        double _step          = 0.1;
        bool   _positive      = false;

        // If this parameter is synced to be equal to another
        bool   _synced        = false;
        int    _sync_to       = -1; // Which param its synced to

        // For iterative fits, we should save the vector of iterated values
        std::vector<double> _itval, _iterr;

        static inline std::string default_label(int i)
        {
            return "par[" + std::to_string(i) + "]";
        };
    };

    // ---------------------------------------------------------------------------
    // Actual fitter object
    // This is templated because it requires an implementation of the chi2 function to be fit
    // The template F should contain the details of the fit and the following static functions
    // fcn(amplitude, std::vector<data_set>&) [function to be minimized, e.g. chi2]
    // 
    template<class F>
    class fitter
    {
        public: 

        // Basic constructor, only requires amplitude to be fit 
        // uses default settings for minuit
        fitter(isobar iso_to_fit, trajectory alpha_to_fit)
        : _isobar(iso_to_fit), _trajectory(alpha_to_fit),
          _minuit(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined"))
        {
            reset_parameters();
        };

        // Parameterized constructor 
        // with explicit choice of minimization strategy and tolerance of minuit routines
        fitter(isobar iso_to_fit, trajectory alpha_to_fit, std::string strategy, double tolerance = 1.E-6)
        : _isobar(iso_to_fit), _trajectory(alpha_to_fit), _tolerance(tolerance),
          _minuit(ROOT::Math::Factory::CreateMinimizer("Minuit2", strategy))
        {
            reset_parameters();
        };

        // -----------------------------------------------------------------------
        // Methods to add data to be fit against

        inline void add_data(data_set data){ _Ndata += data._N; _data.push_back(data); };
        inline void add_data(std::vector<data_set> data){ for (auto datum : data) add_data(datum); };
        inline void clear_data(){ _data.clear(); _Ndata = 0; };

        // -----------------------------------------------------------------------
        // Set limits, labels, and fix parameters

        // Reset labels, limits and options on all parameters
        inline void reset_parameters()
        {
            _pars.clear(); 
            

            // Count number of parameters from amplitude and from trajectory
            std::vector<std::string> iso_labels = {};
            if (_isobar != nullptr)
            {
                _Niso        = _isobar->Npars();
                iso_labels   = _isobar->parameter_labels();
            };
            std::vector<std::string> traj_labels = {};
            if (_trajectory != nullptr)
            {
                _Ntraj      = _trajectory->Npars();
                traj_labels = _trajectory->parameter_labels();
            };

            // Accumulate the labels 
            _Nfree = _Niso + _Ntraj;
            for (int i = 0; i < _Nfree; i++) _pars.push_back(i);

            std::vector<std::string> labels = iso_labels;
            labels.insert( labels.end(), std::make_move_iterator(traj_labels.begin()), std::make_move_iterator(traj_labels.end()));
            set_parameter_labels(labels);
        };

        // Give each parameter a label beyond their default par[i] name
        inline void set_parameter_labels(std::vector<std::string> labels)
        {
            if (labels.size() != _pars.size())
            {
                warning("fitter::set_parameter_labels", "Labels vector does not match number of parameters!");
                return;
            }
            for (int i = 0; i < _pars.size(); i++) _pars[i]._label = labels[i];
        };

        // Set limits and/or a custom stepsize
        inline void set_parameter_limits(parameter& par, std::array<double,2> bounds, double step = 0.1)
        {
            par._custom_limits = true;
            par._lower         = bounds[0];
            par._upper         = bounds[1];
            par._step          = step;
        };

        inline void set_parameter_limits(std::string label, std::array<double,2> bounds, double step = 0.1)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            return set_parameter_limits(_pars[index], bounds, step);
        };

        inline void set_parameter_posdef(parameter& par, bool x = true){ par._positive = x; };
        inline void set_parameter_posdef(std::string label, bool x = true)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            set_parameter_posdef(_pars[index], x);
        };
            
        inline void fix_parameter(parameter& par, double val)
        {
            // If parameter is already fixed, just update the fixed val
            // otherwise flip the fixed flag and update the number of free pars
            if (!par._fixed) _Nfree--;
            par._fixed = true;
            par._value = val;
            _fit = false;
        };

        inline void fix_parameter(std::string label, double val)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            return fix_parameter(_pars[index], val);
        };

        inline void free_parameter(parameter& par)
        {
            // if not fixed, this does nothing
            if (!par._fixed) return;
            par._fixed = false;
            _fit = false;
            _Nfree++;
        };

        inline void free_parameter(std::string label)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            free_parameter(_pars[index]);
        };

        inline void sync_parameter(std::string par, std::string synced_to)
        {
            int i            = find_parameter(par);
            int i_sync_to    = find_parameter(synced_to);
            _pars[i]._synced  = true;
            _pars[i]._sync_to = i_sync_to;

            // Also save present value in case there is one
            _pars[i]._value = _pars[i_sync_to]._value;
            if (!_pars[i_sync_to]._fixed)_Nfree--;
        };

        inline void unsync_parameter(std::string par){ _pars[ find_parameter(par) ]._synced = false; };

        // -----------------------------------------------------------------------
        // Methods related to fit options

        // Set the maximum number of calls minuit will do
        inline void set_max_calls(int n){ print(n); _max_calls = n; };
        
        // Message level for minuit (0-4)
        inline void set_print_level(int n){ _print_level = n; };

        // Change tolerance
        inline void set_tolerance(double tol){ _tolerance = tol; };

        // Change the guess range for initializing parameters
        inline void set_guess_range(std::array<double,2> bounds){ _guess_range = bounds; };

        // Actually do the fit given a vector of size Npars as starting values
        // Prints results to command line but also returns the best-fit chi2 value
        inline void do_fit(std::vector<double> starting_guess, bool show_data = true)
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

        // Same as above (do single fit) but initialize random parameters
        inline void do_fit()
        {            
            // Initial guess
            std::vector<double> guess;  

            // Initialize the guess for each parameter
            std::vector<int> synced_pars;
            for (auto par : _pars)
            {
                if (par._fixed) continue;
                if (par._synced){ synced_pars.push_back(par._i); continue; } // Flag any synced parameters, we'll come back to them

                if (par._custom_limits) guess.push_back(_guesser->Uniform(par._lower, par._upper)); 
                else                    guess.push_back(_guesser->Uniform(_guess_range[0], _guess_range[1]));
                par._value = guess.back();
            };  

            // After all randomized parameters have been set, rego through the synced ones
            for (int i : synced_pars)
            {
                int j = _pars[i]._sync_to;
                _pars[i]._value = _pars[j]._value;
            };

            do_fit(guess);
        };

        // Do N fits with random parameters and return the best fit found
        inline void do_fit(int N)
        {
            
            divider();
            std::cout << std::left << "Finding best fit out of N = " + std::to_string(N) + " attempts." << std::endl;
            if (N <= 0){ line(); return; };

            // Initial guess
            std::vector<double> guess;

            for (int i = 1; i <= N; i++)
            {
                guess.clear();
                // Whether this is the first iteration
                bool first_fit =  !(i-1);
                // Initialize the guess for each parameter
                
                std::vector<int> synced_pars;
                for (auto par : _pars)
                {
                    if (par._fixed) continue;
                    if (par._synced){ synced_pars.push_back(par._i); continue; } // Flag any synced parameters, we'll come back to them

                    if (par._custom_limits) guess.push_back(_guesser->Uniform(par._lower, par._upper)); 
                    else                    guess.push_back(_guesser->Uniform(_guess_range[0], _guess_range[1]));
                    par._value = guess.back();
                };

                // After all randomized parameters have been set, rego through the synced ones
                for (int i : synced_pars)
                {
                    int j = _pars[i]._sync_to;
                    _pars[i]._value = _pars[j]._value;
                };

                // Do out fit with this random guess
                if (!first_fit) std::cout << std::left << "Fit (" + std::to_string(i) + "/" + std::to_string(N) + ")" << std::endl;
                do_fit(guess, first_fit);


                // Compare with previous best and update
                if ( first_fit || fcn() < _best_fcn)
                {
                    _best_fcn_dof = _minuit->MinValue() / (_Ndata - _minuit->NFree());
                    _best_fcn     = _minuit->MinValue();
                    _best_pars    = convert(_minuit->X());
                    _best_errs    = convert(_minuit->Errors());
                };
            };
            
            // After looping, set the best_pars to the amplitude
            allocate_parameters(_best_pars);

            // And set the global saved pars to the best_fit
            std::cout << std::left << "Best fit found after N = " + std::to_string(N) + " iterations" << std::endl;
            line();
            print_results(false);
        };

        inline void do_iterative_fit(int N, std::string filename = "")
        {
            using std::cout; using std::left; using std::endl; using std::setw;
            cout << std::setprecision(8);
            cout << left;

            std::shared_ptr<iterable> iter_trajectory = std::dynamic_pointer_cast<iterable>(_trajectory);
            if (iter_trajectory == nullptr)
            {
                warning("fitter::do_iterative_fit()", "Trajectory (" + _trajectory->id() + ") is not an iterable model! Returning...");
                return;
            };  
            if (N <= 0) return;

            // Do out fit with random guess
            divider(); print("Doing initial fit with uniterated trajectory...");
            do_fit();
            

            // For each subsequent fit we grab the previous best fit parameters
            // and use them as the seed for the new fit after calling trajectory::iterate()
            divider(); print("Commencing iterative refitting..."); line();

            std::vector<double> fcns, last_pars, last_errors, times; 
            for (int i = 0; i <= N; i++)
            {
                auto start = std::chrono::high_resolution_clock::now();

                last_pars   = convert(_minuit->X());
                last_errors = convert(_minuit->Errors());
                std::vector<double> next_pars;

                // We have to filter fixed parameters
                for (auto& par : _pars)
                {
                    if (!par._fixed && !par._synced) next_pars.push_back( last_pars[par._i] );

                    // save each iteration inside the parameter for calling later in summary
                    par._itval.push_back( last_pars[   par._i ] );
                    par._iterr.push_back( last_errors[ par._i ] );
                };
                fcns.push_back(_minuit->MinValue());

                if (i != N)
                {
                    iter_trajectory->iterate();

                    // Refit
                    std::cout << std::left << "Iteration (" + std::to_string(i+1) + "/" + std::to_string(N) + "):" << std::endl;
                    do_fit(next_pars, false);
                }

                // Timing info
                auto stop     = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast< std::chrono::seconds>(stop - start);
                times.push_back(duration.count());
            };

            print_iteration_results(fcns, times);

            if (filename == "") return;

            // Save the summary and the parameters from each iteration to files
            std::string summary = filename + "_summary.txt";
            freopen(summary.c_str(),"w", stdout);
            print_iteration_results(fcns, times);
            fclose(stdout);

            // Also print a file tabulating the parameters of the trajectory
            // This allows us to calculate the "path dependent" iteration
            std::string results = filename + "_pars.txt";
            freopen(results.c_str(),"w", stdout);
            
            for (int i = _Niso; i < _pars.size(); i++) 
            {
                if (_pars[i]._i == _Niso) cout << setw(20) << "#" + _pars[i]._label;
                else                  cout << setw(20) << _pars[i]._label;
            };
            cout << endl;

            for (int n = 0; n < fcns.size(); n++)
            {
                for (int i = _Niso; i < _pars.size(); i++) 
                {
                    if (_pars[i]._fixed) { cout << setw(20) << _pars[i]._value; continue; }
                    if (_pars[i]._synced)
                    { 
                        if (!_pars[_pars[i]._sync_to]._fixed) cout << setw(20) << _pars[_pars[i]._sync_to]._itval[n];
                        else                                  cout << setw(20) << _pars[_pars[i]._sync_to]._value;
                        continue;
                    };
                    cout << setw(20) << _pars[i]._itval[n];
                };
                cout << endl;
            };
            fclose(stdout);
        };

        // Return a vector of best-fit parameters from last fit
        inline std::vector<double> pars(){ return (_fit) ? _fit_pars : std::vector<double>(); };

        // Return value of fit function from last fit
        inline double fcn()     { return (_fit) ? _fcn     : NaN<double>(); };
        inline double fcn_dof() { return (_fit) ? _fcn_dof : NaN<double>(); };

        // Return the pointer to the amp being fit at its current state
        inline isobar     get_isobar()    { return _isobar; };
        inline trajectory get_trajectory(){ return _trajectory; };

        private:

        // These ptrs should point to the amplitude and trajectory to be fit
        isobar     _isobar     = nullptr;
        trajectory _trajectory = nullptr;

        // -----------------------------------------------------------------------
        // Data handling

        int _Ndata = 0;  // Total number of data points
        std::vector<data_set> _data;  // Contain all data

        // -----------------------------------------------------------------------
        // MINUIT handling 

        int _print_level   = 0;     // Error code for MINUIT
        int _max_calls     = 1E6;   // Max calls allowed for minimization fcn
        double _tolerance  = 1.E-6; // Minimization tolerance

        ROOT::Math::Minimizer * _minuit;
        ROOT::Math::Functor _wfcn;

        // Random number generator for creating initial guesses;
        TRandom *_guesser = new TRandom(0);

        // Initialize minuit with all our parameter options etc
        inline void set_up(std::vector<double> starting_guess)
        {
            _minuit->Clear();
            _minuit->SetTolerance(_tolerance);
            _minuit->SetPrintLevel(_print_level);
            _minuit->SetMaxFunctionCalls(_max_calls);

            // Iterate over each _par but also keep track of the index in starting_guess 
            // because parameters might be fixed, these indexes dont necessarily line up
            int i = 0;
            for (auto par : _pars)
            {   
                if (par._fixed || par._synced) continue;
                _minuit->SetVariable(i, par._label, starting_guess[i], par._step);
                if (par._positive)      _minuit->SetVariableLowerLimit(i, 0.);
                if (par._custom_limits) _minuit->SetVariableLimits(i, par._lower, par._upper);
                i++; // move index up
            };
        
            _wfcn = ROOT::Math::Functor(this, &fitter::fit_fcn, _Nfree);
            _minuit->SetFunction(_wfcn);
        };

        // -----------------------------------------------------------------------
        // Calcualtions of chi-squared 
        
        inline void allocate_parameters( std::vector<double> pars)
        {
            std::vector<double> ipars(pars.begin(), pars.begin() + _Niso);
            std::vector<double> tpars(pars.begin() + _Niso, pars.begin() + _Niso + _Ntraj);

            if (_isobar     != nullptr) _isobar->set_parameters(ipars);
            if (_trajectory != nullptr) _trajectory->set_parameters(tpars);
        };

        // // This is the actual function that gets called by minuit
        inline double fit_fcn(const double *cpars)
        { 
            // First convert the C string to a C++ vector
            std::vector<double> pars = convert(cpars);

            // Pass parameters to the amplitude and trajectory
            allocate_parameters(pars);

            // Pass both this and data to fit function
            return F::fcn(_data, _isobar, _trajectory); 
        };

        // Save of the last fit run
        bool _fit = false;      // Whether a fit has already been done or not yet
        double _fcn, _fcn_dof;  // Last saved value of fcn function
        std::vector<double> _fit_pars, _errors;

        // Save of the best fits found if running multiple times
        double _best_fcn, _best_fcn_dof;
        std::vector<double> _best_pars, _best_errs;

        // -----------------------------------------------------------------------
        // Parameter handling

        // Store of parameter info
        std::vector<parameter> _pars;

        // Number of parameters belonging to the isobar and trajectory respectively
        int _Niso  = 0, _Ntraj = 0;

        // Number of free parameters
        int _Nfree = 0;

        // Default guess_range to initalize parameters
        std::array<double,2> _guess_range = {-5, 5};

        // Given a parameter label, find the corresponding index
        inline int find_parameter(std::string label)
        {
            for (auto par : _pars) if (par._label == label) return par._i;
            return error("fitter::find_parameter", "Cannot find parameter labeled " + label + "!", -1);
        };

        // Given a C-style array of size _Nfree
        // Convert to a C++ style std::vector and populate
        // fixed value parameters in the expected order
        inline std::vector<double> convert(const double * cpars)
        {
            std::vector<double> result;

            // Move along the pars index when a parameter is not fixed
            int i = 0;

            std::vector<int> synced_pars;
            for (auto par : _pars)
            {
                if (par._synced) synced_pars.push_back(par._i); 
                if (par._fixed || par._synced) result.push_back(par._value);
                else { result.push_back(cpars[i]); i++; };
            };

            for (int j : synced_pars)
            {
                _pars[j]._value = _pars[_pars[j]._sync_to]._value;
                result[j] = result[_pars[j]._sync_to];
            };

            if (i != _Nfree) warning("fitter::convert", "Something went wrong in converting parameter vector.");
            return result;
        };

        // -----------------------------------------------------------------------
        // Methods to print out status to command line

        // Summary of data sets that have been recieved
        inline void data_info()
        {
            using std::cout; using std::left; using std::endl; using std::setw;
            
            if (_data.size() == 0)
            {
                warning("fitter::data_info", "No data found!"); 
                return;
            };

            cout << left;
            divider();

            std::string iso_name, traj_name, conjunction;

            if (_isobar != nullptr)
            {
                iso_name = "isobar (\"" + _isobar->id() + "\")";
            };
            if (_trajectory != nullptr)
            {
                traj_name = "trajectory (\"" + _trajectory->id() + "\")";
            };
            if (_isobar != nullptr && _trajectory != nullptr) conjunction = " and ";

            cout << "Fitting " << iso_name + conjunction + traj_name << endl;

            line();
            cout << setw(25) << "DATA SET"         << setw(25) << "TYPE     "      << setw(5) << "POINTS" << endl;
            cout << setw(25) << "----------------" << setw(25) << "--------------" << setw(5) << "-------" << endl;
            for (auto data : _data)
            {
                cout << setw(25) << data._id  << setw(25)  << F::data_type(data._type)  << setw(5) << data._N << endl;  
            };
        };

        // Similar summary for parameters
        // Display alongside a vector of current parameter values
        // bool start is whether this is the starting guess vector or the 
        // best fit results
        inline void parameter_info(std::vector<double> starting_guess)
        {
            using std::cout; using std::left; using std::endl; using std::setw;

            cout << std::setprecision(8);
            cout << left;

            line(); divider();
            // Print message at the beginning of the fit
            cout << "Fitting " + std::to_string(_Nfree) << " (of " << std::to_string(_pars.size()) << ") parameters" << endl;
            line();

            // Moving index from the guess vector
            int i = 0;
            
            std::vector<int> synced_pars;
            std::vector<double> vals;
            for (auto &par : _pars)
            {
                if (par._synced)
                {
                    vals.push_back(0);
                    synced_pars.push_back(par._i);
                    par._message = "[= " + std::to_string(par._sync_to) + "]";
                    continue;
                };

                // Or is fixed
                if (par._fixed)
                {
                    vals.push_back(par._value);
                    par._message  =  "[FIXED]";
                    continue;
                }

                par._value = starting_guess[i];
                vals.push_back(par._value);
                if (par._custom_limits)
                {   
                    std::stringstream ss;
                    ss << std::setprecision(5) << "[" << par._lower << ", " << par._upper << "]";
                    par._message = ss.str();
                };
                i++;
            };

            for (int j : synced_pars)
            {
                int k = _pars[j]._sync_to;
                _pars[j]._value = _pars[k]._value;
                vals[j] = _pars[k]._value;
            };

            if (_isobar != nullptr)
            { 
                line(); centered(4, "Isobar parameters"); divider(4); 
                cout << left << setw(10) << "id"     << setw(17) << "PARAMETER"  << setw(20) << "START VALUE"  << endl;
                cout << left << setw(10) << "-----" << setw(17) << "----------" << setw(20) << "------------" << endl;
            };
            for (int k = 0; k < _Niso; k++)
            {
                cout << left << setw(10) << _pars[k]._i << setw(17) << _pars[k]._label << setw(20) << vals[_pars[k]._i] << setw(20) << _pars[k]._message << endl;
            };
            
            if (_trajectory != nullptr)
            { 
                line(); centered(4, "Trajectory parameters"); divider(4); 
                cout << left << setw(10) << "id"     << setw(17) << "PARAMETER"  << setw(20) << "START VALUE"  << endl;
                cout << left << setw(10) << "-----" << setw(17) << "----------" << setw(20) << "------------" << endl;
            }; 
            for (int k = _Niso; k < _Niso + _Ntraj; k++)
            {
                cout << left << setw(10) << _pars[k]._i << setw(17) << _pars[k]._label << setw(20) << vals[_pars[k]._i] << setw(20) << _pars[k]._message << endl;
            };

            line(); divider(); line();
        };

        // After a fit return a summary of fit results
        // At the end of a fit, print out a table sumarizing the fit results
        // if last_fit == true, we grab the results from the most recent fit in _minuit
        // else we print out the ones saved in _best_fit
        inline void print_results(bool last_fit = true)
        {
            using std::cout; using std::left; using std::endl; using std::setw;

            cout << std::setprecision(8);
            cout << left;

            int dof                  = _Ndata - _minuit->NFree();
            double fcn               = (last_fit) ? _minuit->MinValue()               : _best_fcn;
            double fcn_dof           = (last_fit) ? _minuit->MinValue() / double(dof) : _best_fcn_dof;
            std::vector<double> pars = (last_fit) ? convert(_minuit->X())             : _best_pars;
            std::vector<double> errs = (last_fit) ? convert(_minuit->Errors())        : _best_errs;

            divider();
            std::cout << std::left << std::setw(5)  << "fcn = "       << std::setw(15) << fcn     << std::setw(5) << "";
            std::cout << std::left << std::setw(10) << "fcn/dof = "  << std::setw(15) << fcn_dof << "\n";

            for (int i = 0; i < pars.size(); i++) _pars[i]._value = pars[i];

            if (_isobar != nullptr)
            {
                line(); centered(4, "Isobar parameters"); divider(4);            
                cout << left << setw(10) << "id"     << setw(16) << "PARAMETER"  << setw(18) << "FIT VALUE"    << setw(18) << "ERROR"        << endl;
                cout << left << setw(10) << "-----" << setw(16) << "----------" << setw(18) << "------------" << setw(18) << "------------" << endl;
            };
            for (int k = 0; k < _Niso; k++)
            {
                parameter par = _pars[k];
                std::stringstream ss;
                ss << std::setprecision(8) << errs[par._i];
                std::string err = !(par._fixed || par._synced) ? ss.str() : par._message;
                cout << left << setw(10) << par._i << setw(17) << par._label << setw(20) << par._value << setw(20) << err << endl;
            };

            if (_trajectory != nullptr)
            {
                line(); centered(4, "Trajectory parameters"); divider(4);            
                cout << left << setw(10) << "id"     << setw(17) << "PARAMETER"  << setw(20) << "FIT VALUE"    << setw(20) << "ERROR"        << endl;
                cout << left << setw(10) << "-----"  << setw(17) << "----------" << setw(20) << "------------" << setw(20) << "------------" << endl;
            };
            for (int k = _Niso; k < _Niso + _Ntraj; k++)
            {
                parameter par = _pars[k];
                std::stringstream ss;
                ss << std::setprecision(8) << errs[par._i];
                std::string err = !(par._fixed || par._synced) ? ss.str() : par._message;
                cout << left << setw(10) << par._i << setw(17) << par._label << setw(20) << par._value << setw(20) << err << endl;
            };


            line(); divider(); line();
            
            // At the end update the amplitude parameters to include the fit results
            allocate_parameters(pars);

            _fcn      = fcn;
            _fcn_dof  = fcn_dof;
            _fit_pars = pars;

            // Let rest of the fitter that a fit result has been saved
            _fit = true;
        };

        inline void print_iteration_results(std::vector<double>& fcns, std::vector<double>& times)
        {
            using std::cout; using std::left; using std::endl; using std::setw;
            cout << std::setprecision(8);
            cout << left;


            // Print a summary of all the results
            divider(); line();
            centered(4, "SUMMARY OF ITERATIVE FIT RESULTS");
            line(); divider(); line();

            double avg = 0;
            for (auto t : times){ avg += t; };
            avg /= times.size();
            cout << std::setprecision(2);
            cout << setw(30) << "Number of iterations: " << fcns.size()-1 << endl;
            cout << setw(30) << "Average time per iteration: " << avg << " s" << endl; line();
            cout << std::setprecision(8);

            data_info();
            line();

            divider(); print("Parameters excluded from fits:");
            print("----------------------------------------");
            for (auto par : _pars)
            {
                if (par._fixed)  cout << setw(15) << "-- " + par._label  << setw(15) << " fixed to "    << setw(15) << std::to_string(par._value) << endl;
                if (par._synced) cout << setw(15) << "-- " + par._label  << setw(15) << " synced with " << setw(15) << _pars[par._sync_to]._label << endl;
            };
            line();

            divider();
            double dof = _Ndata - _minuit->NFree();
            cout << setw(15) << "iter" << setw(25) << "fcn" << setw(25) << "fcn/dof" << endl; 
            cout << setw(15) << "-----" << setw(25) << "-------------" << setw(25) << "-------------" << endl; 
            for (int n = 0; n < fcns.size(); n++) cout << setw(15) << n << setw(25) << fcns[n] << setw(25) << fcns[n]/dof << endl; 
            line(); divider(); line();

            for (auto par : _pars)
            {
                if (par._fixed || par._synced) continue;

                centered(4, par._label); divider();
                cout << setw(15) << "iter" << setw(25) << "FIT VALUE" << setw(25) << "ERROR" << endl; 
                cout << setw(15) << "-----" << setw(25) << "-------------" << setw(25) << "-------------" << endl; 
                for (int j = 0; j < fcns.size(); j++) cout << setw(15) << j << setw(25) << par._itval[j] << setw(25) << par._iterr[j] << endl; 
                line(); divider(); line();
            };
        };
    };
};

#endif