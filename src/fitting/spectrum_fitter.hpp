// Class for taking in mass/width data and performing a fit a given trajectory model
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef SPECTRUM_FIT_HPP
#define SPECTRUM_FIT_HPP

#include "fitter.hpp"
#include "trajectory.hpp"

namespace analyticRT
{
    class spectrum_fitter : public fitter
    {
        public: 

        // Basic constructor, only requires amplitude to be fit 
        // uses default settings for minuit
        spectrum_fitter(trajectory to_fit)
        : _trajectory(to_fit)
        {
            // Populate the _pars vector with the appropriate sized array
            reset_parameters();
        };

        // Parameterized constructor 
        // with explicit choice of minimization strategy and tolerance of minuit routines
        spectrum_fitter(trajectory to_fit, std::string strategy, double tolerance = 1.E-6)
        : _trajectory(to_fit), fitter(strategy, tolerance)
        {
            reset_parameters();
        };

        // -----------------------------------------------------------------------
        // Set limits, labels, and fix parameters

        // Reset labels, limits and options on all parameters
        void reset_parameters();

        // -----------------------------------------------------------------------
        // Methods to add data to be fit against

        // Parse a data set and add it to the fitting pool
        // Currently supported data_types are integrated and differential x-sections
        void add_data(data_set data);

        // Remove all saved data 
        void clear_data();

        protected:

        // This ptr should point to the trajectory to be fit
        trajectory _trajectory = nullptr;

        inline void pass_pars(std::vector<double> pars){ _trajectory->set_parameters(pars); };
        inline void iterate(){ _trajectory->iterate(); };

        // -----------------------------------------------------------------------
        // Data handling

        // Masses and widths
        std::vector<data_set> _timelike_data, _spectrum_data; 
        
        // Summary of data sets that have been recieved
        void data_info();

        // -----------------------------------------------------------------------
        // Calcualtions of chi-squared 
        
        // Calculate the chi2 for a given set of parameters, pars, 
        // from a given timelike or spectrum data set
        double chi2_timelike(data_set data);
        double chi2_spectrum(data_set data);

        // This is the actual function that gets called by minuit
        // Combined chi2 from all observables and data sets
        double fit_chi2(const double *pars);


    };
};

#endif