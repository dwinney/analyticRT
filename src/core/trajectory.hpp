// Abstract class for a general Regge trajectory defined by an once-subtracted
// dispersion relation across RHC and optionally also a LHC
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include "constants.hpp"

#include <memory>
#include <cmath>
#include <complex>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/differentiation/finite_difference.hpp>

namespace analyticRT
{
    // Easiest interface between different models of trajectories
    // is if everything is a pointer of the abstract class
    class raw_trajectory;
    using trajectory = std::shared_ptr<raw_trajectory>;

    // Instead of the constructors of raw_trajectory we use this function
    // to create trajectory pointers which can then be put into amplitudes
    template<class A>
    trajectory new_trajectory(std::string id = "trajectory"){ return std::make_shared<A>(4*M2_PION, id); };
    template<class A, class B>
    trajectory new_trajectory(double R_th, B b){ return std::make_shared<A>(R_th, b); };
    template<class A, class B, class C>
    trajectory new_trajectory(double R_th, B b, C c){ return std::make_shared<A>(R_th, b, c); };
    template<class A, class B, class C, class D>
    trajectory new_trajectory(double R_th, B b, C c, D d){ return std::make_shared<A>(R_th, b, c, d); };
    template<class A, class B, class C, class D, class E>
    trajectory new_trajectory(double R_th, B b, C c, D d, E e){ return std::make_shared<A>(R_th, b, c, d, e); };

    // Each model must be derived from this virtual class
    class raw_trajectory
    {
        // -------------------------------------------------------------------
        public:

        // Constructor with only a RHC
        raw_trajectory(double R_threshold, std::string id)
        : _sRHC(R_threshold),  _hasLHC(false), _id(id)
        {};

        // Constructor with only a RHC
        raw_trajectory(double R_threshold, int N, std::string id)
        : _sRHC(R_threshold),  _hasLHC(false), _id(id), _Npars(N)
        {};

        // Constructor specifying both thresholds
        raw_trajectory(double R_threshold, double L_threshold, std::string id)
        : _sRHC(R_threshold), _sLHC(L_threshold), _hasLHC(true), _id(id)
        {};

        // Constructor specifying both thresholds
        raw_trajectory(double R_threshold, double L_threshold, int N, std::string id)
        : _sRHC(R_threshold), _sLHC(L_threshold), _hasLHC(true), _id(id), _Npars(N)
        {};

        // Public getters and setters
        inline std::string id(){ return _id; };
        inline void set_id(std::string id){ _id = id; };

        // Use set_parameters to set free variables
        // this checks that the size is the expected size
        void set_parameters( std::vector<double> pars);
        inline int Nfree(){ return _Npars; };

        // This function actually distributes paramters to the model and must be specified
        virtual inline void allocate_parameters(std::vector<double> pars){ return; };

        // Assign each free parameter a label to identify it in fitters and such
        virtual std::vector<std::string> parameter_labels();

        // Set an optional int flag to change model
        // By default this does nothing but save the flag value
        virtual inline void set_option( int opt ){ _option = opt; };
        inline int option(){ return _option; };

        // Similar flag for debugging purposes
        virtual inline void set_debug( int opt ){ _debug = opt; };
        inline int debug(){ return _debug; };

        // Set the subtraction location for the dispersion relation
        inline void set_subtraction(double s_sub, double val){ _sSUB = s_sub; _alphaSUB = val; };

        // Many models can form iterative solutions between real and imaginary part
        // here you can set how many iterations to do
        inline  void max_iterations(int x){ _Niters = x; };
        virtual void iterate(){};

        // Output the trajectory from evaluating dispersion relation
        virtual complex evaluate(double s);

        inline void set_initial_guess(std::array<double,2> alpha, double mu2 = 20) 
        {
            _azInitial = alpha[0]; _apInitial = alpha[1]; _mu2 = mu2;
            this->initialize();
        };
        inline double initial_guess(double s){ return (_azInitial + _apInitial * s) / sqrt(1. + s / _mu2); };

        // Output the real or imaginary parts along real line
        double real_part(double s);
        double imaginary_part(double s);

        // In the resonance regime we can approximate the real part as a linear function
        double width(double s);

        // -------------------------------------------------------------------
        protected:

        virtual void initialize(){};

        // Change number of expected parameters
        // We dont want this changing outside of the class itself
        inline void set_Npars(int n){ _Npars = n; };

        // RHC must always be specified by model
        virtual double RHC(double s){ return 0.; };

        // LHC is optional and will default to 0
        inline virtual double LHC(double s){ return 0.; };

        // Threshold energies for LHC and RHC
        double _sRHC, _sLHC;
        bool   _hasLHC = false; // Whether or not to automatically skip LHC integral

        // Once-subtracted DR parameters
        double _sSUB = 0., _alphaSUB = 0.;

        // Maximum number of iterations to calculate
        int _Niters = 0;

        // In iterative solutions we will need to be able to change the inital guess at the real part        
        double _azInitial = 0.5, _apInitial = 0.9, _mu2 = 20;

        // -------------------------------------------------------------------
        private:

        // Values of dispersion relations from each of the two cuts
        // This shouldnt be needed outside of internal functions 
        complex DR_RHC(double s);

        // String identifier to differentiate models at runtime
        std::string _id = "trajectory";

        // Number of free parameters
        int _Npars  = 0;

        // Some flag to change model options
        int _option = 0, _debug = 0;
    };
};

#endif