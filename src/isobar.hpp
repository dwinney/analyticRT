// Define an isobar, which is restricted to a single isospin
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef ISOBAR_HPP
#define ISOBAR_HPP

#include "key.hpp"
#include "constants.hpp"
#include "kinematics.hpp"
#include "legendre_P.hpp"
#include "trajectory.hpp"
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace analyticRT
{
    // Forward declare amplitude for typedef below
    class raw_isobar;

    // We want to be able to add isobars together so never use raw_amplitude directly
    // Instead always work with pointers
    using isobar    = std::shared_ptr<raw_isobar>;

    template<class A> 
    inline isobar new_isobar()             {     return std::make_shared<A>(key());              };
    template<class A, class B> 
    inline isobar new_isobar(B x)          {     return std::make_shared<A>(key(), x);           };
    template<class A, class B, class C> 
    inline isobar new_isobar(B x, C y)     {     return std::make_shared<A>(key(), x, y);        };
    template<class A, class B, class C, class D> 
    inline isobar new_isobar(B x, C y, D z){     return std::make_shared<A>(key(), x, y, z);     };
    template<class A, class B, class C, class D, class E> 
    inline isobar new_isobar(B x, C y, D z, E a){ return std::make_shared<A>(key(), x, y, z, a); };

    class raw_isobar
    {
        public: 

        // Basic constructor
        raw_isobar(key x, unsigned int isospin) 
        : _isospin(isospin)
        { 
            check_isospin(isospin);
        };

        // Save a string to id the amplitude
        raw_isobar(key x, unsigned int isospin, std::string id)
        : _id(id), _isospin(isospin)
        {
            check_isospin(isospin);
        };

        // Constructor to be used if initializing a sum
        raw_isobar(key x, unsigned int isospin, std::vector<isobar> vs, std::string id)
        : _id(id), _isospin(isospin), _isobars(vs)
        {
            check_isospin(isospin);

            int n = 0;
            for (auto v : vs)
            {
                n += v->Npars();
            }
            set_Npars(n);
        };

        // ---------------------------------------------------------------------------
        // Getters

        inline std::string id(){ return _id;    };         // string id 
        inline int Npars()    { return _Npars; };         // # of free parameters
        inline unsigned int isospin(){ return _isospin; }; // fixed isospin

        // If the current isobar pointer describes a sum
        inline bool is_sum(){ return !(_isobars.size() == 0); };

        // ---------------------------------------------------------------------------
        // Setters

        inline void set_id(std::string id){ _id = id; };

        // Script-side call to set the parameters
        void set_parameters(std::vector<double> pars);

        // Set options for different things
        // By default do nothing
        virtual void set_option(int x){};

        inline void use_adaptive(bool x){ _adaptive = x; };

        // ---------------------------------------------------------------------------
        // Virtual methods 

        // Amplitude side call to set the parameters which should be overriden by each implementation
        // By default we do nothing
        virtual void allocate_parameters(std::vector<double> pars);

        // Evaluate isobar function as a function of s and zs
        virtual complex evaluate(double s, double zs);

        // Calculate the projection onto the direct or cross channels numerically
        // These can be overridden if the projections can easily be done analytically
        virtual complex direct_projection(unsigned int j, double s);
        virtual complex cross_projection( unsigned int j, double s);

        // Access the saves trajectory pointer
        virtual trajectory get_trajectory(){ return nullptr; };

        // For fitting print a list of labels to differentiate each free parameter
        // By default just print {p[0], p[1], p[2], ...}
        virtual std::vector<std::string> parameter_labels()
        {
            std::vector<std::string> v;
            for (int i = 0; i < _Npars; i++)
            {
                v.push_back( "p["+std::to_string(i)+"]");
            };
            return v;
        };

        protected: 

        // Change the number of free parameters from inside
        inline void set_Npars(unsigned int n){ _Npars = n; };

        friend std::vector<isobar> get_isobars(isobar a);

        // Saved option
        int _option = 0;
        inline void option_error(){ warning(id() + "::set_option()", "Invalid option recieved! Ignoring..."); };
    
        private:

        // A single isobar can be the sum of multiple terms 
        std::vector<isobar> _isobars;

        int _isospin = -1;
        inline void check_isospin(unsigned int iso)
        {
            if (iso > 2) fatal("raw_isobar", "Unphysical isospin passed (" + std::to_string(iso) + ")");
        };

        // Basic properties of the amplitude
        std::string _id = "raw_amplitude";  // String id
        int _Npars = 0;                     // # of free parameters

        // Wheter to use adaptive integration
        bool _adaptive = false;
    };

    // ------------------------------------------------------------------------------
    // External methods for sums of isobars

    inline std::vector<isobar> get_isobars(isobar a){ return a->is_sum() ? a->_isobars : std::vector<isobar>{{a}}; };

    isobar operator+(isobar a, isobar b);

};

#endif