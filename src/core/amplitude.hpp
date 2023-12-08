// Assemble a bunch of isobars together into a total amplitude with the correct
// crossing structure
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef AMPLITUDE_HPP
#define AMPLITUDE_HPP

#include "key.hpp"
#include "constants.hpp"
#include "kinematics.hpp"
#include "isobar.hpp"

namespace analyticRT
{
    class raw_amplitude;

    using amplitude = std::shared_ptr<raw_amplitude>;

    inline amplitude new_amplitude(std::string id = "amplitude"){ return std::make_shared<raw_amplitude>(key(), id); };

    template<class A>
    inline amplitude new_amplitude(A a, std::string id){ return std::make_shared<raw_amplitude>(key(), a, id); };

    class raw_amplitude
    {
        public: 
        
        // Constructors
        raw_amplitude(key x, std::string id)
        : _id(id)
        {};

        raw_amplitude(key x, isobar iso, std::string id)
        : _isobars({iso}), _id(id)
        {};

        raw_amplitude(key x, std::vector<isobar> isobars, std::string id)
        : _isobars(isobars), _id(id)
        {};

        // Acces sdirect and cross channel pieces seperately
        // Or use evaluate to get the sum
        
        complex direct_channel( unsigned int i, double s, double zs);
        complex cross_channels( unsigned int i, double s, double zs);

        // Evaluate full isopin amplitude by iterating all the isobars
        complex evaluate(unsigned int i, double s, double zs)
        {
            return direct_channel(i, s, zs) + cross_channels(i, s, zs);
        };

        // Access the direct and cross projection of the full amplitude
        complex direct_projection(unsigned int i, unsigned int j, double s);
        complex cross_projection( unsigned int i, unsigned int j, double s);

        // Evaluate the full partial wave projection
        inline complex partial_wave(unsigned int i, unsigned int j, double s)
        {
            return direct_projection(i, j, s) + cross_projection(i, j, s);
        };

        // Manipulate which isobars contribute to the amplitude
        inline void add_isobar(isobar x) { _isobars.push_back(x); };
        inline void add_isobars(std::vector<isobar> v){ for (auto i : v) _isobars.push_back(i); };
        inline void clear_isobars(){ _isobars.clear(); };

        // Whether or not to set all cross-channel contributions to zero
        inline void ignore_cross(bool x) { _ignore_cross = x; };

        // Access the string id
        inline std::string id(){ return _id; };
        
        private:

        // String identifier
        std::string _id = "amplitude";

        // Constituent isobars
        std::vector<isobar> _isobars;

        // Whether or not we ignore cross channe contributions
        bool _ignore_cross = false;
    };
};

#endif