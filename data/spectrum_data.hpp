// Methods to get data related to the trajectory itself
// This includes particle spectra for spacelike energies
//  and effective trajectory data for timelike energies
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef SPECTRUM_DATA_HPP
#define SPECTRUM_DATA_HPP

#include "data_set.hpp"
#include "elementwise.hpp"

namespace analyticRT
{
    // --------------------------------------------------------------------------]
    // Filter spectra to isolate only odd or even spin resonances

    template<int N>
    inline std::array<std::vector<double>,N> filter(bool odd, std::array<std::vector<double>,N> data)
    {
        std::array<std::vector<double>,N> out;
        for (int i = 0; i < data[0].size(); i++)
        {
            int J = int(data[0][i]);
            if (J%2 != odd) continue;

            for (int k = 0; k < N; k++)
            {
                out[k].push_back(data[k][i]);
            };
        };

        return out;
    };

    // --------------------------------------------------------------------------]
    // Particle spectra

    std::string isoscalar_file = "/data/isoscalar_spectrum.dat";
    std::string isovector_file = "/data/isovector_spectrum.dat";

    inline data_set EXD_spectrum(std::string input, std::string id)
    {
        auto raw      = import_data<6>(input);
        int N         = check<6>(raw, id);

        data_set output;
        output._id   = id;
        output._N    = N;
        output._type = 0;

        output._x   = raw[2];
        output._dx  = raw[3];
        output._y   = raw[0];
        output._dy  = raw[1];
        output._z   = raw[4];
        output._dz  = raw[5];

        return output;
    };
    inline data_set isovector_spectrum(){ return EXD_spectrum(isovector_file, "Isovector Spectrum");};
    inline data_set isoscalar_spectrum(){ return EXD_spectrum(isoscalar_file, "Isoscalar Spectrum");};

    inline data_set filtered_spectrum(bool odd_parity, std::string input, std::string id)
    {
        auto raw      = import_data<6>(input);
        auto filtered = filter<6>(odd_parity, raw);
        int N         = check<6>(filtered, id);

        data_set output;
        output._id   = id;
        output._N    = N;
        output._type = 0;

        output._x   = filtered[2];
        output._dx  = filtered[3];
        output._y   = filtered[0];
        output._dy  = filtered[1];
        output._z   = filtered[4];
        output._dz  = filtered[5];

        return output;
    };
    inline data_set rho_spectrum()  { return filtered_spectrum(true,  isovector_file, "rho Mesons");};
    inline data_set a_spectrum()    { return filtered_spectrum(false, isovector_file, "a Mesons");  };
    inline data_set omega_spectrum(){ return filtered_spectrum(true,  isoscalar_file, "omega Mesons");};
    inline data_set f_spectrum()    { return filtered_spectrum(false, isoscalar_file, "f Mesons");  };
};

#endif