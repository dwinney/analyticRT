// Functions for strings and printing things to the commandline
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef PRINT_HPP
#define PRINT_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <complex> 

namespace analyticRT
{
    // Default values
    const int TEXT_WIDTH       = 62;
    const int PRINT_SPACING    = 15;
    const int PRINT_PRECISION  = 9;    
    const int STRING_PRECISION = 3;
    const std::string UNIT_DIV = std::string(PRINT_SPACING, '-');

    // ---------------------------------------------------------------------------   
    // Output an empty line to the terminal
    inline void line()
    {
        std::cout << std::endl;
    };

    // Print out a horizontal line
    inline void divider()
    {
        std::cout << std::string(TEXT_WIDTH, '-') << std::endl;
    };

    inline void divider(int n)
    {
        std::string div;
        for (int i = 0; i < n; i++)
        {
            div = div + UNIT_DIV;
        }
        std::cout << div << std::endl;
    };
    
    inline void dashed_divider()
    {
        std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
    };

    template<typename T>
    inline void print(T x)
    {
        std::cout << std::boolalpha << std::left << std::setprecision(9);  
        std::cout << std::setw(PRINT_SPACING) << x << std::endl;
    };

    template <typename First, typename... Rest>
    inline void print(First first, Rest... rest)
    {
        std::cout << std::boolalpha << std::left << std::setprecision(9);  
        std::cout << std::setw(PRINT_SPACING) << first;
        print(rest...);
    } 

    template<typename T>
    inline void print(std::vector<T> v)
    {
        std::cout << std::boolalpha << std::setprecision(9);  
        for (auto vi : v)
        {
            std::cout << std::left << std::setw(PRINT_SPACING) << vi << std::endl;
        };
    };

    // ---------------------------------------------------------------------------
    // String operations

    // Produce a string with the format "name = value units"

    template <typename T>
    inline std::string var_def(std::string name, T value, std::string units = "")
    {
        std::stringstream ss;
        ss << std::setprecision(STRING_PRECISION) << name + " = " << value << " " + units;
        return ss.str();
    };

    // Print a string centered on the terminal 
    inline void centered(int n, std::string words)
    {
        int x = words.length();
        int gap_width = (n * PRINT_SPACING - x)/2;
        std::cout << std::left << std::setw(gap_width) << "" << std::setw(x) << words << std::setw(gap_width) << "" << std::endl;
    };

    // ---------------------------------------------------------------------------
    // Print N columns of data to file
    // We assume theyre all the same size, if not this breaks and its not our fault

    template<int N>
    inline void print_to_file(std::string outname, std::array<std::vector<double>,N> data)
    {
        std::ofstream out;
        out.open(outname);

        for (int j = 0; j < data[0].size(); j++)
        {
            out << std::left;
            for (int i = 0; i < N; i++)
            {
                out << std::setw(PRINT_SPACING) << data[i][j];
            }
            out << std::endl;
        };

        out.close();
        return;
    };

    template<int N>
    inline void print_to_file(std::string outname, std::array<std::string,N> headers, std::array<std::vector<double>,N> data)
    {
        std::ofstream out;
        out.open(outname);

        out << std::left << std::setw(PRINT_SPACING) << "#" + headers[0];
        for (int i = 1; i < N; i++)
        {
            out << std::setw(PRINT_SPACING) << headers[i]; 
        };
        out << std::endl;

        for (int j = 0; j < data[0].size(); j++)
        {
            out << std::left;
            for (int i = 0; i < N; i++)
            {
                out << std::setw(PRINT_SPACING) << data[i][j];
            }
            out << std::endl;
        };

        out.close();
        return;
    };
};

#endif