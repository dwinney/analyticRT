// Extend std::vector<double> with elementwise operations.
// This is useful for fast data manipulation
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef ELEMENTWISE_HPP
#define ELEMENTWISE_HPP

#include <vector>
#include <array>
#include "debug.hpp"

namespace analyticRT
{
    
    // ---------------------------------------------------------------------------
    // Element-wise operations on data vectors

    // Given two vector<double>s of the same size, calculate the average element wise
    inline std::vector<double> operator*( std::vector<double> lhs, double c)
    {
        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]*c );
        };
        return result;
    };

    inline std::vector<double> operator*(double c, std::vector<double> rhs)
    {
        std::vector<double> result;
        for (int i = 0; i < rhs.size(); i++)
        {
            result.push_back( c*rhs[i] );
        };
        return result;
    };

    inline std::vector<double> operator/( std::vector<double> lhs, double c)
    {
        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]/c );
        };
        return result;
    };

    inline std::vector<double> operator-(const std::vector<double> & x)
    {
        return -1 * x;
    };

    
    inline std::vector<double> operator+(std::vector<double> lhs, std::vector<double> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i] + rhs[i] );
        };
        return result;
    };

    inline std::vector<double> operator-(std::vector<double> lhs, std::vector<double> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i] - rhs[i] );
        };
        return result;
    };

    // Add a constant to all elements of a vector
    inline std::vector<double> operator+(double lhs, std::vector<double> rhs)
    {
        std::vector<double> result;
        for (int i = 0; i < rhs.size(); i++)
        {
            result.push_back( lhs + rhs[i] );
        };
        return result;
    };   
    inline std::vector<double> operator-(double lhs, std::vector<double> rhs) 
    {
        return lhs + (-rhs);
    };

    inline std::vector<double> operator+(std::vector<double> lhs, double rhs)
    {
        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( rhs + lhs[i] );
        };
        return result;
    };
    inline std::vector<double> operator-(std::vector<double> lhs, double rhs) 
    {
        return lhs + (-rhs);
    };

    inline std::vector<double> multiply_elementwise(std::vector<double> in1, std::vector<double> in2)
    {
        if (in1.size() != in2.size()) warning("multiply_elementwise()", "Input vectors not the same size!");

        std::vector<double> out;
        for (int i = 0; i < in1.size(); i++)
        {
            out.push_back(in1[i]*in2[i]);
        }
        return out;
    };
    inline std::vector<double> square_elementwise(std::vector<double> in){ return multiply_elementwise(in, in); };
};

#endif