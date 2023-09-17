// TColor inex definitions for the JPAC color palatte used for plotting
// These are not actually defined here instead 
// they are initialized by the plotter object in plotter.hpp
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef COLORS_HPP
#define COLORS_HPP

#include <TColor.h>

namespace analyticRT
{
    // Give each jpacColor a ROOT TColor index name so it may be called globally after its initialized
    enum class jpacColor: Int_t { Blue   = 2001, Red      = 2002, Green = 2003,
                                  Orange = 2004, Purple   = 2005, Brown = 2006, 
                                  Pink   = 2007, Gold     = 2008, Aqua  = 2009, 
                                  Grey   = 2010, DarkGrey = 2011                };

    // Convert from jpacColor to its underlying int
    inline constexpr Int_t operator+(jpacColor x)
    {
        return static_cast<Int_t>(x);
    };

    constexpr std::array<jpacColor,11> JPACCOLORS = {jpacColor::Blue,   jpacColor::Red,    jpacColor::Green, 
                                                     jpacColor::Orange, jpacColor::Purple, jpacColor::Brown, 
                                                     jpacColor::Pink,   jpacColor::Gold,   jpacColor::Aqua, 
                                                     jpacColor::Grey,   jpacColor::DarkGrey };
};


#endif 