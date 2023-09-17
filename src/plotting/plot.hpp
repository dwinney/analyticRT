// Reworking of the jpacStyle as internal library special for jpacPhoto
// This defines a single plot object, which is what actually prints out the files
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef PLOT_HPP
#define PLOT_HPP

#include <array>
#include <vector>
#include <deque>
#include <iostream>   
#include <sstream> 
#include <functional>

#include <TROOT.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLine.h>

#include "data_set.hpp"
#include "elementwise.hpp"
#include "colors.hpp"

namespace analyticRT
{
    class plotter;

    struct entry_style
    {
        jpacColor _color      = jpacColor::DarkGrey;  // Color code 
        int  _style           = 0;                    // Either linestyle or markerstyle code
        bool _add_to_legend   = true;                 // Whether to add this curve to the legend
        std::string _label    = "";                   // Label to add to Legend
        std::string _draw_opt = "L";                // string which enters ROOT::Draw() 
    };

    enum curve_type { 
                      sigma_s,      sigma_w,      sigma_Egam,       // Integrated cross-sections as functions of s, W, and Egam
                      dsigmadt_s,   dsigmadt_w,   dsigmadt_Egam,    // Differential x-sections as function of t at fixed s, W, Egam
                    };

    // Each entry represents a curve to draw as a TGraph
    struct plot_entry 
    {
        // Constructor with just the graph
        plot_entry(TGraph* graph, bool is_data)
        : _graph(graph), _isdata(is_data)
        {};

        // Constructor with an already set up style settings
        plot_entry(TGraph* graph, entry_style xstyle, bool is_data)
        : _graph(graph), _style(xstyle), _isdata(is_data)
        {
            apply_style();
        };

        inline void set_style(entry_style xstyle)
        {
            _style = xstyle; 
            apply_style();
        }

        // Drawing options
        entry_style _style;              

        // Whether this is a data points entry or a curve
        bool _isdata = false; 

        // Graphics object
        TGraph * _graph = NULL;      

        // Apply relevant settings from _style to _graph
        inline void apply_style()
        {
            if (_isdata) _graph->SetLineWidth(_default_markerwidth);
            else         _graph->SetLineWidth(_default_linewidth);
            
            _graph->SetLineColorAlpha(+_style._color, 0.9);
            _graph->SetLineStyle(_style._style);
            _graph->SetMarkerColorAlpha(+_style._color, 0.9);
            _graph->SetMarkerStyle(_style._style);
        };

        // Default line-width
        static constexpr double _default_linewidth   = 3;
        static constexpr double _default_markerwidth = 2;
    };  

    // This class contains the entries, data, and options of producing a single plot/file
    // These can be generated from the plotter->make_plot() method which applies
    // all global settings
    class plot 
    {
        public:

        // Outputting function, generates plot and saves it to file
        void save(); 

        // Adding an optional string will override the saved filename and then save
        inline void save(std::string filename){ _filename = filename; save(); };

        // Add a plot entry
        inline void add_entry(plot_entry curve)
        {
            _entries.push_back(curve);
        };
        
        // -----------------------------------------------------------------------
        // Methods to add data points to your plot

        // Convert a data_set object to a plot_entry
        void add_data(data_set data);

        // This second function can be used if you want the data_set to have
        // a different string id in the legend than the one saved in the data_set
        inline void add_data(data_set data, std::string different_id)
        {
            data_set copy(data);
            copy._id = different_id;
            add_data(copy);
        };

        // Add a small offset to change the running color index
        inline void color_offset(unsigned n)
        {
            _Ncurve += n;
        };

        // -----------------------------------------------------------------------
        // Methods to add curves to your plot

        // Basic function which uses the raw vectors and a string id
        void add_curve(std::vector<double> x, std::vector<double> fx, entry_style style);
        void add_curve(std::vector<double> x, std::vector<double> fx, std::string id = "");
        
        // Take in a lambda an evaluation range to get the vectors
        void add_curve(std::array<double,2> bounds, std::function<double(double)> F, entry_style style);
        void add_curve(std::array<double,2> bounds, std::function<double(double)> F, std::string id = "");

        // Curves added by these functions appear as dashed, not on the legend, and synced with the
        // colors of the "full" curves
        void add_dashed(std::vector<double> x, std::vector<double> fx);
        void add_dashed(std::array<double,2> bounds, std::function<double(double)> F);

        // -----------------------------------------------------------------------
        // Add an error band

        void add_band(std::vector<double> x, std::array<std::vector<double>,2> band, int fill = 1001);

        // -----------------------------------------------------------------------
        // Draw other things such as a vertical line at a certain x value

        inline void add_vertical(double x_val, std::array<int, 2> style = {kBlack, kDashed})
        {
            vertical new_vert;
            new_vert._xvalue    = x_val;
            new_vert._color     = style[0];
            new_vert._linestyle = style[1];   
            _lines.push_back(new_vert);
        };

        // -----------------------------------------------------------------------
        // OPTION SETTERS

        // Add string labels to the axes, follows TLatex 
        inline void set_labels(std::string x, std::string y){ _xlabel = x; _ylabel = y;};
        
        // Set if the x and/or y axes is in logscale
        inline void set_logscale(bool x, bool y){ _xlog = x; _ylog = y; };

        // Set custom bounds for both axes
        inline void set_ranges( std::array<double,2> x,  std::array<double,2> y)
        { _xbounds = x; _ybounds = y; _customranges = true; };

        inline void set_legend(double x, double y){ _addlegend = true; _legendxcord = x; _legendycord = y; };
        inline void set_legend(bool x){_addlegend = x;};
        inline void add_header(std::string x){ _header = x; _addheader = true; };
        inline void add_header(std::string variable, double value, std::string units = "")
        {
            std::ostringstream ss;
            ss << std::setprecision(3) << variable + " " << value << " " + units;
            _header = ss.str(); 
            _addheader = true; 
        };

        inline void add_logo(bool x, std::array<double, 2> coords = {0.93, 0.885}, double scale = 1)
        {
            _add_logo = x; _logo_coords = coords; _logo_scale = scale;
        };
        inline void reset_logo()
        {
            _add_logo = true; _logo_coords =  {0.93, 0.885}; _logo_scale = 1;
        };

        inline void set_legend_spacing(double x)
        {
            _legendyscale = x;
        };

        inline void preliminary(bool x)
        {
            _prelim = x;
        };

        // -----------------------------------------------------------------------
        
        private: 
        
        // Constructor is private, only creatable through plotter
        plot(TCanvas* canvas, std::string file = "")
        : _canvas(canvas), _filename(file)
        {};

        friend class plotter;

        // Change linewidth and propagate it to all the curves in the plot
        inline void scale_linewidth(double x)
        {
            for (auto entry : _entries)
            {
                if (entry._isdata) entry._graph->SetLineWidth(plot_entry::_default_markerwidth * x);
                else               entry._graph->SetLineWidth(plot_entry::_default_linewidth   * x);
            };
        };

        inline void reset_linewidth()
        {
            for (auto entry : _entries)
            {
                if (entry._isdata) entry._graph->SetLineWidth(plot_entry::_default_markerwidth);
                else               entry._graph->SetLineWidth(plot_entry::_default_linewidth);
            };
        };

        // Canvas that the plot actually gets drawn on
        TCanvas* _canvas;

        // Filename of where to produce the desired plot
        std::string _filename;

        bool _add_logo = true;
        std::array<double,2> _logo_coords = {0.93, 0.885};
        double _logo_scale = 1;
        inline void add_logo()
        {
            int red  = +jpacColor::Red;
            int blue = +jpacColor::Blue;
            
            std::string JPAC = "#scale[1.3]{#font[32]{#color[" + std::to_string(blue) + "]{J}}"
                      + "^{#scale[0.8]{#font[32]{" + "#color[" + std::to_string(blue) + "]{P}"
                                                   + "#color[" + std::to_string(red) +  "]{A}"
                                                   + "#color[" + std::to_string(blue) + "]{C}}}}}";

            TLatex *logo = new TLatex(_logo_coords[0], _logo_coords[1], JPAC.c_str());

            logo->SetNDC();
            logo->SetTextSize(2/30. * _logo_scale);
            logo->SetTextAlign(32);
            logo->Draw();
        };

        bool _prelim = false;
        inline void add_watermark()
        {
            TLatex *watermark = new TLatex(0.33,0.15,"PRELIMINARY");
            watermark->SetNDC();
            watermark->SetTextColorAlpha(16, 0.5);
            watermark->SetTextSize(0.1);
            watermark->SetLineWidth(2);
            watermark->Draw();
        };

        // Load all the settings into the graphs and draw them
        void draw();

        // -----------------------------------------------------------------------
        // AXIS SETUP 

        // Whether to use logscale of either axis
        bool _xlog = false, _ylog = false;
            
        // Axis labels
        std::string _xlabel = "", _ylabel = "";
        
        // Custom bounds for the different axes. 
        bool _customranges = false;
        std::array<double,2> _xbounds, _ybounds;

        // -----------------------------------------------------------------------
        // LEGEND SET UP

        bool   _addlegend     = true;
        double _legendxcord   = 0.3, _legendycord   = 0.7;
        double _legendxoffset = 0.3, _legendyoffset = 0.15, _legendyscale = 0.036;

        // Number of entries to expect on the legend
        // used to calculate the offset to be visually appealing
        int _Nlegend = 0; 
        
        bool _addheader = false;
        std::string _header = "";

        // ------------------------------------------------------------------------
        // ENTRY MANAGEMENT

        // Number of entries which are data and which are curves
        // Count used for choosing colors, etc
        int _Ndata = 0, _Ncurve = -1; 

        // When generating curves from functions, evaluate this many points
        int _Npoints = 100;

        // List of all entries (theoretical curves) to be plotted
        std::deque<plot_entry> _entries;

        // ------------------------------------------------------------------------
        // LINE MANAGEMENT

        struct vertical 
        {
            double _xvalue = 0;
            int _color     = kBlack;
            int _linestyle = kDashed;
        };

        std::vector<vertical> _lines;
    };
};

#endif