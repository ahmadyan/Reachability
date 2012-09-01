/*
 Copyright (c) 2012 Seyed Nematollah Ahmadyan
 All rights reserved.
 
 Developed by: 		Seyed Nematollah Ahmadyan [ahmadyan@gmail.com]
 Center of Reliability and High-Performance Computing,
 Coordinated Science Lab,
 Electrical and Computer Engineering Department,
 University of Illinois at Urbana-Champaign
 
 http://netfiles.uiuc.edu/ahmadya2/www
 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal with the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 
 Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
 Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimers in the documentation and/or other materials provided with the distribution.
 Neither the names of <Name of Development Group, Name of Institution>, nor the names of its contributors may be used to endorse or promote products derived from this Software without specific prior written permission.
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
 */


#pragma once

#include <string>
#include <iostream>
#include <vector>
#include "system.h"
using namespace std;

#ifdef _WIN32
	const string gnuPlot = "C:\\opt\\gnuplot\\bin\\gnuplot.exe -persist" ;
#elif  __linux__
	const string gnuPlot = "/usr/bin/gnuplot" ;
#elif  __APPLE__
	//const string gnuPlot = "/Applications/Gnuplot.app/Contents/Resources/bin/gnuplot" ; //Old gnuplot path (v4.0) on my Mac, does not support drawing a polygin, on aqua.
	const string gnuPlot = "/usr/local/bin/gnuplot"; //current gnuPlot path, v4.6 on x11
#else
	const string gnuPlot = "/usr/local/bin/gnuplot";
#error Undefined operating system/Unknown gnuplot path
#endif

//


class Plotter{
    int dim;
    string gnuPlotPath ;
    FILE *gnuplotPipe;
    bool closed;
   
    //EmptyPlots are called in the constructor
    void emptyPlot(string title, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int angleX, int angleY, string xlabel, string ylabel,string zlabel);
    void emptyPlot(string title, double xmin, double xmax, double ymin, double ymax );
    
public:
    //creates empty 2D or 3D plot, the path to gnuPlot must be provided
    Plotter(string path, reachability::System*);
    ~Plotter();
    void execute(string str);
    
    void drawArray(vector<pair<double, double> > trace);
    void drawLine(double iFromX, double iFromY, double iToX, double iToY);
    void drawLine(double iFromX, double iFromY, double iFromZ, double iToX, double iToY, double iToZ);
    void saveToPdf(string path);
    void waitForKey();
    void close();
};

