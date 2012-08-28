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


#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
#include <conio.h>   //for getch(), needed in wait_for_key()
#include <windows.h> //for Sleep()
void sleep(int i) { Sleep(i*1000); }
#endif


#define SLEEP_LGTH 2  // sleep time in seconds
#include "system.h"
#include "plot.h"
#include <sstream>
#include <cstdio>

Plotter::Plotter(string path, reachability::System* system){
    gnuPlotPath = path ;
    dim = system->getDimension();
    gnuplotPipe = popen(gnuPlotPath.c_str(),"w");
    string buffer ;
    fflush(gnuplotPipe);
    if(dim==2)
        emptyPlot("test", system->getMin(0), system->getMax(0), system->getMin(1), system->getMax(1)  );
    //else
    //    emptyPlot("test", -5, 5, -5, 5, -5, 5, 0, 0, "a", "b", "c")
    closed=false;
        
}

Plotter::~Plotter(){
    if(!closed) close();
}

void Plotter::close(){
    string buffer = "replot\n";
    fprintf(gnuplotPipe, buffer.c_str());
    fflush(gnuplotPipe);
    
    waitForKey();
    
    pclose(gnuplotPipe);
    closed=true;
}

void Plotter::emptyPlot(string title, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int angleX, int angleY, string xlabel, string ylabel, string zlabel){
    stringstream cmdstr;
    cmdstr << "splot [" << xmin << ":" << xmax << "][" << ymin << ":" << ymax << "][" << zmin << ":" << zmax <<"] 0 with linespoints lt \"white\" pt 0.01" ;
    cmdstr  << " title \"" << title << "\"  \n";
    cmdstr << "set xlabel \"$" << xlabel << "$\" \n";
    cmdstr << "set ylabel \"$" << ylabel << "$\" \n";
    cmdstr << "set zlabel \"$"<< zlabel<< "$\" \n";
    
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}

void Plotter::emptyPlot(string title, double xmin, double xmax, double ymin, double ymax ){
    stringstream cmdstr;
    cmdstr << "plot [" << xmin << ":" << xmax << "][" << ymin << ":" << ymax << "] 0 with linespoints lt \"white\" pt 0.01" ;
    cmdstr  << " title \"" << title << "\"  \n";
    
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}


void Plotter::drawLine(double iFromX, double iFromY, double iFromZ, double iToX, double iToY, double iToZ){
    stringstream cmdstr;
    cmdstr << " set arrow from " << iFromX << "," << iFromY << "," << iFromZ <<  "     to     " << iToX << "," << iToY  << "," << iToZ << "        nohead  lc rgb \"blue\" lw 2 \n" ;
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}

void Plotter::drawLine(double iFromX, double iFromY, double iToX, double iToY){
    stringstream cmdstr;
    cmdstr << " set arrow from " << iFromX << "," << iFromY << " to " << iToX << "," << iToY << " nohead  lc rgb \"blue\" lw 2 \n" ;
    //cout  << " set arrow from " << iFromX << "," << iFromY << " to " << iToX << "," << iToY << " nohead  lc rgb \"blue\" lw 2 " << endl ;
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}


void Plotter::drawArray(vector<pair<double, double> > trace){
    for(int i=0;i<trace.size()-1;i++){
        pair<double,double> point1 = trace[i] ;
        pair<double,double> point2 = trace[i+1] ;
        drawLine( point1.first, point1.second, point2.first, point2.second ) ;
        //cout << point1.first << "," << point1.second <<" ->" << point2.first <<","<< point2.second <<endl;
    }
}

void Plotter::waitForKey(){
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;
    
    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;
    
    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}

void Plotter::saveToPdf(string path){
    stringstream cmdstr;
    cmdstr << " set output \"" << path << "\"\n" ;
    cmdstr << " replot \n" ;
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}

void Plotter::execute(string str){
    fprintf(gnuplotPipe, str.c_str());
    fflush(gnuplotPipe);
}
