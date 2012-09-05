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


#include <iostream>
#include <vector>
#include <cstdlib>

#include "system.h"
#include "plot.h"
#include "algorithm.h"
#include "tree.h"
#include "polytope.h"
#include "point.h"

const int simulationSampleSize = 2;
const double errorBound=0.1;

int main (void){    
    srand((unsigned)time(NULL));
    //Setting-up system
    reachability::System* system = new reachability::System(reachability::vanderpol);
    Plotter* plotter = new Plotter(gnuPlot, system);
    //reachability::Algorithm* alg = new reachability::Algorithm(reachability::Alg_hyperbox, system, errorBound);
    reachability::Algorithm* alg = new reachability::Algorithm(reachability::Alg_polytope, system, errorBound);
    alg->reachabilityAlgorithm();

    reachability::Tree* tree = (reachability::Tree*)alg->getTree();

    cout << ((reachability::Polytope*)(tree->getNodeByID(11)))->dumpNeighbors() << endl ;
    cout << "isReachable 11,12 " << system->isReachable(tree->getNodeByID(12), tree->getNodeByID(11)) << endl ;
    cout << "isReachable 11,14 " << system->isReachable(tree->getNodeByID(14), tree->getNodeByID(11))<< endl ;

    plotter->execute(tree->drawTree());
    //delete tree;
    

    //Code for drawing a single transient
    for(int i=0;i<simulationSampleSize;i++){
    	vector<pair<double, double> > trace = system->simulate();
    	plotter->drawArray(trace);
    }
    
    plotter->execute(system->generateVectorField());
    plotter->saveToPdf("test.ps");
    plotter->close();
    
    return 0;
}


