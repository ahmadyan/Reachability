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
#include <vector>
#include <string>

namespace reachability{
    class Node;
    class Polytope;
    class Hyperbox;
    class Point;
    int func (double t, const double y[], double f[], void *params);
    int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);
    
    enum SystemType {vanderpol};
    enum functionSign { positive, negative, posnegative, undefined };
    
    class System{
        int dim ;
        double* min ;
        double* max ;
        Node* init; 
        Polytope* poly;
        double mu;
        
    public:
        System(SystemType);
        
        double getMin(int);
        double getMax(int);
        double* getMinn();
        double* getMaxx();
        int getDimension();
        double* getInitialState();
        Node* getInitialNode();
        Polytope* getInitialPoly();
        double random(double, double);

        
        double* getVector(double* y);
        std::string generateVectorField();
        
        //transient simulation method, for a single state
        std::vector<std::pair<double, double> > simulate();
        
        bool isReachable(Node*, Node*);
        bool isReachable(Polytope* _source, Polytope* _sink);
        bool isReachable(Hyperbox* _source, Hyperbox* _sink);
        functionSign reachability(Point* p1, Point* p2, int axis);
        functionSign reachabilityUsingSampling(Point* p1, Point* p2);
        functionSign fx(double x, double min, double max);
        functionSign fy(double x, double min, double max);
    };
}
