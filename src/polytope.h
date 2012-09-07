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
#include "node.h"
#include "point.h"
using namespace std;

namespace reachability {
	enum PointType {PolyOriginalPoint, PolyExtensionPoint};
    class Polytope : public Node{
    public:
        vector<Point*> points ;
        vector<Node*> neighbors ;
        
        double volume ;

        vector<Polytope*> nodes ;

        Polytope(int);
        Polytope(int _dim, vector<Point*> v);
        ~Polytope();
        void push_back(Point*);
        int getSize();
        Point* getPoint(int);
        double getVolume();
        Point* getCentroid();
        int getDimension();
        
        bool pointInPoly(Point* p);
        virtual std::vector<Node*> getNeighbors(Node*);
        virtual bool contain(Node*);
        virtual void divide(System*);
        virtual std::string toString();
        virtual bool isAdjacent(Node*);        
        virtual std::string draw();
        virtual std::string dump();
        virtual std::string draw(int);
        std::vector<Node*> getChildren();
        vector<Point*> findIntersectionWithBorders(double* u2, System* system, Point* center);
        void test();
        pair<bool, vector<Point*> > findCommonEdge(Polytope* p1, Polytope* p2);
        vector<Point*> getSharedPoints(Polytope* p1, Polytope* p2);
        bool IsSharingAnEdge(Polytope* p1, Polytope* p2);
        void removeNeighbor(Polytope* p);
        void addNeighbor(Polytope* p);
        void updateNeighbors();
        std::string dumpNeighbors();
    };
}
