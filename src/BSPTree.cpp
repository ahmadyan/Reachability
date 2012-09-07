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
#include <sstream>
#include "BSPTree.h"

namespace reachability {
    BSPTree::BSPTree(System* system){
        //TODO: this only works for two dimension.
        vector<Point*> points ;
        Point* p0 = new Point(2);   p0->setData(0, system->getMin(0)); p0->setData(1, system->getMin(1)); points.push_back(p0);
        Point* p1 = new Point(2);   p1->setData(0, system->getMax(0)); p1->setData(1, system->getMin(1)); points.push_back(p1);
        Point* p2 = new Point(2);   p2->setData(0, system->getMax(0)); p2->setData(1, system->getMax(1)); points.push_back(p2);
        Point* p3 = new Point(2);   p3->setData(0, system->getMin(0)); p3->setData(1, system->getMax(1)); points.push_back(p3);
        root = new Polytope(2, points);
        ((Polytope*)root)->getVolume();
        root->setAsInitialState(); 
    }
    
    BSPTree::~BSPTree(){
    }
    
    int BSPTree::getTreeType(){
        return 2;
    }
    //calling the draw function by id will fill that specific poly in red.
    std::string BSPTree::draw(int x){
        std::stringstream str ;
        str << root->draw(x) << endl ;
        return str.str();
    }

    std::string BSPTree::draw(){
    	return draw(-1);
    }

    std::string BSPTree::dump(){
        std::stringstream str ;
        str << root->dump() << endl ;
        return str.str();
    }
    
    
    
}
