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
#include "system.h"

namespace reachability {
    //    class System;
    class Node{
    public:
        int id;
        static int nodeCount ;
        int dim;       
        Node* parent ;
        bool isReachable ;
        bool containInitialState ;
        bool isInitialState ;
        bool divided;
        
        Node();
        
        static int generateID();
        int getID();
        void setAsInitialState();
        bool isDivided();
        
        virtual std::vector<Node*> getNeighbors(Node*)=0;
        virtual bool contain(Node*)=0;
        virtual void divide(System*)=0;
        virtual std::string toString()=0;
        virtual bool isAdjacent(Node*)=0;        
        virtual std::vector<Node*> getChildren()=0;
        Node* getNodeByID(int _id);
        virtual std::string draw(int)=0;
        virtual std::string draw()=0;
        virtual std::string dump()=0;
    };
}
