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
#include <queue>
#include <limits>
#include "algorithm.h"

#include "BSPTree.h"

#include "kdtree.h"
#include "hyperbox.h"
namespace reachability{
    int count = 0;
    int countMax = 10 ;
    //Constructor, Sets the system & create the space-partitioning tree data structure
    Algorithm::Algorithm(ReachabilityAlgorithm m,System* _system, double _errorBound){
        mode = m;
        errorBound = _errorBound;
        system = _system;
        currentError= system->getMax(0) - system->getMin(0) ;
        switch (mode) {
            case Alg_hyperbox:
                tree = new KDTree(system);
                break;
            case Alg_polytope:
                tree = new BSPTree(system);
                break;
            default:
                break;
        }
    }
    
    Algorithm::~Algorithm(){}
    
    void Algorithm::divideANode(Node* _node){
        if( mode==Alg_hyperbox ){
            Hyperbox* node = (Hyperbox*) _node ;
            node->divide(system);
            if(node->isInitialState){
                if(node->upLeft->contain(system->getInitialNode())) node->upLeft->setAsInitialState() ;
                if(node->upRight->contain(system->getInitialNode())) node->upRight->setAsInitialState() ;
                if(node->downLeft->contain(system->getInitialNode())) node->downLeft->setAsInitialState() ;
                if(node->downRight->contain(system->getInitialNode())) node->downRight->setAsInitialState() ;
            }
            
            double err = node->upLeft->max[0] - node->upLeft->min[0];
            if( err<currentError ) currentError = err ;
            if( currentError <= errorBound ) node->divided=false ;
            checkReachability(node->upLeft);
            checkReachability(node->upRight);
            checkReachability(node->downLeft);
            checkReachability(node->downRight);    
        }else if( mode==Alg_polytope ){
            Polytope* node = (Polytope*) _node ;
            node->divide(system);
            vector<Node*> kids= node->getChildren();
            if(node->isInitialState){
                for(unsigned int i=0;i<kids.size();i++){
                    if(kids[i]->contain(system->getInitialPoly())) kids[i]->setAsInitialState() ;
                }
            }
        
            double err=numeric_limits<double>::max( );
            for(unsigned int i=0;i<kids.size();i++){
                double volume = ((Polytope*)kids[i])->getVolume();
                if(err<volume) err=volume ;
            }
            if( err<currentError ) currentError = err ;
            if( currentError <= errorBound ) node->divided=false ;

            for(unsigned int i=0;i<kids.size();i++){
                checkReachability(kids[i]);
            }
        }
    }
    
    void Algorithm::checkReachability(Node* _node){
        if(mode==Alg_hyperbox){
            Hyperbox* node = (Hyperbox*) _node ;
            bool reachable = false; 
            if(node->isInitialState) reachable = true ;
            vector<Node*> neighbors = ((KDTree*)tree)->getNeighbors(node);
            for(unsigned int i=0; i<neighbors.size(); i++){
                if(neighbors[i]->isReachable){
                    if( system->isReachable(neighbors[i], node) )
                        reachable = true ;
                }
            }
            node->isReachable = reachable;            
        }else if( mode == Alg_polytope){
            Polytope* node = (Polytope*) _node ;
            bool reachable = false;
            if(node->isInitialState) reachable = true ;
            vector<Node*> neighbors = node->getNeighbors(node);
            cout << "********************************" << endl ;
            cout << _node->toString() << endl ;
            cout << "Checking neighbors for reachability " << neighbors.size() <<  endl ;
            for(int i=0; i<neighbors.size(); i++){
                cout << "Neighbors #" << neighbors[i]->toString() << endl ;
                if(neighbors[i]->isReachable){
                    if( system->isReachable(neighbors[i], node) ){
                        reachable = true ;
                        cout << "Node is reachable from " << neighbors[i]->toString() << endl ;
                    }
                        
                }
            }
            node->isReachable = reachable;
            cout << "Reachability decision: " << reachable << endl;
            cout << "********************************" << endl ;
        }
    }
    
    Tree* Algorithm::getTree(){
        return tree;
    }
    
    System* Algorithm::getSystem(){
        return system;
    }
    
    //Decide if we need to divide this box any further, the decision is based on 
    //error control parameter and the reachability decision.
    bool Algorithm::doesNodeNeedsDivision(Node* _node){
        if(mode==Alg_hyperbox){
            Hyperbox* node = (Hyperbox*) _node ;
            if( currentError <= errorBound ) return false ;
            if( !node->isReachable ) return false;
            
            
            bool nodeIsAdjacentToUnreachableRegion = false; 
            //We assume that the outside of the defined space is unreachable. We need this assumption for the initialization of the algorithm
            if(node->min[0]==system->getMin(0)) nodeIsAdjacentToUnreachableRegion = true ;
            if(node->min[1]==system->getMin(1)) nodeIsAdjacentToUnreachableRegion = true ;    
            if(node->max[0]==system->getMax(0)) nodeIsAdjacentToUnreachableRegion = true ;
            if(node->max[1]==system->getMax(1)) nodeIsAdjacentToUnreachableRegion = true ;
            
            vector<Node*> neighbors = tree->getNeighbors(node);
            for(unsigned int i=0; i<neighbors.size(); i++){
                if(!neighbors[i]->isReachable){
                    nodeIsAdjacentToUnreachableRegion = true ;
                }
            }
            
            if(nodeIsAdjacentToUnreachableRegion) 
                return true ;
            else 
                return false;
   
        }else if(mode==Alg_polytope){
            Polytope* node = (Polytope*) _node ;
            count++;
            if(count>countMax) return false;
            if( currentError <= errorBound ) return false ;
            if( !node->isReachable ) return false;

            bool nodeIsAdjacentToUnreachableRegion = false;
            //We assume that the outside of the defined space is unreachable. We need this assumption for the initialization of the algorithm
            //check if any of polytope points are at the perimeters
            for(int i=0;i<node->getSize();i++){
                    if(node->getPoint(i)->getData(0)==system->getMin(0)  ) nodeIsAdjacentToUnreachableRegion = true ;
                    if(node->getPoint(i)->getData(0)==system->getMax(0)  ) nodeIsAdjacentToUnreachableRegion = true ;
                    if(node->getPoint(i)->getData(1)==system->getMin(1)  ) nodeIsAdjacentToUnreachableRegion = true ;
                    if(node->getPoint(i)->getData(1)==system->getMax(1)  ) nodeIsAdjacentToUnreachableRegion = true ;
            }
            
            for(unsigned int i=0;i<node->getNeighbors(NULL).size(); i++){
                if( !node->getNeighbors(NULL)[i]->isReachable){
                    nodeIsAdjacentToUnreachableRegion = true ;
                }
            }
            if(nodeIsAdjacentToUnreachableRegion)
                return true ;
            else
                return false;
        }
        return false;
    }
    
    //	Main reachability algorithm
    void Algorithm::reachabilityAlgorithm(){
        queue<Node*> q;
        q.push(tree->getRoot());
        while(!q.empty()){
            Node* node = q.front(); q.pop();
            if(doesNodeNeedsDivision(node)){
                //If nodes needs a division, divide it into finer nodes,
                //function divideANode also sets the currentError.
                divideANode(node);
                vector<Node*> kids= node->getChildren();
                
                cout << "Listing kids " <<  kids.size() << endl ;
                for(int i=0;i<kids.size();i++){
                    Polytope* p = (Polytope*)kids[i];
                    cout << p->toString() << endl ;
                    cout << "Listing neighbors for " << p->neighbors.size() << endl ;
                    for(int j=0;j<p->neighbors.size();j++ ){
                        cout << p->neighbors[j]->toString() << endl ;
                    }
                    cout << "----131-2312----123123----" << endl << endl << endl ;
                }
                
                for(int i=0;i<kids.size();i++){
                    q.push(kids[i]);
                }
            }
        }
    }
    
}
