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
    int countMax = 249  ;
    int polyCount = 0;
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
    	count++;
    	polyCount += 4;
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
            for(int i=0; i<neighbors.size(); i++){
                if(neighbors[i]->isReachable){
                    if( system->isReachable(neighbors[i], node) ){
                        reachable = true ;
                    }
                        
                }
            }
            node->isReachable = reachable;
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
                for(int i=0;i<kids.size();i++){
                    q.push(kids[i]);
                }
            }
        }
        if(true){
        cout << "some statistics about SPT:" << endl ;
        vector<Node*> *v = tree->getNodes();

        double levelMin = 1000 ;
        double levelMax = 0    ;
        double levelTotal = 0 ;
        double volumeMin = 10000;
        double volumeMax = 0    ;
        double volumeTotal = 0 ;
        double leafVolumeMin = 10000 ;
        double leafVolumeMax = 0 ;
        double leafVolumeTotal = 0 ;
        double reachableLeafVolumeMin = 10000 ;
        double reachableLeafVolumeMax = 0 ;
        double reachableLeafVolumeTotal = 0 ;
        int leaf=0;
        int reachableLeaf=0;
        int boundryPoly=0;
        double boundryReachableLeafVolumeMin = 10000 ;
        double boundryReachableLeafVolumeMax = 0 ;
        double boundryReachableLeafVolumeTotal = 0 ;
        for(int i=0;i<v->size();i++){
        	if( v->at(i)->isDivided() == false ){//this is a leaf node.
        		leaf++;
        		leafVolumeTotal += ((Polytope*)(v->at(i)))->getVolume() ;
        		if(leafVolumeMin>((Polytope*)(v->at(i)))->getVolume()) leafVolumeMin = ((Polytope*)(v->at(i)))->getVolume() ;
        		if(leafVolumeMax<((Polytope*)(v->at(i)))->getVolume()) leafVolumeMax = ((Polytope*)(v->at(i)))->getVolume() ;
        		if( v->at(i)->isReachable == true){ //this is a reachable leaf node
        			reachableLeaf++;
        			reachableLeafVolumeTotal += ((Polytope*)(v->at(i)))->getVolume() ;
        			if(reachableLeafVolumeMin>((Polytope*)(v->at(i)))->getVolume()) reachableLeafVolumeMin = ((Polytope*)(v->at(i)))->getVolume() ;
        			if(reachableLeafVolumeMax<((Polytope*)(v->at(i)))->getVolume()) reachableLeafVolumeMax = ((Polytope*)(v->at(i)))->getVolume() ;
        			vector<Node*> c = v->at(i)->getNeighbors(v->at(i));
        			bool isAtTheBoundry=false;
        			for(int j=0;j<c.size(); j++  ){
        				if (c[j]->isReachable ==false){ //this is a reachable polytope at the boundry
        					isAtTheBoundry=true;

        				}
        			}
        			if(isAtTheBoundry){
        				boundryPoly++;
        				boundryReachableLeafVolumeTotal += ((Polytope*)(v->at(i)))->getVolume() ;
        				if(boundryReachableLeafVolumeMin>((Polytope*)(v->at(i)))->getVolume()) boundryReachableLeafVolumeMin = ((Polytope*)(v->at(i)))->getVolume() ;
        				if(boundryReachableLeafVolumeMax<((Polytope*)(v->at(i)))->getVolume()) boundryReachableLeafVolumeMax = ((Polytope*)(v->at(i)))->getVolume() ;
        			}
        		}
        	}

        	levelTotal += v->at(i)->level;
        	if(levelMin>v->at(i)->level) levelMin = v->at(i)->level ;
        	if(levelMax<v->at(i)->level) levelMax = v->at(i)->level ;
        	volumeTotal += ((Polytope*)(v->at(i)))->getVolume() ;
        	if(volumeMin>((Polytope*)(v->at(i)))->getVolume()) volumeMin = ((Polytope*)(v->at(i)))->getVolume() ;
            if(volumeMax<((Polytope*)(v->at(i)))->getVolume()) volumeMax = ((Polytope*)(v->at(i)))->getVolume() ;

        }

        double levelAverage = levelTotal/v->size();
        double volumeAverage = volumeTotal / v->size();
        double leafVolumeAverage = leafVolumeTotal / leaf;
        double reachableLeafVolumeAverage=reachableLeafVolumeTotal / reachableLeaf ;
        double boundryReachableLeafVolumeAverage=boundryReachableLeafVolumeTotal / boundryPoly ;
        cout << "1. Number of divisions=" << count << endl ;
        cout << "2. Number of hyperplanes=" << count*2 << endl ;
        cout << "3. Number of Polytopes=" << polyCount << endl ;
        cout << "4. Polytope Volume (min)=" << volumeMin << endl ;
        cout << "5. Polytope Volume (max)=" << volumeMax << endl ;
        cout << "6. Polytope Volume (avg)=" << levelAverage << endl ;
        cout << "7. SPT depth (max)=" << levelMax << endl ;
        cout << "8. SPT depth (avg)=" << volumeAverage << endl ;
        cout << "9. Polytopes @boundry=" << 0 << endl ;
        cout << "9. V.size=" << v->size() << endl ;
        cout << "10. Leafs=" << leaf << endl ;
        cout << "11. leaf Volume (min)=" << leafVolumeMin << endl ;
        cout << "12. leaf Volume (max)=" << leafVolumeMax << endl ;
        cout << "13. leaf Volume (avg)=" << leafVolumeAverage << endl ;
        cout << "14. Reachable Leafs=" << reachableLeaf << endl ;
        cout << "15. Reachable leaf Volume (min)=" << reachableLeafVolumeMin << endl ;
        cout << "16. Reachable leaf Volume (max)=" << reachableLeafVolumeMax << endl ;
        cout << "17. Reachable leaf Volume (avg)=" << reachableLeafVolumeAverage << endl ;
        cout << "18. Boundry Polys=" << boundryPoly << endl ;
        cout << "19. Boundry Reachable leaf Volume (min)=" << boundryReachableLeafVolumeMin << endl ;
        cout << "20. Boundry Reachable leaf Volume (max)=" << boundryReachableLeafVolumeMax << endl ;
        cout << "21. Boundry Reachable leaf Volume (avg)=" << boundryReachableLeafVolumeAverage << endl ;
        }
    }
}
