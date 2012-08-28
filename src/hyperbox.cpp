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

#include "hyperbox.h"
#include <iostream>
#include <sstream>
using namespace std;

namespace reachability{
    Hyperbox::Hyperbox(Hyperbox* p, int _dim, double* _min, double* _max){
        parent=p;
        dim = _dim ;
        min = new double[dim];
        max = new double[dim];
        for(int i=0;i<dim;i++){
            min[i] = _min[i];
            max[i] = _max[i];
        }
    }
    
    Hyperbox::~Hyperbox(){
        delete min;
        delete max;
        if (isDivided()) {
            delete upLeft;
            delete upRight;
            delete downLeft;
            delete downRight;
        }
    }
    
    
    string Hyperbox::toString(){
        stringstream str;
        str << " Hyperbox#" << id << " : " ;
        for(int i=0;i<dim;i++){
            str << " [" << min[i] << "," << max[i] << "]" ;
        }
        //if(parent!=NULL)
        //    str << " parent->" << parent->toString() ;
        return str.str();
    }
    
    
    //This function (& the entire class is not scalable, needs revision)
    void Hyperbox::divide(System* sys){
        divided=true;
        
        double *minUpLeft = new double[dim];
        double *minUpRight = new double[dim];
        double *minDownLeft = new double[dim];
        double *minDownRight = new double[dim];
        double *maxUpLeft = new double[dim];
        double *maxUpRight = new double[dim];
        double *maxDownLeft = new double[dim];
        double *maxDownRight = new double[dim];
        
        //  p1=(a,d)      p2=(b,d)
        //  .-------------.
        //  |  UL  |  UR  |
        //  |-------------|
        //  |  DL  |  DR  |
        //  .-------------.
        //  p3=(a,c)      p4=(b,c)
        
        //Warning: dirty code ahead, viewer discresion is advised. 
        minUpLeft[0] = min[0];
        maxUpLeft[0] = (min[0]+max[0])/2;
        minUpLeft[1] = (min[1]+max[1])/2;
        maxUpLeft[1] = max[1];
        
        minUpRight[0]= (min[0]+max[0])/2;
        maxUpRight[0]= max[0];
        minUpRight[1]= (min[1]+max[1])/2;
        maxUpRight[1]= max[1];
        
        minDownLeft[0] = min[0];
        maxDownLeft[0] = (min[0]+max[0])/2;
        minDownLeft[1] = min[1];
        maxDownLeft[1] = (min[1]+max[1])/2;
        
        minDownRight[0] = (min[0]+max[0])/2;
        maxDownRight[0] = max[0];
        minDownRight[1] = min[1];
        maxDownRight[1] = (min[1]+max[1])/2;
        
        upLeft = new Hyperbox(this, dim, minUpLeft, maxUpLeft);
        upRight = new Hyperbox(this, dim, minUpRight, maxUpRight);
        downLeft = new Hyperbox(this, dim, minDownLeft, maxDownLeft); 
        downRight = new Hyperbox(this, dim, minDownRight, maxDownRight);    
    }
    
    //Checks if two Hyperboxs have any common edge between them.
    //This function is used in list neighbor function in Tree class
    bool Hyperbox::isAdjacent(Node* obj){
        Hyperbox* node = (Hyperbox*) obj;
        //checking adjacency for one edge: (if one parameter is the same and the other param is inside the other interval.
        
        //checking p1-p2, p1'-p2' adjacency. //same y & different x
        if( 
           (max[1]==node->min[1]) //same y 
           && ((min[0]<=node->min[0] && max[0]>=node->max[0]) || //decides if Hyperbox is within us or we're bounded within Hyperbox
               (min[0]>=node->min[0] && max[0]<=node->max[0])  ) //different x
           ){
            return true ;
        }
        
        
        if(
           (min[1]==node->max[1]) //same-y
           &&((min[0]<=node->min[0] && max[0]>=node->max[0])||
              (min[0]>=node->min[0] && max[0]<=node->max[0]))
           ){
            return true ;
        }
        
        if(
           (min[0]==node->max[0])
           &&((min[1]<=node->min[1] && max[1]>=node->max[1])||
              (min[1]>=node->min[1] && max[1]<=node->max[1]))
           ){
            return true ;
        }
        
        if(
           (max[0]==node->min[0])
           &&((min[1]<=node->min[1] && max[1]>=node->max[1])||
              (min[1]>=node->min[1] && max[1]<=node->max[1]))
           ){
            return true;
        }
        
        return false;
    }
    
    
    string Hyperbox::draw(){
        stringstream str ;
        if(isDivided()){
            str << upLeft->draw();
            str << upRight->draw();
            str << downLeft->draw();
            str << downRight->draw();
        }else{
            if(isInitialState){
                str << "set object " << id << " rect from " << min[0] << "," << min[1] << " to " << max[0] << "," << max[1] 
                << " fc rgb \"purple\"" << endl ;            
            }else if(isReachable){
                str << "set object " << id << " rect from " << min[0] << "," << min[1] << " to " << max[0] << "," << max[1] 
                << " fc rgb \"gold\"" << endl ;
            }else{
                str << "set object " << id << " rect from " << min[0] << "," << min[1] << " to " << max[0] << "," << max[1] 
                << " fc lt 2" << endl ;
                str << "set obj " << id << " fillstyle empty border -1 front " << endl ;  
            }
        }
        return str.str();
    }
    
    string Hyperbox::dump(){
        stringstream str ;
        str << "Hyperbox# " << id << " from " << min[0] << "," << min[1] << " to " << max[0] << "," << max[1] ;
        if(parent!=0) str << " parent#" << parent->id ;
        str << endl ;
        if(isDivided()){
            str << upLeft->dump();
            str << upRight->dump();
            str << downLeft->dump();
            str << downRight->dump();
        }
        return str.str();
    }
    
    //this code can be further optimized, 
    //we don't need to check all the kids for adjancency, only those that might be adjacent according to their position.
    vector<Node*> Hyperbox::getNeighbors(Node* object){
        Hyperbox* node = (Hyperbox*)object; 
        vector<Node*> v ;
        //Recursively checks for neighboring Hyperboxs
        if(isDivided()){//check for the kids
            vector<Node*> v1 = upLeft->getNeighbors(node);
            vector<Node*> v2 = upRight->getNeighbors(node);
            vector<Node*> v3 = downLeft->getNeighbors(node);
            vector<Node*> v4 = downRight->getNeighbors(node);
            v.insert( v.end(), v1.begin(), v1.end() );
            v.insert( v.end(), v2.begin(), v2.end() );
            v.insert( v.end(), v3.begin(), v3.end() );
            v.insert( v.end(), v4.begin(), v4.end() );
        }else{
            //cout << "checking " << toString() << endl ;
            if( isAdjacent(node) ){
                v.push_back(this);
                //  cout << "***** It is a neighbor" << endl ;
            }else{
                //  cout << "Nope, It is not a neighbor" << endl ;
            }
        }
        return v;
    }
    
    //This function checks whether this Hyperbox contains another Hyperbox. The target Hyperbox is usually the initial set. 
    bool Hyperbox::contain(Node* object){
        Hyperbox* node = (Hyperbox*) object;
        //Hyperbox has four points, for each of them we check whether they are inside or not.
        bool doContain = false; 
        //point1 = Hyperbox->min[0], Hyperbox->min[1]
        //point2 = Hyperbox->min[0], Hyperbox->max[1]
        //point3 = Hyperbox->max[0], Hyperbox->min[1]
        //point4 = Hyperbox->max[0], Hyperbox->max[1]
        if( (min[0]<node->min[0] && node->min[0] < max[0]) && (min[1]< node->min[1] && node->min[1] < max[1]) ){
            doContain = true;
        }  
        if( (min[0]<node->min[0] && node->min[0]  < max[0]) && (min[1]< node->max[1] && node->max[1] < max[1]) ){
            doContain = true;
        }  
        if( (min[0]<node->max[0] && node->max[0] < max[0]) && (min[1]< node->min[1] && node->min[1] < max[1]) ){
            doContain = true;
        } 
        if( (min[0]<node->max[0] && node->max[0] < max[0]) && (min[1]< node->max[1] & node->max[1] < max[1]) ){
            doContain = true;
        }   
        
        if( (node->min[0]<min[0] && min[0] < node->max[0]) && (node->min[1]< min[1] && min[1] < node->max[1]) ){
            doContain = true;
        }  
        if( (node->min[0]<min[0] && min[0]  < node->max[0]) && (node->min[1]< max[1] && max[1] < node->max[1]) ){
            doContain = true;
        }  
        if( (node->min[0]<max[0] && max[0] < node->max[0]) && (node->min[1]< min[1] && min[1] < node->max[1]) ){
            doContain = true;
        } 
        if( (node->min[0]<max[0] && max[0] < node->max[0]) && (node->min[1]< max[1] & max[1] < node->max[1]) ){
            doContain = true;
        }   
        
        return doContain;
    }
    
    vector<Node*> Hyperbox::getChildren(){
        vector<Node*> v ;
        v.push_back(upLeft);
        v.push_back(upRight);
        v.push_back(downRight);
        v.push_back(downLeft);
        return v;
    }
}
