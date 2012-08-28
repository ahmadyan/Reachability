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
#include <cmath>
#include <algorithm>
#include "point.h"
#include "geometryUtil.h"
#include "utility.h"
#include "polytope.h"
using namespace std;
namespace reachability {
    Polytope::Polytope(int _dim){
        dim=_dim;
    }
    
    Polytope::Polytope(int _dim, vector<Point*> v){
        dim=_dim;
        points=v;
    }
    
    
    Polytope::~Polytope(){
        for(int i=0;i<points.size();i++) delete points[i] ;
    }
    
    void Polytope::push_back(Point* p){
        points.push_back(p);
    }
    
    int Polytope::getSize(){
        return (int)(points.size());
    }
    
    Point* Polytope::getPoint(int i){
        if(i>=points.size()) throw 0x00001 ;
        else return points[i];
    }
    
    double Polytope::getVolume(){
        if(volume==-1){
            if(dim==2){
                double sum0=0, sum1=0, sum2=0;
                for(int i=0;i<points.size()-1;i++){
                    sum0 += points[i]->getData(0) * points[i+1]->getData(1);
                }
                sum0 += points[points.size()-1]->getData(0) * points[0]->getData(1);
                
                for(int i=0;i<points.size()-1;i++){
                    sum1 += points[i]->getData(1) * points[i+1]->getData(0);
                }
                sum1 += points[points.size()-1]->getData(1) * points[0]->getData(0);
                
                sum2=sum0-sum1;
                sum2/=2;
                volume=abs(sum2);
            }else{
                cout << "Volume undefined for more than 2-dimension" << endl;
            }
        }
        return volume;
    }
    
    Point* Polytope::getCentroid(){
            Point* p = new Point(dim);
            for(int i=0;i<dim;i++){
                double sum=0;        
                for(int j=0;j<points.size();j++){
                    sum += points[j]->getData(i);
                }
                sum /= points.size();
                p->setData(i, sum);
            }
            
        /*    double* min = new double[dim];
            double* max = new double[dim];
            for(int i=0;i<dim;i++){
                min[i]=9999;
                max[i]=-9999;
                
                for(int j=0; j<points.size();j++){
                    double x = points[j]->getData(i);
                    if(x>max[i]) max[i]=x;
                    if(x<min[i]) min[i]=x; 
                }
                p->setData(i, utility::unifRand(min[i], max[i]));
            }
            delete min;
            delete max;    
        */
        
        return p;
    }
    
    int Polytope::getDimension(){
        return dim;
    }
    
    std::vector<Node*> Polytope::getNeighbors(Node*){
        return neighbors;
    }
    
    //This code is adapted from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
    bool Polytope::pointInPoly(Point* p){
        bool result = false;
        int i, j;
        for (i = 0, j = getSize()-1; i < getSize(); j = i++) {
            if ( ((points[i]->getData(1) >p->getData(1) ) != (points[j]->getData(1)>p->getData(1) )) &&
                (p->getData(0)  < (points[j]->getData(0)-points[i]->getData(0)) * (p->getData(1) -points[i]->getData(1) ) / (points[j]->getData(1)-points[i]->getData(1) ) + points[i]->getData(0)) )
                result = !result;
        }
        return result;
    }
    
    //Polytope would contain another polytope if at least one of their edges intersect OR all of their vertex is inside another one.
    bool Polytope::contain(Node* _node){
        bool doContain = false;
        Polytope* node = (Polytope*)_node;
        
        //Test 1: edge intersection
        for(int i=0;i<getSize();i++){
            for(int j=0;j<node->getSize();j++){
                if(geometry::intersection(points[i], points[(i==getSize()-1? 0: i+1)], node->points[j], node->points[(j==node->getSize()-1? 0: j+1)]) ) doContain=true;
            }
        }
        
        
        //Test 2: all vertex inside polytope
        for(int i=0;i<node->getSize();i++){
            if(pointInPoly(node->points[i]) ) doContain=true;
        }
        
        //Test 3: polytope inside the node
        for(int i=0;i<getSize();i++){
            if(node->pointInPoly(points[i]) ) doContain=true;
        }        
        return doContain;        
    }
    
    vector<Point*> Polytope::findIntersectionWithBorders(double* u2, System* system, Point* center){
        vector<Point*> v ;

        double m = u2[1]/u2[0] ;   
        double b = center->getData(1) - m*center->getData(0) ;
        double xmin = system->getMin(0);
        double xmax = system->getMax(0);
        double ymin = system->getMin(1);
        double ymax = system->getMax(1);
        
        double x ;
        double y ;
        
        //             I
        //        .---------.
        //     IV |         |   II
        //        .---------.
        //            II
        
        //I y=max
        y=ymax;
        x = (y-b)/m;
        if(xmin<=x&&xmax>=x){
            Point* p = new Point(dim);
            p->setData(0, x); p->setData(1,y);
            v.push_back(p);
        }
        //II
        x=xmax;
        y=m*x+b;
        if(ymin<=y&&ymax>=y){
            Point* p = new Point(dim);
            p->setData(0, x); p->setData(1,y);
            v.push_back(p);
        }
        
        //III
        y=ymin;
        x = (y-b)/m;
        if(xmin<=x&&xmax>=x){
            Point* p = new Point(dim);
            p->setData(0, x); p->setData(1,y);
            v.push_back(p);
        }
        
        //IV
        x=xmin;
        y=m*x+b;
        if(ymin<=y&&ymax>=y){
            Point* p = new Point(dim);
            p->setData(0, x); p->setData(1,y);
            v.push_back(p);
        }
        
        if(v.size()<2) cout<< "something is wrong at polytope divide function!" << endl ; 
        return v;
    }
        
    void Polytope::removeNeighbor(Polytope* p){
        int index=-1;
        for(int i=0;i<neighbors.size();i++){
            if(neighbors[i]==p) index=i;
        }
        neighbors.erase(neighbors.begin()+index);
    }
    
    void Polytope::addNeighbor(Polytope* p){
        neighbors.push_back(p);
    }
    
    
    //  the update neighbor function has two parts.
    //  first we need to create a list of neighbors for the kids, this is a subset of neighbors from the parent
    //  then we have to update the list of neighbors of the neighbors of the parent.
    void Polytope::updateNeighbors(){
        node1->addNeighbor(node2);
        node2->addNeighbor(node1);
        for(int i=0;i<neighbors.size();i++){
            ((Polytope*)(neighbors[i]))->removeNeighbor(this);
            //Update each of these neighbors and add the new kids
            if( IsSharingAnEdge(((Polytope*)(neighbors[i])), node1) ){
                node1->addNeighbor(((Polytope*)(neighbors[i])));
                ((Polytope*)(neighbors[i]))->addNeighbor(node1);
            }
            if( IsSharingAnEdge(((Polytope*)(neighbors[i])), node2) ){
                node2->addNeighbor(((Polytope*)(neighbors[i])));
                ((Polytope*)(neighbors[i]))->addNeighbor(node2);
            }
        }
    }
    
    pair<bool, vector<Point*> > Polytope::findCommonEdge(Polytope* p1, Polytope* p2){
        //debugging:
        cout << "dumping p1's points:" << endl ;
        for(int i=0;i<p1->points.size();i++){
            cout << p1->points[i]->toString() ;
        }
        
        cout << "dumping p2's points:" << endl ;
        for(int i=0;i<p2->points.size();i++){
            cout << p2->points[i]->toString() ;
        }
        
        
        vector<Point*> v;
        for(int i=0;i<p1->points.size(); i++){
            for(int j=0;j<p2->points.size();j++){
                v.clear();
                Point* a = p1->points[i];
                Point* b = p1->points[((i==p1->points.size()-1)?0:i+1)];
                Point* c = p2->points[j];
                Point* d = p2->points[((j==p2->points.size()-1)?0:j+1)];
                cout << i << " " << ((i==p1->points.size()-1)?0:i+1) << " " << j << " " << ((j==p2->points.size()-1)?0:j+1) << endl ;
                if(geometry::lineSegmentsOverlap(a, b, c, d)){
                    vector<Point*> v ;
                    if( geometry::eq(distance(a,b), distance(a,c)+ distance(c,b))){//c is between a-b
                        if(std::find(v.begin(), v.end(), d)== v.end())
                            v.push_back(c);
                    }
                    if( geometry::eq(distance(a,b), distance(a,d)+ distance(d,b))){// d is between a-b
                        if(std::find(v.begin(), v.end(), d)== v.end())
                            v.push_back(d);
                    }
                    if( geometry::eq(distance(c,d), distance(c,a)+ distance(a,d))){// a is between c-d
                        if(std::find(v.begin(), v.end(), a)== v.end())
                            v.push_back(a);
                    }
                    if( geometry::eq(distance(c,d), distance(c,b)+ distance(b,d))){// b is between c-d
                        if(std::find(v.begin(), v.end(), b)== v.end())
                            v.push_back(b);
                    }
                    cout << "v.size()="  << v.size() << endl ;
                    for(int i=0;i<v.size();i++){
                        cout << v[i]->toString() << endl ;
                    }
                    cout << "-----" << endl ;
                    if(v.size()==2){
                        return make_pair(true, v);
                    }
                }
            }
        }
        if(v.size()!=2){
            return make_pair(false, NULL);
        }
    }
    
    vector<Point*> Polytope::getSharedPoints(Polytope* p1, Polytope* p2){
        pair<bool, vector<Point*> > result = findCommonEdge(p1,p2);
        return result.second;
    }
    
    
    bool Polytope::IsSharingAnEdge(Polytope* p1, Polytope* p2){
        pair<bool, vector<Point*> > result = findCommonEdge(p1,p2);
        return result.first;
    }
        
    //This function is a pain in the ass!
    //TODO: This function is only valid for 2-dimensional systems
    void Polytope::divide(System* system){
        divided=true;
        //First find the center point of the polytope
        Point* center = getCentroid() ;
        
        //Compute direction of vector flow at the center point using system
        double* v1 = system->getVector(center->getData());
        double* u2 = geometry::GramSchmidt(dim, v1) ;
        
        
        //Todo: using points will increase the pointID, we should not use points here. use some other thing in here later.
        vector<Point*> v = findIntersectionWithBorders(u2, system, center);
        
        vector<int> w ;
        for(int i=0;i<points.size()-1;i++){
            if(geometry::intersection(v[0], v[1], points[i], points[i+1])){
                w.push_back(i);
            }
        }
        if(geometry::intersection(v[0], v[1], points[points.size()-1], points[0])){
            w.push_back(points.size()-1);
        }
        
        Point* pi = geometry::intersectionPoint(points[w[0]], (w[0]==points.size()-1?points[0]:points[w[0]+1]), v[0], v[1]);
        Point* pj = geometry::intersectionPoint(points[w[1]], (w[1]==points.size()-1?points[0]:points[w[1]+1]), v[0], v[1]);
                
        
        //Constructing the divided polyople, a and b
        //Collecting points
        //First polytope
        vector<Point*> va ;
        vector<Point*> vb ;
        
        for(int i=0;i<=w[0];i++){
            va.push_back(points[i]);
        }
        va.push_back(pi);
        va.push_back(pj);
        for(int i=w[1]+1; i<points.size();i++){
            va.push_back(points[i]);
        }
        
        vb.push_back(pj);
        vb.push_back(pi);
        
        //2nd Polytope
        for(int i=w[0]+1;i<=w[1];i++){
            vb.push_back(points[i]);
        }
                
        node1 = new Polytope(dim, va);
        node2 = new Polytope(dim, vb);
        divided=true;
        
        updateNeighbors();
        
        //Cleaning up
        delete v1;
        delete u2;
        delete center;
    }
    
    std::string Polytope::toString(){
        stringstream ss ;
        ss << id << " % " ;
        for(int i=0;i<points.size();i++){
            ss << points[i]->toString() << " " ;
        }
        cout << endl ;
        return ss.str(); 
    }
    
    // Checks the neighbors list, If node is inside the neighbors list, return's true.
    bool Polytope::isAdjacent(Node* node){
        for(int i=0;i<neighbors.size();i++){
            if(neighbors[i]==node) return true;
        }
        return false;
    }
    
    vector<Node*> Polytope::getChildren(){
        vector<Node*> v ;
        if(isDivided()){
            v.push_back(node1);
            v.push_back(node2);
        }
        return v;
    }
    
    std::string Polytope::draw(){
        stringstream str ;
        if(isDivided()){
            vector<Node*> kids = getChildren() ;
            for(int i=0;i<kids.size();i++){
                str << kids[i]->draw();
            }
        }else{
            if(isInitialState){
                str << "set object " << id << " polygon from " << points[0]->getData(0) << "," << points[0]->getData(1) ;
                for(int i=1;i<points.size();i++){
                    str << " to " << points[i]->getData(0) << "," << points[i]->getData(1) ;  
                }
                                str << " to " << points[0]->getData(0) << "," << points[0]->getData(1) ;  
                str << endl;
                str << "set object " << id <<" fc rgb \"purple\" fillstyle solid 1.0 border lt -1" << endl ;
            }else if(isReachable){
                
                str << "set object " << id << " polygon from " << points[0]->getData(0) << "," << points[0]->getData(1) ;
                for(int i=1;i<points.size();i++){
                    str << " to " << points[i]->getData(0) << "," << points[i]->getData(1) ;  
                }
                str << " to " << points[0]->getData(0) << "," << points[0]->getData(1) ;  
                str << endl;
                str << "set object " << id <<" fc rgb \"gold\" fillstyle solid 1.0 border lt -1" << endl ;
            }else{
                str << "set object " << id << " polygon from " << points[0]->getData(0) << "," << points[0]->getData(1) ;
                for(int i=1;i<points.size();i++){
                    str << " to " << points[i]->getData(0) << "," << points[i]->getData(1) ;  
                }
                str << " to " << points[0]->getData(0) << "," << points[0]->getData(1) ;  
                str << endl;
                str << "set object " << id <<" fc rgb \"white\" fillstyle solid 1.0 border lt -1" << endl ;
                
            }
        }
        return str.str();   
    }
    
    std::string Polytope::dump(){
        stringstream str ;
        if(isDivided()){
            vector<Node*> kids = getChildren() ;
            for(int i=0;i<kids.size();i++){
                str << kids[i]->dump() << endl ;
            }
        }else{
            str << toString()  << endl ;
        }
        return str.str();    
    }


}
