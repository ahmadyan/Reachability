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
#include <vector>
#include <cstdlib>
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
    	for(int i=0;i<nodes.size();i++){
    		nodes[i]->addNeighbor(nodes[(i==nodes.size()-1?0:i+1)]);
            nodes[i]->addNeighbor(nodes[(i==0?nodes.size()-1:i-1)]);
    	}
        for(int i=0;i<neighbors.size();i++){
            ((Polytope*)(neighbors[i]))->removeNeighbor(this);
            //Update each of these neighbors and add the new kids
            for(int j=0;j<nodes.size();j++){
            	if( IsSharingAnEdge(((Polytope*)(neighbors[i])), nodes[j]) ){
            		nodes[j]->addNeighbor(((Polytope*)(neighbors[i])));
            		((Polytope*)(neighbors[i]))->addNeighbor(nodes[j]);
            	}
            }
        }
    }
    
    pair<bool, vector<Point*> > Polytope::findCommonEdge(Polytope* p1, Polytope* p2){
        vector<Point*> v;
        for(int i=0;i<p1->points.size(); i++){
            for(int j=0;j<p2->points.size();j++){
                v.clear();
                Point* a = p1->points[i];
                Point* b = p1->points[((i==p1->points.size()-1)?0:i+1)];
                Point* c = p2->points[j];
                Point* d = p2->points[((j==p2->points.size()-1)?0:j+1)];
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


    //divide will partition a given polytope into 2^n new convex polytopes.
    //The newly generated polytope will be added to the kids nodes collection.
    //The criteria of partitioning is based on the direction of the trajectories at the centeroid of polytope.
    //This function is a pain in the ass!
    //TODO: This function is only valid for 2-dimensional systems, extend this to higher dimensions
    void Polytope::divide(System* system){
        divided=true;

        // Part 1: Find the vertexes of new polytopes.

        //First find the center point of the polytope
        Point* center = getCentroid() ;
        
        //Compute direction of vector flow at the center point using system
        double* trajectoryVector1 = system->getVector(center->getData());
        if(trajectoryVector1[0]==0 && trajectoryVector1[1]==0){
        	trajectoryVector1[0]= 0.27 ; // ((double) rand() / (RAND_MAX+1)) ;
        	trajectoryVector1[1]= 0.86  ;
        }

        //trajectory2 is perpendicular vector to trajectory 1 at center point.
        //keep in mind that there are infinite orthogonal vector to any vector,
        //in GS, we first create the plane by generating a random vector, then return the result.
        double* trajectoryVector2 = geometry::GramSchmidt(dim, trajectoryVector1) ;
        
        //Todo: using points will increase the pointID, we should not use points here. use some other thing in here later.
        //To find intersection points with the polytope, first we need to create a line. Since this line cannot be unlimited
        //and must be bigger than the polytope, we create a line from the edge-to-edge of the state-space with the selected slope.
        //checking where those lines intersects with edges of space,


        vector<Point*> v1 = findIntersectionWithBorders(trajectoryVector1, system, center);
        vector<Point*> v2 = findIntersectionWithBorders(trajectoryVector2, system, center);


        vector< pair<Point*, PointType> > p;
        for(int i=0;i<points.size(); i++){
        	p.push_back(make_pair(points[i], PolyOriginalPoint));
        }
        
        vector<int> w1; //w1 is the array of new intersection points of v1 lines with the polytopes
        vector<int> w2; //w2 is the array of new intersection points of v2 lines with the polytopes
        for(unsigned int i=0;i<points.size()-1;i++){
            if(geometry::intersection(v1[0], v1[1], points[i], points[i+1])){
                w1.push_back(i);
            }
            if(geometry::intersection(v2[0], v2[1], points[i], points[i+1])){
            	w2.push_back(i);
            }
        }
        if(geometry::intersection(v1[0], v1[1], points[points.size()-1], points[0])){
            w1.push_back(points.size()-1);
        }
        if(geometry::intersection(v2[0], v2[1], points[points.size()-1], points[0])){
            w2.push_back(points.size()-1);
        }

        p.push_back(make_pair( geometry::intersectionPoint(points[w1[0]], (w1[0]==points.size()-1?points[0]:points[w1[0]+1]), v1[0], v1[1]) , PolyExtensionPoint));
        p.push_back(make_pair( geometry::intersectionPoint(points[w1[1]], (w1[1]==points.size()-1?points[0]:points[w1[1]+1]), v1[0], v1[1]) , PolyExtensionPoint));
        p.push_back(make_pair( geometry::intersectionPoint(points[w2[0]], (w2[0]==points.size()-1?points[0]:points[w2[0]+1]), v2[0], v2[1]) , PolyExtensionPoint));
        p.push_back(make_pair( geometry::intersectionPoint(points[w2[1]], (w2[1]==points.size()-1?points[0]:points[w2[1]+1]), v2[0], v2[1]) , PolyExtensionPoint));


        //Sort p in clock-wise order based on center point
        for(int i=0;i<p.size();i++){
        	for(int j=0;j<p.size();j++){
        		double angle1 = atan2( p[i].first->getData(1) - center->getData(1), p[i].first->getData(0) - center->getData(0)) ;
        		double angle2 = atan2( p[j].first->getData(1) - center->getData(1), p[j].first->getData(0) - center->getData(0)) ;
        		if(angle1>angle2){
        			swap(p[i], p[j]);
        		}
        	}
        }


        //Constructing the partitioned polytopes,
        //All the points are in sorted in p collection, the idea is that I start with a point from intersection of w1 to one of the faces of parent poly,
        //then I traverse it until I find the point from the intersection of w2 to parent poly. then add center point to this and we have a new winner!
        int i=0;
        while(p[i].second!=PolyExtensionPoint) i++;
        for(int polyCount=0;polyCount< (2<<(dim-1));polyCount++){
        	vector<Point*> poly;
        	poly.push_back(center);
        	poly.push_back(p[i].first);
        	i=(i==p.size()-1)?0:i+1;	//circular-array
        	while(p[i].second!=PolyExtensionPoint){
        		poly.push_back(p[i].first);
        		i=(i==p.size()-1)?0:i+1;	//circular-array
        	}
        	poly.push_back(p[i].first);
        	nodes.push_back( new Polytope(dim, poly) ) ;
        }
        
        updateNeighbors();

        //Cleaning up
        delete trajectoryVector1;
        delete trajectoryVector2;
    }
    
    std::string Polytope::toString(){
        stringstream ss ;
        ss << id << " % " ;
        for(int i=0;i<points.size();i++){
            ss << points[i]->toString() << " " ;
        }
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
        	//upcasting from Polytope to Node class
        	for(int i=0;i<nodes.size();i++){
        		v.push_back(nodes[i]);
        	}
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
