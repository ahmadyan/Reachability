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


#include "geometryUtil.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

namespace geometry{
    
    pair<double, double> solveLine(Point* pi, Point* pj){
        double xi = pi->getData(0);
        double yi = pi->getData(1);
        double xj = pj->getData(0);
        double yj = pj->getData(1);
        double m = ( yj-yi )/( xj-xi );
        double b = yj-m*xj;
        return make_pair(m , b);
    }
    
    //TODO: Scale this to higher dimensions
    double crossProduct(Point* p1, Point* p2){
        if(p1->getDimension()==2){
            //p1*p2 = det[x1 x2; y1 y2] = x1y2-x2y1
            return p1->getData(0)*p2->getData(1) - p2->getData(0)*p1->getData(1);
        }
        return 0;
    }
    
    bool onSegment(Point* pi, Point* pj, Point* pk){
        double xi = pi->getData(0);
        double yi = pi->getData(1);
        double xj = pj->getData(0);
        double yj = pj->getData(1);
        double xk = pk->getData(0);
        double yk = pk->getData(1);
        if(  ((min(xi, xj)<=xk)&&(xk<= max(xi,xj))) && ((min(yi, yj)<=yk)&&(yk<= max(yi,yj))) ) {
            return true;
        }else{
            return false;
        }
    }
    
    //returns (pk-pi) * (pj-pi)
    double direction(Point* pi, Point* pj, Point* pk){
        double x1 = pk->getData(0) - pi->getData(0) ;
        double y1 = pk->getData(1) - pi->getData(1) ;
        double x2 = pj->getData(0) - pi->getData(0) ;
        double y2 = pj->getData(1) - pi->getData(1) ;
        double d  = x1*y2 - x2*y1 ;
        return d;
    }
    
    bool intersection(Point* p1, Point* p2, Point* p3, Point* p4){
        double d1 = direction(p3,p4,p1);
        double d2 = direction(p3,p4,p2);
        double d3 = direction(p1,p2,p3);
        double d4 = direction(p1,p2,p4);
        
        if( 
           (((d1>0)&&(d2<0)) || ((d1<0)&&(d2>0))) && 
           (((d3>0)&&(d4<0)) || ((d3<0)&&(d4>0)))
           ){
            return true;
        }else if(d1==0 && onSegment(p3, p4, p1)){
            return true;
        }else if(d2==0 && onSegment(p3, p4, p2)){
            return true;
        }else if(d3==0 && onSegment(p1, p2, p3)){
            return true;
        }else if(d4==0 && onSegment(p1, p2, p4)){
            return true;
        }else{
            return false;
        }
    }
    
    void unit_test(){
        vector<Point*> points ;
        Point* p0 = new Point(2);   p0->setData(0, 0); p0->setData(1, 0); points.push_back(p0);
        Point* p1 = new Point(2);   p1->setData(0, 0); p1->setData(1, 4); points.push_back(p1);
        Point* p2 = new Point(2);   p2->setData(0, 4); p2->setData(1, 4); points.push_back(p2);
        Point* p3 = new Point(2);   p3->setData(0, 4); p3->setData(1, 0); points.push_back(p3);
        Point* p4 = new Point(2);   p4->setData(0, -4); p4->setData(1, 9); points.push_back(p4);
        
        
        bool b = intersection(p0, p2, p1, p3);
        cout << "Intersection result: " << b << endl ;
        
        b = intersection(p0, p3, p1, p3);
        cout << "Intersection result: " << b << endl ;
        
        b = intersection(p0, p3, p1, p2);
        cout << "Intersection result: " << b << endl ;
        
    }
    
    Point* intersectionPoint(Point* p1, Point* p2, Point* p3, Point* p4) {
        // Store the values for fast access and easy
        // equations-to-code conversion
        float x1 = p1->getData(0), x2 = p2->getData(0), x3 = p3->getData(0), x4 = p4->getData(0);
        float y1 = p1->getData(1), y2 = p2->getData(1), y3 = p3->getData(1), y4 = p4->getData(1);
        
        float d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        // If d is zero, there is no intersection
        if (d == 0){
            cout << "Point intersection, d==0 --> No intersection " << endl ;
            return NULL;   
        }
        
        // Get the x and y
        float pre = (x1*y2 - y1*x2), post = (x3*y4 - y3*x4);
        float x = ( pre * (x3 - x4) - (x1 - x2) * post ) / d;
        float y = ( pre * (y3 - y4) - (y1 - y2) * post ) / d;
        
        // Check if the x and y coordinates are within both lines
       // if ( x < min(x1, x2) || x > max(x1, x2) ||
       //     x < min(x3, x4) || x > max(x3, x4) ) return NULL;
       // if ( y < min(y1, y2) || y > max(y1, y2) ||
       //     y < min(y3, y4) || y > max(y3, y4) ) return NULL;
        // Return the point of intersection
        Point* ret = new Point(2);
        ret->setData(0, x);
        ret->setData(1, y);
        return ret;
    }
    
    double* GramSchmidt(int dim, double* v1){
        double* u1 = new double[dim];
        double* u2 = new double[dim];
        double* v2 = new double[dim];
        double* proj_u1_v2 = new double[dim];
        
        //Todo: this is only for 2d spaces, extend GS for higher dimensions.
        //special cases:
        
        if((v1[0]==0) && (v1[1]==0)){ //equilibrium point, assign a path randomly
            u2[0] = 0.567 ; // ((double) rand() / (RAND_MAX+1)) ;
            u2[1] = 0.46 ;//((double) rand() / (RAND_MAX+1)) ; //0.46 ;
        }else if(v1[0]==0){
            // The vector trajectory is paraxial.
            u2[0] = 1;
            u2[1] = 0;
        }else if(v1[1]==0){ //This is a special case, needs carefull consideration later.
            u2[0] = 0;
            u2[1] = 1;            
        }else{ //Gram–Schmidt process for orthogonalization
            //Find orthogonal vector using gram-schmidt process
            //we have   1) u1=v1
            //          2)  u2=v2-proj_u1(v2), v2 is a random vector
            for(int i=0;i<dim;i++){
                u1[i]=v1[i];
                v2[i]=((double) rand() / (RAND_MAX+1)) ;
                u2[i]=0;
            }
            double u1_dot_v2 = u1[0]*v2[0] + u1[1]*v2[1] ;
            double u1_dot_u1 = u1[0]*u1[0] + u1[1]*u1[1] ;
            for(int i=0;i<dim;i++){
                proj_u1_v2[i] = (u1_dot_v2/u1_dot_u1)*u1[i];
                u2[i] = v2[i] - proj_u1_v2[i];
            }
            
            //Generating a hyperplane orthogonal to the direction of vector flow at center point of polyope            
        }
        delete proj_u1_v2;
        delete u1;
        delete v2;
        return u2;
    }
    
    bool eq(double a, double b){
        if( abs(a-b)<1e-6) return true;
        else return false;
    }
    
    //Wow! there was a bug in my code, apparantly there is a function in std names space that returns the distance between iterators.
    //and gcc took that function instead of using my geometry::distance (without using the namespace geometry).
    //so this caused a lot of problems!
    //Now I changed it to euclideanDistance just to be sure. I hope there is no more function named euclideanDistance in c++ anymore. s
    double euclideanDistance(Point* a, Point* b){
        double x1 = a->getData(0);  
        double y1 = a->getData(1);
        double x2 = b->getData(0);
        double y2 = b->getData(1);
        double d  = sqrt( (y2-y1)*(y2-y1) + (x2-x1)*(x2-x2)) ;
        return d;
    }
    
    
    //For Polytope Adjacency:
    //This is O(n^2) alg, for a more efficient (ON(nlogn):
    //Sort all edges of polygons. In sorting function taking couples of edges, first rule is by direction, and second rule (apllied if directions are the same) is by shift. This way you'll get sets of edges lying on the same line in O(n*log(n)). In such set you can make task one-dimensional and find overlappings in linear (of the number of overlappings) time - just sort ends of line segments and support set of segments owning current position while moving along the line). Resulted complexity is O(n*log(n)+m), n - number of edges, m -number of couples of adjecent edjes. To prevent problems with extremal angles of edges don't divide enything - use vector and scalar products and compare things according to them.
    //From: http://forums.codeguru.com/showthread.php?344401-Determining-polygon-adjacency
   bool lineSegmentsOverlap(Point* a, Point* b, Point* c, Point* d){
        bool pointsAreOnTheSameLine = false;
        double ax = a->getData(0);        double ay = a->getData(1);
        double bx = b->getData(0);        double by = b->getData(1);
        double cx = c->getData(0);        double cy = c->getData(1);
        double dx = d->getData(0);        double dy = d->getData(1);

        //determining line from a to b
        if( (eq(bx,ax))&&eq(cx,dx)&&eq(cx,ax)  ){ //infinity steep
            //checking if c-d has any common segments with a-b or vice versa.
            //if( (min(cy,dy) < min(ay,by) && min(ay,by) < max(cy,dy)) || (min(ay,by)<min(cy,dy) && min(cy,dy)<max(cy,dy)) ){
            if( eq(distance(a,b), distance(a,c) + distance(c,b)) ) pointsAreOnTheSameLine = true;
            if( eq(distance(a,b), distance(a,d) + distance(d,b)) ) pointsAreOnTheSameLine = true;
            if( eq(distance(c,d), distance(c,a) + distance(a,d)) ) pointsAreOnTheSameLine = true;
            if( eq(distance(c,d), distance(c,b) + distance(b,d)) ) pointsAreOnTheSameLine = true;
        }else{
            if(eq(ax, bx)){
                if(ax!=cx) return false ;
            }
            double ab_m = (by-ay)/(bx-ax);
            double ab_b = ay - ab_m*ax;
            //checking if point c and d are on the same line
            if( eq(dy,ab_m*dx+ab_b) && eq(cy,ab_m*cx+ab_b) ){
                //checking if c-d has any common segments with a-b or vice versa.
                // Next, we need to see if the lines actually overlap, which would make the two polygons "adjacent." To do this, I simply use a little point-length test that takes 3 arguments -- the endpoints of one edge in one polygon and a single point from the edge in polygon B as determined by 1-4. If the length from A1 to A2 is equal to the length from A1 to B1 + B1 to A2, then the point is on the line.
                // I then repeat this for all four points, for a total of four "point-on-line" tests... from A1 to A2 testing B1, from A1 to A2 testing B2, from B1 to B2 testing A1 and from B1 to B2 testing A2. (A* are endpoints of an edge from A, likewise for B* points)
                if( eq(distance(a,b), distance(a,c) + distance(c,b)) ) pointsAreOnTheSameLine = true;
                if( eq(distance(a,b), distance(a,d) + distance(d,b)) ) pointsAreOnTheSameLine = true;
                if( eq(distance(c,d), distance(c,a) + distance(a,d)) ) pointsAreOnTheSameLine = true;
                if( eq(distance(c,d), distance(c,b) + distance(b,d)) ) pointsAreOnTheSameLine = true;
                //if(                        
                //   ( (min(cx,dx) < min(ax,bx) && min(ax,bx) < max(cx,dx)) || (min(ax,bx)<min(cx,dx) && min(cx,dx)<max(cx,dx)) ) 
                //   &&
                //   ( (min(cy,dy) < min(ay,by) && min(ay,by) < max(cy,dy)) || (min(ay,by)<min(cy,dy) && min(cy,dy)<max(cy,dy)) ))
                //{
                //    pointsAreOnTheSameLine = true;                            
                // }
            }
        }
    return pointsAreOnTheSameLine;
    }
    
    pair<double,double> rotate(double x, double y, double theta){
        double x_rotated = x*cos(theta) - y*sin(theta) ;
        double y_rotated = x*sin(theta) + y*cos(theta) ;
        return make_pair(x_rotated, y_rotated);
    }
    //Compute the position of the point M w.r.t line thorugh P1-P2
    //It is 0 on the line, and +1 on the left side, -1 on the right side 
    int position(Point* p1, Point* p2, Point* M){
        double Ax = p1->getData(0);
        double Ay = p1->getData(1);
        double Bx = p2->getData(0);
        double By = p2->getData(1);
        double X  = M->getData(0);
        double Y  = M->getData(1);
        double det = (Bx-Ax)*(Y-Ay) - (By-Ay)*(X-Ax) ;
        int pos = (det > 0) - (det < 0) ;
        return pos ;
    }
    
    double rad2deg(double rad){
        return rad * 180.0 / M_PI ;
    }
    
    double deg2rad(double deg){
        return deg * M_PI / 180.0 ;
    }

    //Radially compare point a with point b with respect to point origin.
    //This function is used for sorting points in clock-wise or counter-clockwise order
    //source: http://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
    bool comparePoint(Point* a, Point* b, Point* o){
    	double ax=a->getData(0);
    	double ay=a->getData(1);
    	double bx=b->getData(0);
    	double by=b->getData(1);
    	double ox=o->getData(0);
    	double oy=o->getData(1);
    	if (ax >= 0 && bx < 0)
    		return true;
    	if (eq(ax,0) && eq(bx,0))
    		return ay > by;

    	// compute the cross product of vectors (center -> a) x (center -> b)
    	double det = (ax-ox) * (by-oy) - (bx - ox) * (ay - oy);
    	if (det < 0)
    		return true;
    	if (det > 0)
    		return false;

    	// points a and b are on the same line from the center
    	// check which point is closer to the center
    	double d1 = (ax-ox) * (ax-ox) + (ay-oy) * (ay-oy);
    	double d2 = (bx-ox) * (bx-ox) + (by-oy) * (by-oy);
    	return d1 > d2;
    }
}
