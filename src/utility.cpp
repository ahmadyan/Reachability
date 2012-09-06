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

#include "utility.h"
#include <iostream>
#include <ctime>
#include <cstdlib>
namespace utility {

    using namespace std;
    double initial_time = clock(), final_time=0, total_time=0;
    time_t t1=0, t2=0, t3=0;
    
    //
    // Generate a random number between 0 and 1
    // return a uniform number in [0,1].
    double unifRand(){
        return rand() / double(RAND_MAX);
    }
    //
    // Generate a random number in a real interval.
    // param a one end point of the interval
    // param b the other end of the interval
    // return a inform rand numberin [a,b].
    double unifRand(double a, double b){
        return (b-a)*unifRand() + a;
    }
    //
    // Generate a random integer between 1 and a given value.
    // param n the largest value 
    // return a uniform random value in [1,...,n]
    long unifRand(long n){
        
        if (n < 0) n = -n;
        if (n==0) return 0;
        /* There is a slight error in that this code can produce a return value of n+1
         **
         **  return long(unifRand()*n) + 1;
         */
        //Fixed code
        long guard = (long) (unifRand() * n) +1;
        return (guard > n)? n : guard;
    }
    //
    // Reset the random number generator with the system clock.
    void seed(){
        srand(time(0));
    }

    void tick(){
    	initial_time=clock();
    	time(&t1);
    }

    void tock(){
    	final_time = clock();
    	total_time = (double)(final_time-initial_time) / (double) CLOCKS_PER_SEC ;
    	time(&t2);
    	t3=t2-t1;
    	if(t3>1000) total_time=t3 ;
    	cout << total_time << " seconds " << endl ;
    }

}
