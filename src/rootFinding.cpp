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

#include "rootFinding.h"
#include <iostream>

const int __gsl_exception = 0x33627 ;

void gsl_handler(const char * reason, const char * file, int line,  int gsl_errno){
    printf("Some error happend at GSL\n");
    throw __gsl_exception;
}

double polynomial(double x, void *params){
    return polynomial(x, (struct PolyParams *) params);
}

double polynomial(double x, struct PolyParams *params){
    return polynomial(x, params->coeff ) ;
}

double polynomial(double x, std::vector<double> coeff){
    double result = 0 ;
    for(int i=0;i< coeff.size(); i++){
        result += coeff[i]*pow(x, i);
    }
    return result ;
}


double polynomial_derivative(double x, void *params){
    return polynomial_derivative(x, (struct PolyParams *) params);
}

double polynomial_derivative(double x, struct PolyParams *params){
    return polynomial_derivative(x, params->coeff ) ;
}

double polynomial_derivative(double x, std::vector<double> coeff){
    double result=0;
    for(int i=1;i<coeff.size();i++){
        result += i * coeff[i] * pow(x, i-1) ;  //derivative of polynomial
    }
    return result ;
}

void polynomial_fdf (double x, void *params, double *y, double *dy){
    struct PolyParams *p  = (struct PolyParams *) params;    
    for(int i=0;i< p->coeff.size(); i++){
        *y  += p->coeff[i]*pow(x, i);
        *dy += i * p->coeff[i] * pow(x, i-1);
    }
}

Solution findRoot(std::vector<double> coeff, double min, double max){
    
    Solution solution ;
    solution.existance = false ;
    
//    gsl_set_error_handler (&gsl_handler);
    
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0, r_expected = (min+max)/2;
    double x_lo = min, x_hi = max;
    gsl_function F;
    
    struct PolyParams params;
    params.coeff = coeff;
    
    F.function = &polynomial;
    F.params = &params;
    
    //T = gsl_root_fsolver_brent;
    T = gsl_root_fsolver_falsepos;
    s = gsl_root_fsolver_alloc (T);
    try{
        gsl_root_fsolver_set (s, &F, x_lo, x_hi);
        printf ("using %s method\n", gsl_root_fsolver_name (s));
        printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root",  "err", "err(est)");
        
        do{
            iter++;
            status = gsl_root_fsolver_iterate (s);
            r = gsl_root_fsolver_root (s);
            x_lo = gsl_root_fsolver_x_lower (s);
            x_hi = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
            
            if (status == GSL_SUCCESS){
                solution.existance = true ;
                solution.root = r ;
                printf ("Converged:\n");
                printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
                        iter, x_lo, x_hi,
                        r, r - r_expected, 
                        x_hi - x_lo);
            }
        }
        while (status == GSL_CONTINUE && iter < max_iter);
        
        gsl_root_fsolver_free (s);
       return solution ;        

    }catch(int e){
        printf("GSL throws an exception ... %d\n" , e );
        solution.existance=false;
        return solution;
    }

}


Solution findRoot_fdf(std::vector<double> coeff, double min, double max){
    
    Solution solution ;
    solution.existance = false ;
    
    //    gsl_set_error_handler (&gsl_handler);
    
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    double x0, x = (min+max)/2;
    double r = 0, r_expected = (min+max)/2;
    double x_lo = min, x_hi = max;
    gsl_function_fdf FDF;
    
    struct PolyParams params;
    params.coeff = coeff;
    
    FDF.f = &polynomial;
    FDF.df = &polynomial_derivative;
    FDF.fdf = &polynomial_fdf;
    FDF.params = &params;
    
    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc (T);
    

    try{
        gsl_root_fdfsolver_set (s, &FDF, x);
        printf ("using %s method\n", gsl_root_fdfsolver_name (s));
        printf ("%-5s %10s %10s %10s\n", "iter", "root", "err", "err(est)");
        do{
            status = gsl_root_fdfsolver_iterate (s);
            x0 = x;
            x = gsl_root_fdfsolver_root (s);
            status = gsl_root_test_delta (x, x0, 0, 1e-3);
            
            if (status == GSL_SUCCESS){
                //double checking if the result is within the given interval 
                if( min<=x && x<=max){
                    solution.existance = true ;
                    solution.root = x ;                    
                    printf ("Converged:\n");
                    printf ("%5d %10.7f %+10.7f %10.7f\n", iter, x, x - r_expected, x - x0);
                }else{
                    solution.existance = false ;
                    printf("converged, but the root is not within the input interval\n");
                }
            }
        }
        while (status == GSL_CONTINUE && iter < max_iter);
        
        gsl_root_fdfsolver_free (s);
        return solution ;        
        
    }catch(int e){
        printf("GSL throws an exception ... %d\n" , e );
        solution.existance=false;
        return solution;
    }
    
}