#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include "integration.h"

using namespace std;

double integrate(double left, double right, double (*func)(double),
                            double (*Q1)(double, double, double (*f)(double)),
                            double (*Q2)(double, double, double (*f)(double))
                            ){

    // compute the approximation for error checking
    double I = integral_estimate(left, right, func);
    // perform quadrature  
    return adaptive_quadrature(left, right, func, Q1, Q2, I);
    };


double adaptive_quadrature(double left, double right, double (*func)(double),
                            double (*Q1)(double, double, double (*f)(double)),
                            double (*Q2)(double, double, double (*f)(double)),
                            double I){

    double I1 = Q1(left, right, func);
    double I2 = Q2(left, right, func);
    double m = left + (right - left)/2.0;

    if(m <= left || m >= right){
        cout<<"computation terminated because of a lack of machine numbers."<<endl;
        cout<<"the precision condition was not met"<<endl;
        return I2;
    };
    if(I + (I2 - I1) == I) return I2;
    else return (adaptive_quadrature(left, m, func, Q1, Q2, I) 
                + adaptive_quadrature(m, right, func, Q1, Q2, I));

};


double trapezoid_rule(double left, double right, double (*func)(double)){
    return 0.5*(right-left)*(func(left) + func(right));
};


double simpsons_rule(double left, double right, double (*func)(double)){
    return (right-left)*(func(left) + 4.0*func(0.5*(left+right)) + func(right))/6;
};


double gauss2_rule(double left, double right, double (*func)(double)){

};


double fast_integrate(double left, double right, double (*func)(double)){

};


double fast_quadrature(double left, double right, double (*func)(double)){

};


double integral_estimate(double left, double right, double (*func)(double)){

};