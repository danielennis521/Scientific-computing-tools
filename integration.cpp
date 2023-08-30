#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include "integration.h"

using namespace std;


double adaptive_quadrature(double left, double right, double (*func)(double),
                            double (*Q1)(double, double, double (*f)(double)),
                            double (*Q2)(double, double, double (*f)(double))
                            ){

    double I1 = Q1(left, right, func);
    double I2 = Q2(left, right, func);
    double m = left + (right - left)/2.0;
    cout<<"difference: "<<(I2-I1)<<endl;
    cout<<"left: "<<left<<" right: "<<right<<endl;

    if(m <= left || m >= right){
        //cout<<"computation terminated because of a lack of machine numbers."<<endl;
        //cout<<"the precision condition was not met"<<endl;
        return I2;
    } else if(abs(I1 - I2) < 1e-12) {
        return I2;
    }else{
        cout<<"recurse"<<endl; 
        return (adaptive_quadrature(left, m, func, Q1, Q2) 
                + adaptive_quadrature(m, right, func, Q1, Q2));
    };

};


double trapezoid_rule(double left, double right, double (*func)(double)){
    return 0.5*(right-left)*(func(left) + func(right));
};


double midpoint_rule(double left, double right, double (*func)(double)){
    return (right-left)*func((right+left)/2.0);
};


double simpsons13_rule(double left, double right, double (*func)(double)){
    return (right-left)*(func(left) + 4.0*func(0.5*(left+right)) + func(right))/6.0;
};


double simpsons38_rule(double left, double right, double (*func)(double)){
    return (right-left)*(func(left) + 3.0*func((2.0*left + right)/3.0) 
            + 3.0*func((left + 2.0*right)/3.0) + func(right))/8.0;
};


double fast_integrate(double left, double right, double (*func)(double)){

};


double turbo_quadrature(double left, double right, double (*func)(double)){

};


double integral_estimate(double left, double right, double (*func)(double)){

// sample a series of points across the domain of integration
};
