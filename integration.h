#pragma once
#include<cmath>
#include<iostream>
#include<vector>
#include<string>

using namespace std;

double integrate(double left, double right, double (*func)(double),
                            double (*Q1)(double, double, double (*f)(double)),
                            double (*Q2)(double, double, double (*f)(double))
                            );

double adaptive_quadrature(double left, double right, double (*func)(double),
                            double (*Q1)(double, double, double (*f)(double)),
                            double (*Q2)(double, double, double (*f)(double)),
                            double I);

double trapezoid_rule(double left, double right, double (*func)(double));

double simpsons_rule(double left, double right, double (*func)(double));


/* fast_quadrature doesnt allow for the selection of quadrature rules but the 
values of the function that are computed are stored and referenced in order to 
achieve greater efficiency */
double fast_integrate(double left, double right, double (*func)(double));

double fast_quadrature(double left, double right, double (*func)(double));


/* gives the value of an integral within an order of magnitude, used for error 
checking in adaptive quadratures */
double integral_estimate(double left, double right, double (*func)(double));