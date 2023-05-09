#pragma once
#include<cmath>
#include<iostream>
#include<vector>
#include "complex.h"
using namespace std;

class polynomial
{
private:
    vector<double> coefficients;
    int deg;

public:
    polynomial(vector<double> coef);

    void disp();

    double at(double x);

    complex at(complex z);

    double deriv_at(double x);

    complex deriv_at(complex z);

    double nth_deriv_at(double x, int n);

    complex nth_deriv_at(complex z, int n);

    double find_real_root(double guess);

    complex find_root(complex guess);

    vector<complex> newton_roots();

};