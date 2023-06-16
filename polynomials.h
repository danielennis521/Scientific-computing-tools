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

    polynomial(polynomial &p);

    polynomial();

    void disp();

    int get_deg();

    vector<double> get_coeff();

    double at(double x);

    complex at(complex z);

    double deriv_at(double x);

    complex deriv_at(complex z);

    double nth_deriv_at(double x, int n);

    complex nth_deriv_at(complex& z, int n);

    double find_real_root(double guess);

    complex find_root(complex& guess);

    vector<complex> newton_roots();

    polynomial operator+(polynomial& p);

    polynomial operator+(double x);

    polynomial operator-(polynomial& p);

    polynomial operator-(double x);

    polynomial operator*(polynomial& p);

    polynomial operator*(double x);

    double& operator[](const int i);

    int counter_find_root(complex guess);

};
