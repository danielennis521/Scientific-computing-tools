#pragma once
#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include <complex>
using namespace std;

class polynomial
{
private:
    vector<complex<double>> coefficients;
    int deg;

    void mon_div(complex<double> x);

public:

    polynomial(vector<complex<double>> coef);

    polynomial(vector<double> coef);

    polynomial(polynomial& p);

    polynomial();

    void disp();

    int get_deg();

    vector<complex<double>> get_coeff();

    complex<double> at(complex<double> z);

    complex<double> deriv_at(complex<double> z);

    complex<double> nth_deriv_at(complex<double>& z, int n);

    complex<double> find_root(complex<double>& guess);

    complex<double> min_root();

    vector<complex<double>> jt_roots();

    vector<complex<double>> newton_roots();

    polynomial deriv();

    polynomial nth_deriv(int n);

    polynomial operator+(polynomial& p);

    polynomial operator+(complex<double> x);

    polynomial operator-(polynomial& p);

    polynomial operator-(complex<double> x);

    polynomial operator*(polynomial& p);

    polynomial operator*(complex<double> x);

    // vector<polynomial> operator/(polynomial& p);

    polynomial operator/(complex<double> x);

    complex<double>& operator[](const int i);

    void push_head(complex<double> x);

    void push_tail(complex<double> x);

    int counter_find_root(complex<double> guess);

    void pop_tail();

    void pop_head();

};

