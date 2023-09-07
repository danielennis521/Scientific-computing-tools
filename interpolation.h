#pragma once
#include<cmath>
#include<iostream>
#include<vector>
#include "polynomials.h"
#include "Matricies.h"
using namespace std;

// these functions will provide you with the interpolating polynomial given some set of points
polynomial lagrange_coeff(vector<double> &x, vector<double> &y);

polynomial newton_coeff(vector<double> &x, vector<double> &y);

vector<vector<double>> cubic_spline_interp(vector<double> &x, vector<double> &y);

vector<vector<double>> nth_spline_interp(vector<double> &x, vector<double> &y);

vector<vector<double>> quick_cubic_spline(vector<double> &x, vector<double> &y);

// these functions will interpolate the data to the given point "t"
double lagrange_interp(double t, vector<double> &x, vector<double> &y);

double newton_interp(double t, vector<double> &x, vector<double> &y);


/* 
this is meant to provide an easy way for you to interact with a cubic spline interpolation
without worrying about the correct interval, set of coefficients etc
*/
class spline_interpolation
{
private:
    vector<double> t;
    vector<polynomial> p;
    int n;

public:
    spline_interpolation(vector<double> &x, int a);

    double at(double x);

    double deriv_at(double x);

    double second_deriv_at(double x);

    void switch_spline_degree(int a);

};
