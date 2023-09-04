#pragma once
#include<cmath>
#include<iostream>
#include<vector>
#include "polynomials.h"
using namespace std;

// these functions will provide you with the interpolating polynomial given some set of points
polynomial lagrange_coeff(vector<double> &x, vector<double> &y);

polynomial newton_coeff(vector<double> &x, vector<double> &y);

vector<polynomial> cubic_spline_interp(vector<double> &x, vector<double> &y);

// these functions will interpolate the data to the given point "t"
double lagrange_interp(double t, vector<double> &x, vector<double> &y);

double newton_interp(double t, vector<double> &x, vector<double> &y);
