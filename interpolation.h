#pragma once
#include<cmath>
#include<iostream>
#include<vector>
#include "polynomials.h"
using namespace std;


polynomial lagrange_coeff(vector<double> x, vector<double> y);

polynomial newton_coeff(vector<double> x, vector<double> y);
