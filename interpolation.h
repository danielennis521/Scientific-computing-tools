#pragma once
#include<cmath>
#include<iostream>
#include<vector>
#include "polynomials.h"
using namespace std;


polynomial lagrange_interp(vector<double> x, vector<double> y);

polynomial newton_interp(vector<double> x, vector<double> y);