#pragma once
#include<cmath>
#include<iostream>
#include<vector>
#include<complex>

using namespace std;

/* The convention Im using here is that each vector in the matrix is a column.
If you have a list of vectors you want to use as rows in a matrix just construct the matrix 
and then use the transpose function */
class matrix
{
public:
    matrix(vector<vector<complex<double>>> M);

    matrix(vector<vector<double>> M);

    vector<complex<double>>& operator[](int i);

    int matrix::get_dim();

    void disp();

    vector<complex<double>> solve(vector<complex<double>> b);

    vector<complex<double>> solve(vector<double> b);

    vector<complex<double>> transform(vector<complex<double>> v);

    vector<complex<double>> max_eigen();

private:

    vector<complex<double>> l_sqrs(vector<complex<double>> b);

    vector<complex<double>> lin_sys(vector<complex<double>> b);

    void lu_decomp();

    vector<complex<double>> qr_decomp();

    vector<vector<complex<double>>> A;         // the actual entries of the matrix
    vector<vector<complex<double>>> LU;        // the triangular decomposition of the matrix
    vector<int> order;                     // track order of rows rather than actually interchange
    int row_dim;
    int col_dim;  
    bool square;  
    bool decomp_current;

};
