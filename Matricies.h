#pragma once
#include<cmath>
#include<iostream>
#include<vector>
#include<complex>

using namespace std;


class matrix
{
public:
    matrix(vector<vector<complex<double>>> M);

    matrix(vector<vector<double>> M);

    vector<complex<double>>& operator[](int i);

    int get_dim();

    void disp();

    vector<complex<double>> solve(vector<complex<double>> b);

    vector<complex<double>> solve(vector<double> b);

    vector<complex<double>> map(vector<complex<double>> v);

    vector<complex<double>> max_eigen();

    vector<complex<double>> all_eigen_vals();

    vector<vector<complex<double>>> all_eigen_vec();

    void transform();

    matrix operator*(matrix B);

    matrix operator*(complex<double> z);

    matrix operator+(matrix B);

    matrix operator-(matrix B);

private:

    vector<complex<double>> l_sqrs(vector<complex<double>> b);

    vector<complex<double>> lin_sys(vector<complex<double>> b);

    void lu_decomp();

    void qr_decomp();

    vector<vector<complex<double>>> A;         // the actual entries of the matrix
    vector<vector<complex<double>>> LU;        // the triangular decomposition of the matrix
    vector<vector<complex<double>>> Q;
    vector<vector<complex<double>>> R;
    vector<int> order;                     // track order of rows rather than actually interchange
    int row_num;
    int col_num;  
    bool square;  
    bool lu_current, qr_current;

};
