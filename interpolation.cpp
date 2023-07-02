#include<cmath>
#include<iostream>
#include<vector>
#include "interpolation.h"


polynomial lagrange_coeff(vector<double> x, vector<double> y){
    if (x.size() != y.size()) throw invalid_argument( "vectors must be of equal size" );
    int n=x.size();
    vector<double> u;
    for (int i=0; i<n; i++) u.push_back(0.0);
    double c = 1.0;
    vector<double> v = {0.0, 1.0};
    vector<double> w = {1.0};
    polynomial t(v), s(w), l(w), res(u);

    for(int i=0; i<n; i++){
        l = s;
        
        for(int j=0; j<n; j++){     // mutiply by monomial factors to get terms of lagrange polynomial
            if (i!=j){
                t[0] = -1.0 * x[j];
                l = l*t;
            };
        };
        
        c = 1.0;
        for(int j=0; j<n; j++) if (i!=j) c *= (x[i] - x[j]); //scaling factor
        for(int k=0; k<n; k++) l[k] = y[i]*l[k]/c;

        res = res + l;
    };
    return res;
};


polynomial newton_coeff(vector<double> x, vector<double> y){

    if (x.size() != y.size()) throw invalid_argument( "vectors must be of equal size" );


};


double lagrange_interp(double t, vector<double> x, vector<double> y){

    if (x.size() != y.size()) throw invalid_argument( "vectors must be of equal size" );
    int n=x.size();
    double l, c, res = 0.0;

    for(int i=0; i<n; i++){
        l = 1.0;
        c = 1.0;
        for(int j=0; j<n; j++){
             if(i!=j){ 
                l *= (t - x[j]);
                c *= (x[i] - x[j]);
             };
        };
        res += y[i]*l/c;
    };
    return res;
};


double newton_interp(double t, vector<double> x, vector<double> y){
    int n=x.size();
    double res, p;
    vector<double> c;
    if (x.size() != y.size()) throw invalid_argument( "vectors must be of equal size" );

    c.push_back(y[0]);
    for (int i=1; i<n; i++){
        c.push_back(y[i] - c[0]);
        for (int j=1; j<i; j++){
            p=1.0;
            for (int k=0; k<j; k++) p *= (x[i] - x[k]);
            c[i] -= p*c[j];
        };
        p=1.0;
        for (int k=0; k<i; k++) p *= (x[i] - x[k]);
        c[i] /= p;
    };

    res = c[n-2] + c[n-1]*(t-x[n-2]);
    for (int i=n-3; i>=0; i--) res = res*(t - x[i]) + c[i];

    return res;
};
