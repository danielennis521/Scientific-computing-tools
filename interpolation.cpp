#include<cmath>
#include<iostream>
#include<vector>
#include "interpolation.h"
#include "polynomials.h"
#include "Matricies.h"


polynomial lagrange_coeff(vector<double> &x, vector<double> &y){
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


polynomial newton_coeff(vector<double> &x, vector<double> &y){
    int n=x.size();
    double p;
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

    vector<double> v = {c[n-2] - c[n-1]*x[n-2], c[n-1]};
    polynomial res(v), t(v);

    t[1] = 1;
    for (int i=n-3; i>=0; i--){
        t[0] = -x[i];
        res = res*t + c[i];
    };

    return res;
};


vector<vector<double>> cubic_spline_interp(vector<double> &x, vector<double> &y){
    if (x.size() != y.size()) throw invalid_argument( "vectors must be of equal size" );
    int dim = 4*(x.size()-1);
    vector<vector<double>> rows, res;
    vector<double> v, t={0.0, 0.0, 0.0, 0.0};
    vector<complex<double>> u;

    for(int i=0; i<dim; i++) v.push_back(0.0);

    for(int i=0; i<(x.size()-1); i++){
        v[4*i] = 1.0;             // left endpoint continuity
        v[4*i+1] =  x[i]; 
        v[4*i+2] = x[i]*x[i];
        v[4*i+3] = v[4*i+2]*x[i];
        rows.push_back(v);

        v[4*i+1] =  x[i+1];        // right endpoint continuity
        v[4*i+2] = x[i+1]*x[i+1];
        v[4*i+3] = v[4*i+2]*x[i+1];
        rows.push_back(v);

        v[4*i] = v[4*i+1] = v[4*i+2] = v[4*i+3] = 0.0;
    };

    for(int i=0; i<(x.size()-2); i++){
        v[4*i+1] = 1.0;           // continuity of first derivative
        v[4*i+2] = 2*x[i+1];
        v[4*i+3] = 3*x[i+1]*x[i+1];

        v[4*i+5] = -1.0;
        v[4*i+6] = -v[4*i+2];
        v[4*i+7] = -v[4*i+3];
        rows.push_back(v);

        v[4*i+1] = v[4*i+5] = 0.0;
        v[4*i+2] = 2.0;           // continuity of second derivative
        v[4*i+3] = 6*x[i+1];

        v[4*i+6] = -2.0;
        v[4*i+7] = -v[4*i+3];
        rows.push_back(v);

        v[4*i+1] = v[4*i+2] = v[4*i+3] = v[4*i+5] = v[4*i+6] = v[4*i+7] = 0.0;
    };

    v[2] = 2.0;                 // natural edge conditions (second derivative zero)
    v[3] = 6*x[0];
    rows.push_back(v);
    v[2] = v[3] = 0.0;

    v[dim-2] = 2.0;
    v[dim-1] = 6*x[x.size()-1];
    rows.push_back(v);
    v[dim-2] = v[dim-1] = 0.0;    

    matrix M(rows);
    M.transform();
    
    v[0] = y[0];                // construct the rhs vector
    for(int i=1; i<y.size()-1; i++){ 
        v[2*i-1] = y[i];
        v[2*i] = y[i];
    };
    v[2*y.size()-3] = y[y.size()-1];

    // M.disp();
    // cout<<endl;
    // for(int i=0; i<v.size(); i++) cout<<v[i]<<' ';
    // cout<<endl;

    u = M.solve(v);
    for(int i=0; i<u.size(); i++) cout<<u[i]<<' ';
    cout<<endl;
    for(int i=0; i<u.size(); i++) v[i] = u[i].real();

    for(int i=0; i<x.size()-1; i++){
        res.push_back({v[i*4], v[i*4+1], v[i*4+2], v[i*4+3]});
    };

    return res;
};


vector<vector<double>> quick_cubic_spline(vector<double> &x, vector<double> &y){

};


double lagrange_interp(double t, vector<double> &x, vector<double> &y){

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



double newton_interp(double t, vector<double> &x, vector<double> &y){
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
