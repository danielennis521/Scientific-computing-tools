#include<cmath>
#include<iostream>
#include<vector>
#include "interpolation.h"


polynomial lagrange_interp(vector<double> x, vector<double> y){
    if (x.size() != y.size()) throw invalid_argument( "vectors must be of equal size" );
    int n=x.size();
    vector<double> u;
    for (int i=0; i<n; i++) u.push_back(0.0);
    double c = 0.0;
    vector<double> v = {0.0, 0.0};
    vector<double> w = {1.0};
    polynomial t(v);
    polynomial s(w);
    polynomial l(w);
    polynomial res(u);

    // need to add code to handle edge case
    for(int i=1; i<n; i++){
        l = s;
        c=x[i] - x[0];
        for(int j=1; j<n; j++) if (i!=j) c *= (x[i] - x[j]); 
        
        for(int j=1; j<n; j++){     // mutiply by monomial factors to get terms of lagrange polynomial
            if (i!=j){

            };
        };

        res = res + l;
    };
};


polynomial newton_interp(vector<double> x, vector<double> y){

    if (x.size() != y.size()) throw invalid_argument( "vectors must be of equal size" );


};