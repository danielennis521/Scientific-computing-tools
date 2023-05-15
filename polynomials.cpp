#include<cmath>
#include<iostream>
#include<vector>
#include "complex.h"
#include "polynomials.h"
using namespace std;


polynomial::polynomial(vector<double> coef){
    coefficients = coef;
    deg = coef.size() - 1;
};

double polynomial::find_real_root(double guess){    // find real root from starting point via newtons method
    double prev = guess-1.0;
    double cur = guess;

    while(abs(prev-cur) > 0.0000000000001){
        prev = cur;
        if (deriv_at(cur) == 0.0) cur += 0.00001; 
        cur -= at(cur)/deriv_at(cur);
    };

    return cur;
};

complex polynomial::find_root(complex guess){    // find a real root or complex conjugate pair
    complex prev(guess.r-1.0, guess.i-1.0);
    complex cur = guess;

    while(abs(prev.norm()-cur.norm()) > 0.0000000000001){
        prev = cur;
        cur = cur - at(cur)/deriv_at(cur);
    };
    if (abs(cur.i) < 0.00000000001) cur.i =0.0;
    return cur;
};

int polynomial::counter_find_root(complex guess){    // find a root and return the root and iterations 
    complex prev(guess.r-1.0, guess.i-1.0);
    complex cur = guess;
    int counter = 0;

    while(abs(prev.norm()-cur.norm()) > 1.0e-20){
        prev = cur;
        cur = cur - at(cur)/deriv_at(cur);
        counter++;
    };
    return counter;
};

/*
vector<complex> jt_roots(){        // find all roots via the Jenkins-Traub method

};
*/

vector<complex> polynomial::newton_roots(){    // find all roots via newtons method
    complex guess, next, d;
    vector<complex> roots;
    vector<double> temp = coefficients;
    vector<double> q;
    double a, b, r, s;
    int num_found=0;
    int mult;
    bool check=false;

    while (num_found <= deg){

        guess.r = rand()%100 -50;
        guess.i = rand()%100 -50;
        next = find_root(guess);
        roots.push_back(next);

        if(next.i == 0){        // monomial synthetic division
            a = next.r;
            r = coefficients[deg-num_found];
            coefficients.pop_back();
            for(int i=deg-1-num_found; i>=0; i--){
                s = coefficients[i];
                coefficients[i] = r;
                r = s + r*a;
            };
        }else{                  // quadratic synthetic division
            a = next.r*next.r + next.i*next.i;
            b = -2.0 * next.r;
            roots.push_back(next.conj());
            for(int i=0; i<deg-2-num_found; i++) q.push_back(0.0);

            for(int i=deg-2-num_found; i>=0; i--){
                q[i] = coefficients[i+2];
                coefficients[i+1] -= q[i]*b;
                coefficients[i] -= q[i]*a;
            };
            coefficients.pop_back();
            coefficients.pop_back();
            for(int i=0; i<=deg-num_found; i++) coefficients[i] = q[i];
        };
        if (next.i==0) num_found ++;
        else num_found += 2;
    };

    coefficients = temp;
    return roots;
};

double polynomial::at(double x){                    // evalutate at real number
    double result=coefficients[deg];
    for(int i=deg-1; i>=0; i--) result = result*x + coefficients[i];
    return result;
};

double polynomial::deriv_at(double x){              // evaluate deriv at real number
    double result=coefficients[deg]*deg;
    for(int i=deg-1; i>=1; i--) result = result*x + coefficients[i]*i;
    return result;
};

complex polynomial::at(complex z){                // evaluate at complex number
    complex result(coefficients[deg], 0.0);
    for(int i=deg-1; i>=0; i--) result = result*z + coefficients[i];
    return result;
};

complex polynomial::deriv_at(complex z){          // evaluate derivative at complex number
    complex result(coefficients[deg]*deg, 0.0);
    for(int i=deg-1; i>=1; i--) result = result*z + coefficients[i]*i;
    return result;
};

double polynomial::nth_deriv_at(double x, int n){
    double c, result;
    c=1.0;
    for(int j=deg; j>deg-n; j--) c*=j;
    result=c*coefficients[deg];
    for(int i=deg-1; i>=n; i--){
        c = 1.0;
        for(int j=i; j>deg-n; j--) c*=j;
        result += result*x + c*coefficients[i];
    };
};

complex polynomial::nth_deriv_at(complex z, int n){

};

polynomial polynomial::operator+(polynomial p){
    vector<double> s = coefficients;
    vector<double> t = p.get_coeff();
    int m = min(deg, p.get_deg());

    for(int i=0; i<=m; i++) s[i] += t[i];
    if (deg > p.get_deg()){
        for(int i=m+1; i<=deg; i++) s.push_back(t[i]);
    };
    polynomial result(s);
    return result;
};

polynomial polynomial::operator+(double x){
    vector<double> s = coefficients;
    s[deg] += x;
    polynomial result(s);
    return result;
};

polynomial polynomial::operator-(polynomial p){
    vector<double> s = coefficients;
    vector<double> t = p.get_coeff();
    int m = min(deg, p.get_deg());

    for(int i=0; i<=m; i++) s[i] -= t[i];
    if (deg > p.get_deg()){
        for(int i=m+1; i<=deg; i++) s.push_back(-1.0*t[i]);
    };
    polynomial result(s);
    return result;
};

polynomial polynomial::operator-(double x){
    vector<double> s = coefficients;
    s[deg] -= x;
    polynomial result(s);
    return result;
};

void polynomial::disp(){
    for(int i=0; i<=deg; i++) cout<<coefficients[i]<<' ';
    cout<<endl;
}

int polynomial::get_deg(){
    return deg;
};

vector<double> polynomial::get_coeff(){
    return coefficients;
};
