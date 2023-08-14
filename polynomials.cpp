#include<cmath>
#include<iostream>
#include<vector>
#include <complex>
#include "polynomials.h"
using namespace std;



/*
The following section defines the class methods for the polynomial class
*/

polynomial::polynomial(vector<complex<double>> coef){
    coefficients = coef;
    deg = coef.size() - 1;
};

polynomial::polynomial(vector<double> coef){
    complex<double> c;
    deg = coef.size() - 1;

    for (int i=0; i<=deg; i++){
        c = {coef[i], 0.0};
        coefficients.push_back(c);
    };
};

polynomial::polynomial(polynomial &p){
    vector<complex<double>> coef = p.get_coeff();
    coefficients = coef;
    deg = coef.size()-1;
};

polynomial::polynomial(){
    vector<double> coeff;
    deg = 0;
};

complex<double> polynomial::find_root(complex<double>& guess){    // find a real root or complex conjugate pair
    complex<double> prev(guess.real()-1.0, guess.imag()-1.0);
    complex<double> cur = guess;

    while(abs(abs(prev)-abs(cur)) > 1.0e-20){
        prev = cur;
        cur = cur - at(cur)/deriv_at(cur);
    };
    if (abs(cur.imag()) < 1.0e-12) cur = {cur.real(), 0.0};
    return cur;
};

int polynomial::counter_find_root(complex<double> guess){    // find a root and return the root and iterations 
    complex prev(guess.real()-1.0, guess.imag()-1.0);
    complex cur = guess;
    int counter = 0;

    while(abs(norm(prev)-norm(cur)) > 1.0e-20){
        prev = cur;
        cur = cur - at(cur)/deriv_at(cur);
        counter++;
    };
    return counter;
};


complex<double> polynomial::min_root(){        // finds the minimum modulus root of a polynomial via Jenkins Traub

    int L, M;
    complex<double> c(0.0, 0.0), s(0.0, 0.0), t;
    vector<complex<double>> h;
    
    // stage one 
    // for(int i=1; i<=deg; i++){
    //     c = {coefficients[i], };
    //      h.push_back(i*1.0*c);
    // };

    // for(int i=1; i<5; i++){
    //     t = this->at(s)/horner_eval(s, h);
        
    // };

    // stage two
    for(int i=1; i<L; i++){

    };

    // stage three


};


vector<complex<double>> polynomial::jt_roots(){     // finds all roots via Jenkins Traub algorithm



};


vector<complex<double>> polynomial::newton_roots(){    // find all roots via newtons method
    complex<double> guess, next, d;
    vector<complex<double>> roots;
    vector<complex<double>> temp = coefficients;
    vector<complex<double>> q;
    complex<double> a, b, r, s;
    int num_found=0;
    int mult;
    bool check=false;

    while (num_found <= deg){

        guess =  {rand()%100 -50, rand()%100 -50};
        next = find_root(guess);
        roots.push_back(next);

        if(next.imag() == 0){        // monomial synthetic division
            a = next.real();
            r = coefficients[deg-num_found];
            coefficients.pop_back();
            for(int i=deg-1-num_found; i>=0; i--){
                s = coefficients[i];
                coefficients[i] = r;
                r = s + r*a;
            };
        }else{                  // quadratic synthetic division
            a = next.real()*next.real() + next.imag()*next.imag();
            b = -2.0 * next.real();
            roots.push_back(conj(next));
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
        if (next.imag()==0) num_found ++;
        else num_found += 2;
    };

    coefficients = temp;
    return roots;
};

polynomial polynomial::deriv(){                     // returns the derivative as a new polynomial

    vector<complex<double>> c;
    for(int i=1; i<=deg; i++) c.push_back(i*1.0*coefficients[i]);
    polynomial res(c);
    return res;
};

polynomial polynomial::nth_deriv(int n){                 // returns the nth derivative as a new polynomial


};

complex<double> polynomial::at(complex<double> z){                // evaluate at complex number
    complex<double> result = coefficients[deg];
    for(int i=deg-1; i>=0; i--) result = result*z + coefficients[i];
    return result;
};

complex<double> polynomial::deriv_at(complex<double> z){          // evaluate derivative at complex number
    complex<double> result = deg*1.0*coefficients[deg];
    for(int i=deg-1; i>=1; i--) result = result*z + i*1.0*coefficients[i];
    return result;
};

complex<double> polynomial::nth_deriv_at(complex<double>& z, int n){
    complex<double> result;
    double c=1.0;
    for(int j=deg; j>deg-n; j--) c = c*j;
    result = c*coefficients[deg];
    for(int i=deg-1; i>=n; i--){
        c = 1.0;
        for(int j=i; j>deg-n; j--) c = c*j;
        result = result + result*z + c*coefficients[i];
    };
};

polynomial polynomial::operator+(polynomial& p){
    vector<complex<double>> s = coefficients;
    vector<complex<double>> t = p.get_coeff();
    int m = min(deg, p.get_deg());

    for(int i=0; i<=m; i++) s[i] += t[i];
    if (deg > p.get_deg()){
        for(int i=m+1; i<=deg; i++) s.push_back(t[i]);
    };
    polynomial result(s);
    return result;
};

polynomial polynomial::operator+(double x){
    vector<complex<double>> s = coefficients;
    s[deg] += x;
    polynomial result(s);
    return result;
};

polynomial polynomial::operator-(polynomial& p){
    vector<complex<double>> s = coefficients;
    vector<complex<double>> t = p.get_coeff();
    int m = min(deg, p.get_deg());

    for(int i=0; i<=m; i++) s[i] -= t[i];
    if (deg > p.get_deg()){
        for(int i=m+1; i<=deg; i++) s.push_back(-1.0*t[i]);
    };
    polynomial result(s);
    return result;
};

polynomial polynomial::operator-(double x){
    vector<complex<double>> s = coefficients;
    s[deg] -= x;
    polynomial result(s);
    return result;
};

polynomial polynomial::operator*(polynomial& p){

    vector<complex<double>> v;
    int n=deg;
    int m=p.get_deg();

    for(int i=0; i<= n+m; i++) v.push_back(0.0);

    for(int i=n; i>=0; i--){
        for(int j=m; j>=0; j--){
            v[i+j] += p[j]*coefficients[i];
        };
    };

    polynomial res(v);
    return res;
};

polynomial polynomial::operator*(double x){
    vector<complex<double>> v;
    for (int i=0; i<=deg; i++) v.push_back(coefficients[i]*x);
    polynomial result(v);
    return result;
};

complex<double>& polynomial::operator[](const int i){
    return coefficients[i];
};

void polynomial::disp(){
    for(int i=0; i<=deg; i++) cout<<coefficients[i]<<' ';
    cout<<endl;
}

int polynomial::get_deg(){
    return deg;
};

vector<complex<double>> polynomial::get_coeff(){
    return coefficients;
};

void polynomial::reduce(){                  // resize the polynomial by deleting the constant term
    vector<complex<double>> c;
    for(int i=1; i<=deg; i++) c.push_back(coefficients[i]);
    deg--;
    coefficients = c;
};
