#include<cmath>
#include<iostream>
#include<vector>
#include "complex.h"

complex::complex(){
    r = i = 0.0;
};

complex::complex(double real, double imag){
    r = real;
    i = imag;
};

double complex::norm(){
    return sqrt(r*r + i*i);
};

complex complex::pow(int n){
    complex result(r, i);
    if (n==0){ 
        result.r = 1.0;
        result.i = 0.0;
    }else{
            for(int j=1; j<n; j++) result = result*(*this);
    };
    return result;
};

complex complex::conj(){
    complex c(r, -1.0*i);
    return c;
};

complex complex::operator+(complex const& c){
    complex res(r + c.r, i+c.i);
    return res;
};

complex complex::operator+(double x){
    complex res(r + x, i);
    return res;
};

complex complex::operator-(complex const& c){
    complex res(r - c.r, i-c.i);
    return res;
};

complex complex::operator-(double x){
    complex res(r - x, i);
    return res;
};

complex complex::operator*(complex const& c){
    complex res(r*c.r - i*c.i, r*c.i + i*c.r);
    return res;
};

complex complex::operator*(double x){
    complex res(r*x, i*x);
    return res;
};

complex complex::operator/(complex c){
    complex res(r*c.r + i*c.i, i*c.r - r*c.i);
    double n = std::pow(c.norm(),2);
    res.i = res.i/n;
    res.r = res.r/n;
    return res;
};

complex complex::operator/(double x){
    complex res(r/x, i/x);
    return res;
};


std::ostream& operator<<(std::ostream &out, const complex &c){
    if (c.i < 0.0) out << c.r << "-i" << abs(c.i);
    else if (c.i == 0.0) out << c.r;
    else out << c.r << "+i" << c.i;
    return out;
};

        