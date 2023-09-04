#include<cmath>
#include<iostream>
#include<vector>
#include<complex>
#include<string>
#include<cstdlib>
#include "polynomials.h"
using namespace std;



/*
The following section defines the private methods for the polynomial class
*/

void polynomial::mon_div(complex<double> x){      
    complex<double> s, r;

    r = coefficients[deg];
    for(int i=deg-1; i>=0; i--){      // monomial synthetic division
        s = coefficients[i];
        coefficients[i] = r;
        r = s + r*x;
    };
    coefficients.pop_back();
    deg--;
};


/*
The following section defines the public methods for the polynomial class
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

polynomial::polynomial(polynomial& p){
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

    while(norm(prev - cur) > 1.0e-20){
        prev = cur;
        cur = cur - at(cur)/deriv_at(cur);
    };
    if (abs(cur.imag()) < 1.0e-12) cur = {cur.real(), 0.0};
    return cur;
};

int polynomial::counter_find_root(complex<double> guess){    // find a root and return the root and iterations 
    complex<double> prev(guess.real()-1.0, guess.imag()-1.0);
    complex<double> cur = guess;
    int counter = 0;

    while(abs(norm(prev)-norm(cur)) > 1.0e-20){
        prev = cur;
        cur = cur - at(cur)/deriv_at(cur);
        counter++;
    };
    return counter;
};


complex<double> polynomial::min_root(){        // finds the minimum modulus root of a polynomial via Jenkins Traub
    bool pass=false;
    int reset=0, M=5;
    double r;
    complex<double> s = {0.0, 0.0}, sp = {1.0, 1.0}, pn = this->at(s), t0 = {0.0, 0.0},
                    t1 = {0.0, 0.0}, t2 = {0.0, 0.0};
    polynomial h = this->deriv(), t;
    
    // stage one 
    for(int i=1; i<M; i++){
        t = *this*(h.at(s)/pn);
        h = h-t;
        h.pop_tail();
    };  

    // stage two
    for(int i=0; i<=deg; i++) t[i] = {norm(coefficients[i]/coefficients[deg]), 0.0};
    t[0] = -1.0*t[0];
    s = {1.0, 0.0};
    s = t.find_root(s);     // determine the fixed shift value based on lower bound

    while(!pass){
        r = rand()*2.0*3.141592653589793/RAND_MAX;
        s = {s.real()*cos(r), s.real()*sin(r)};

        for(int i=0; i<3; i++){
            t = *this*(h.at(s)/(this->at(s)));
            h = h-t;
            h.mon_div(s);
            t0 = t1;
            t1 = t2;
            t2 = s - h[get_deg()]*((this->at(s))/h.at(s));   
        };

        while(true){
            t = *this*(h.at(s)/(this->at(s)));
            h = h-t;
            h.mon_div(s);
            t0 = t1;
            t1 = t2;
            t2 = s - h[h.get_deg()]*((this->at(s))/h.at(s));  

            reset++;
            if (reset == 50) break;
            else if(norm(t1 - t0)< 0.5*norm(t0) && norm(t2 - t1)< 0.5*norm(t1)){
                pass = true;
                break;
            };
        };
    };

    // stage three
    s = s - h[0]*((this->at(s))/h.at(s));

    while(norm(s - sp) > 1.0e-12){
        t = *this*(h.at(s)/(this->at(s)));
        h = h-t;
        h.mon_div(s);
        sp = s;
        s = s - h[h.get_deg()]*((this->at(s))/h.at(s));
    };

    return s;
};


vector<complex<double>> polynomial::jt_roots(){     // finds all roots via Jenkins Traub algorithm
    vector<complex<double>> temp = coefficients, roots;
    int d = deg;

    for(int i=0; i<d; i++){ 
        roots.push_back(min_root());
        mon_div(roots[i]);
    };

    deg = d;
    coefficients=temp;
    return roots;
};


vector<complex<double>> polynomial::newton_roots(){    // find all roots via newtons method
    complex<double> guess, next;
    vector<complex<double>> roots, q, temp = coefficients;
    complex<double> a, b, r, s;
    int num_found=0, d=deg;
    bool check=false;

    for(int i=0; i<d; i++){

        guess =  {rand()%100 -50, rand()%100 -50};
        next = find_root(guess);
        roots.push_back(next);
        mon_div(roots[i]);

        num_found++;
    };

    coefficients = temp;
    deg = d;
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

polynomial polynomial::anti_deriv(complex<double> c = {0.0, 0.0}){  // return anti derivative of polynomial with const of integration c
    vector<complex<double>> v;
    v.push_back(c);
    for(int i=0; i<=deg; i++) v.push_back(coefficients[i]/(1.0*(i+1)));
    polynomial res(v);
    return res;   

};

complex<double> polynomial::integrate(complex<double> a, complex<double> b){
    complex<double> Fa=a*coefficients[deg]/(deg-1.0), Fb=b*coefficients[deg]/(deg-1.0);
    for(int i=0; i<=deg; i++){
        Fa = Fa*a + coefficients[i]/(i-1.0);
        Fb = Fb*b + coefficients[i]/(i-1.0);
    };
    return (Fb - Fa);
};

complex<double> polynomial::integrate(double a, double b){
    complex<double> ca={a, 0.0}, cb={b, 0.0};
    return integrate(ca, cb);
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
    return result;
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

polynomial polynomial::operator+(complex<double> x){
    vector<complex<double>> s = coefficients;
    s[deg] += x;
    polynomial result(s);
    return result;
};

polynomial polynomial::operator-(polynomial& p){
    vector<complex<double>> s = coefficients;
    int m = min(deg, p.get_deg());

    for(int i=0; i<=m; i++) s[i] -= p[i];
    if (deg < p.get_deg()) for(int i=m+1; i<=p.get_deg(); i++) s.push_back(-1.0*p[i]);
    polynomial result(s);
    return result;
};

polynomial polynomial::operator-(complex<double> x){
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

polynomial polynomial::operator*(complex<double> x){
    vector<complex<double>> v;
    for (int i=0; i<=deg; i++) v.push_back(coefficients[i]*x);
    polynomial result(v);
    return result;
};

polynomial polynomial::operator/(complex<double> x){
    vector<complex<double>> v;
    for (int i=0; i<=deg; i++) v.push_back(coefficients[i]/x);
    polynomial result(v);
    return result;
};

vector<polynomial> polynomial::operator/(polynomial& p){ // second entry of the returned vector is the remainder from division
    polynomial t, q(coefficients);    
    vector<polynomial> r;
    complex<double> c;

    r.push_back(q);
    r.push_back(t);

    for (int i=deg; i>=p.get_deg(); i--){
        c = ((*this)[i])/p[i];
        r[1].push_tail(c);

        for(int j=0; j<i; j++) r[0][j] = r[0][j] - c*p[j];
    };

    for(int i=deg; i>=p.get_deg(); i--) r[0].pop_head();

    return r;
};

complex<double>& polynomial::operator[](const int i){
    return coefficients[i];
};

void polynomial::disp(){
    for(int i=deg; i>=0; i--){
        cout << coefficients[i];
    };
    cout<<endl;
}

int polynomial::get_deg(){
    return deg;
};

vector<complex<double>> polynomial::get_coeff(){
    return coefficients;
};

void polynomial::pop_tail(){                  // resize the polynomial by deleting the constant term
    vector<complex<double>> c;
    for(int i=1; i<=deg; i++) c.push_back(coefficients[i]);
    deg--;
    coefficients = c;
};

void polynomial::pop_head(){                  // resize the polynomial by deleting the leading term
    coefficients.pop_back();
    deg--;
};

void polynomial::push_head(complex<double> x){
    coefficients.push_back(x);
    deg++;
};

void polynomial::push_tail(complex<double> x){
    vector<complex<double>> c;
    c.push_back(x);
    for(int i=0; i<=deg; i++) c.push_back(coefficients[i]);
    deg++;
    coefficients = c;
};
