#include<cmath>
#include<iostream>
#include<vector>
using namespace std;

class complex{

    public:
        complex(){
            r = i = 0.0;
        };

        complex(double real, double imag){
            r = real;
            i = imag;
        };

        double norm(){
            return sqrt(r*r + i*i);
        };

        complex pow(int n){
            complex result(r, i);
            if (n==0){ 
                result.r = 1.0;
                result.i = 0.0;
            }else{
                 for(int j=1; j<n; j++) result = result*(*this);
            };
            return result;
        };

        complex conj(){
            complex c(r, -1.0*i);
            return c;
        };

        complex operator+(complex const& c){
            complex res(r + c.r, i+c.i);
            return res;
        };

        complex operator+(double x){
            complex res(r + x, i);
            return res;
        };

        complex operator-(complex const& c){
            complex res(r - c.r, i-c.i);
            return res;
        };

        complex operator*(complex const& c){
            complex res(r*c.r - i*c.i, r*c.i + i*c.r);
            return res;
        };

        complex operator*(double a){
            complex res(r*a, i*a);
            return res;
        };

        complex operator/(complex c){
            complex res(r*c.r + i*c.i, i*c.r - r*c.i);
            double n = std::pow(c.norm(),2);
            res.i = res.i/n;
            res.r = res.r/n;
            return res;
        };

        friend ostream& operator<<(ostream &out, const complex &c){
            if (c.i < 0.0) out << c.r << "-i" << abs(c.i);
            else if (c.i == 0.0) out << c.r;
            else out << c.r << "+i" << c.i;
            return out;
        };

        double r, i;
    private:
};


class polynomial{

    public:
        polynomial(vector<double> coef){
            coefficients = coef;
            deg = coef.size() - 1;
        };

        double find_real_root(double guess){    // find real root from starting point via newtons method
            double prev = guess-1.0;
            double cur = guess;

            while(abs(prev-cur) > 0.0000001){
                prev = cur;
                if (deriv_at(cur) == 0.0) cur += 0.00001; 
                cur -= at(cur)/deriv_at(cur);
            };

            return cur;
        };

        complex find_root(complex guess){    // find a real root or complex conjugate pair
            complex prev(guess.r-1.0, guess.i-1.0);
            complex cur = guess;

            while(abs(prev.norm()-cur.norm()) > 0.0000001){
                prev = cur;
                cur = cur - at(cur)/deriv_at(cur);
            };
            if (abs(cur.i) < 0.00001) cur.i =0.0;
            return cur;
        };

        vector<complex> jt_roots(){        // find all roots via the Jenkins-Traub method

        };

        vector<complex> newton_roots(){    // find all roots via newtons method
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
                mult=1;
                while(check){            // identify multiplicities and adjust num_found
                    d = nth_deriv_at(next, mult);
                    if(d.norm() < 0.000001 ) mult++;
                    else check=false;
                };
                if (next.i==0) num_found += mult;
                else num_found += 2*mult;
            };

            coefficients = temp;
            return roots;
            };

        double at(double x){                    // evalutate at real number
            double result=coefficients[deg];
            for(int i=deg-1; i>=0; i--) result = result*x + coefficients[i];
            return result;
        };

        double deriv_at(double x){              // evaluate deriv at real number
            double result=coefficients[deg]*deg;
            for(int i=deg-1; i>=1; i--) result = result*x + coefficients[i]*i;
            return result;
        };

        complex at(complex z){                // evaluate at complex number
            complex result(coefficients[deg], 0.0);
            for(int i=deg-1; i>=0; i--) result = result*z + coefficients[i];
            return result;
        };

        complex deriv_at(complex z){          // evaluate derivative at complex number
            complex result(coefficients[deg]*deg, 0.0);
            for(int i=deg-1; i>=1; i--) result = result*z + coefficients[i]*i;
            return result;
        };

        double nth_deriv_at(double x, int n){

        };

        complex nth_deriv_at(complex z, int n){

        };

        void disp(){
            for(int i=0; i<=deg; i++) cout<<coefficients[i]<<' ';
            cout<<endl;
        }

    private:
        vector<double> coefficients;
        int deg;
};


int main(){
    vector<double> meep = {-2.0, 0.0, 1.0};
    polynomial beep(meep);
    cout<<beep.at(1.0)<<endl;
    cout<<beep.deriv_at(1.0)<<endl;
    cout<<beep.find_real_root(1.0)<<endl;

    vector<double> moop = {-1.0, 0.0, 0.0, 1.0};
    polynomial boop(moop);
    complex guess(-1.0, -1.0);
    complex root = boop.find_root(guess);
    cout<<root<<endl;
    vector<complex> roots = boop.newton_roots();
    cout<<'('<<roots[0]<<", "<<roots[1]<<", "<<roots[2]<<')'<<endl;

};
