#include<cmath>
#include<iostream>
#include<vector>
using namespace std;

class complex{

    public:
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
            if (c.i < 0.0) out << c.r << "-i" << abs(c.i) << endl;
            else if (c.i == 0.0) out << c.r << endl;
            else out << c.r << "+i" << c.i << endl;
            return out;
        };

        double r;
        double i;
    private:
};


class polynomial{

    public:
        polynomial(vector<double> coef){
            coefficients = coef;
            deg = coef.size();
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

        vector<complex> find_root(complex guess){    // find a real root or complex conjugate pair
            complex prev(guess.r-1.0, guess.i-1.0);
            complex cur = guess;
            vector<complex> roots;

            while(abs(prev.norm()-cur.norm()) > 0.0000001){
                prev = cur;
                cur = cur - c_at(cur)/c_deriv_at(cur);
            };
            roots.push_back(cur);
            if (abs(cur.i) > 0.00001) roots.push_back(cur.conj());
            else cur.i =0.0;
            return roots;
        };

        vector<complex> jt_roots(){        // find all roots via the Jenkins-Traub method

        };

        vector<complex> newton_roots(){    // find all roots via newtons method
            complex guess(0.0, 0.0);
            vector<complex> next;
            vector<complex> roots;
            vector<double> temp = coefficients;

            while (roots.size() < deg){
                guess.r = rand()%100 -50;
                guess.i = rand()%100 -50;
                next = newton_roots();
                for (int i=0; i<next.size(); i++) roots.push_back(next[i]);

                // todo: perform the synthetic division 

            };
            coefficients = temp;
            return roots;
            };

        double at(double x){                    // evalutate at real number
            double result=coefficients[deg-1];
            for(int i=deg-2; i>=0; i--) result = result*x + coefficients[i];
            return result;
        };

        double deriv_at(double x){              // evaluate deriv at real number
            double result=coefficients[deg-1]*(deg-1.0);
            for(int i=deg-2; i>=1; i--) result = result*x + coefficients[i]*i;
            return result;
        };

        complex c_at(complex z){                // evaluate at complex number
            complex result(coefficients[deg-1], 0.0);
            for(int i=deg-2; i>=0; i--) result = result*z + coefficients[i];
            return result;
        };

        complex c_deriv_at(complex z){          // evaluate derivative at complex number
            complex result(coefficients[deg-1]*(deg-1.0), 0.0);
            for(int i=deg-2; i>=1; i--) result = result*z + coefficients[i]*i;
            return result;
        };

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
    vector<complex> roots = boop.find_root(guess);
    cout<<roots[0];
    cout<<roots[1];

};
