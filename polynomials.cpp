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

        complex pow(int i){

        };

        complex conj(){
            complex c(r, -1*i);
            return c;
        };

        complex operator+(complex const& c){
            complex res(r + c.r, i+c.i);
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
            double n = c.norm();
            res.i = res.i/n;
            res.r = res.r/n;
        };

        friend ostream& operator<<(ostream &out, const complex &c){
            out << c.r << "+i" << c.i << endl;
            return out;
        };

        double r;
        double i;
    private:
};


class polynomial{

    public:
        polynomial(std::vector<double> coef){
            coefficients = coef;

        };

        double find_real_root(double guess){    // find real root from starting point via newtons method
            double prev = guess-1.0;
            double cur = guess;

            while(abs(prev-cur) > 0.000001){
                prev = cur;
                if (deriv_at(cur) == 0.0) cur += 0.00001; 
                cur -= at(cur)/deriv_at(cur);
            };

            return cur;
        };

        complex find_root(complex guess){       // find any root from starting point via newtons method

        };

        std::vector<complex> jt_roots(){        // find all roots via the Jenkins-Traub method

        };

        std::vector<complex> newton_roots(){    // find all roots via newtons method

        };

        double at(double x){                    // evalutate at real number
            double result=0.0;
            for(int i=0; i<deg; i++) result += coefficients[i]*pow(x,i);
            return result;
        };

        double deriv_at(double x){              // evaluate deriv at real number
            double result=0.0;
            for(int i=1; i<deg; i++) result += coefficients[i]*i*pow(x,i-1);
            return result;
        };

        complex c_at(complex z){                // evaluate at complex number
            complex result(0.0, 0.0);
            for(int i=0; i<deg; i++) result = result + z.pow(i) * coefficients[i];
            return result;
        };

        complex c_deriv_at(complex z){          // evaluate derivative at complex number
            complex result(0.0, 0.0);
            for(int i=1; i<deg; i++) result = result + z.pow(i-1) *i*coefficients[i];
            return result;
        };

    private:
        std::vector<double> coefficients;
        int deg;
};


int main(){
    std::vector<double> meep = {-2.0, 0.0, 1.0};
    polynomial boop(meep);
    cout<<boop.at(1.0)<<endl;
    cout<<boop.at(1.0)<<endl;
    cout<<boop.find_real_root(1.0)<<endl;
};