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
        polynomial(std::vector<double> coef){
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

        std::vector<complex> find_root(complex guess){    // find a real root or complex conjugate pair
            complex prev(guess.r-1.0, guess.i-1.0);
            complex cur = guess;
            complex d(0.0, 0.0);
            std::vector<complex> roots;

            while(abs(prev.norm()-cur.norm()) > 0.0000001){
                prev = cur;
                d = c_deriv_at(cur);
                if(d.r == 0 and d.i == 0) cur = cur + 0.01;
                cur = cur - c_at(cur)/d;
            };
            roots.push_back(cur);
            if (abs(cur.i) > 0.00001) roots.push_back(cur.conj());
            return roots;
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
    polynomial beep(meep);
    cout<<beep.at(1.0)<<endl;
    cout<<beep.deriv_at(1.0)<<endl;
    cout<<beep.find_real_root(1.0)<<endl;

    std::vector<double> moop = {-1.0, 0.0, 0.0, 1.0};
    polynomial boop(moop);
    complex guess(-1.0, -1.0);
    std::vector<complex> roots = boop.find_root(guess);
    cout<<roots[0];
    cout<<roots[1];
};
