#include<cmath>
#include<iostream>
#include<vector>
#include<complex>
#include<string>
#include<cstdlib>
#include "polynomials.h"
using namespace std;


int main(){
    vector<polynomial> vp;
    vector<complex<double>> v1 = {{-1.0, -1.0}, {1.0, 0.0}}, v2 = {{-1.0, -1.0}, {1.0, 0.0}}, v3 = {{-4.0, 3.0}, {1.0, 0.0}},
                            v4 = {{-4.0, -3.0}, {1.0, 0.0}}, v5 = {{-3.999, -3.0}, {1.0, 0.0}};
    polynomial p1(v1), p2(v2), p3(v3), p4(v4), p5(v5), t1, t2;
    complex<double> c;
    vector<complex<double>> n_roots, j_roots;

    t1 = p1*p2;
    t1 = t1*p3;
    t1 = t1*p4;
    t1 = t1*p5;
    t1.disp();

    // cout<<"first division test:"<<endl;
    // vp = t1/p1;
    // vp[0].disp();
    // vp[1].disp();

    // cout<<"second devision test:"<<endl;
    // t2 = p2*p3;
    // t2 = t2*p4;
    // t2 = t2*p5;
    // vp = t1/t2;
    // vp[0].disp();
    // vp[1].disp();

    cout<<"Newtons method root test:"<<endl;
    n_roots = t1.newton_roots();
    for (int i=0; i<5; i++) cout<<n_roots[i];
    cout<<endl;

    cout<<"Jenkins-Traub root test:"<<endl;
    j_roots = t1.jt_roots();
    for (int i=0; i<5; i++) cout<<j_roots[i];
    cout<<endl;
}