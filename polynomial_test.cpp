#include<cmath>
#include<iostream>
#include<vector>
#include "complex.h"
#include "polynomials.h"
using namespace std;


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