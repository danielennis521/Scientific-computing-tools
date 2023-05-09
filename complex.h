#pragma once
#include<cmath>
#include<iostream>
#include<vector>


class complex
{
public:
    double r, i;

    complex();

    complex(double real, double imag);

    double norm();

    complex conj();

    complex pow(int n);

    complex operator+(complex const& c);

    complex operator+(double x);

    complex operator-(complex const& c);

    complex operator-(double x);

    complex operator*(complex const& c);

    complex operator*(double a);

    complex operator/(complex c);

    complex operator/(double x);

    friend std::ostream& operator<<(std::ostream &out, const complex &c);

};