#include<cmath>
#include<iostream>
#include<vector>
#include<complex>
#include "Matricies.h"

using namespace std;


matrix::matrix(vector<vector<complex<double>>> M){
    if(M[0].size() == M.size()){ square = true;
    }else{ square = false;};
    row_num = M[0].size();
    col_num = M.size();
    A = M;
    qr_current = false;
    lu_current = false;

    if (row_num != col_num)square = false;

    for(int i=1; i<M.size(); i++){
        if(row_num != M[i].size()) throw invalid_argument("Not a valid matrix");
    };

    for(int i=1; i<M.size(); i++) order.push_back(i);
};


matrix::matrix(vector<vector<double>> M){
    if(M[0].size() == M.size()){ square = true;
    }else{ square = false;};
    row_num = M[0].size();
    col_num = M.size();
    complex<double> c;
    qr_current = false;
    lu_current = false;

    if (row_num != col_num)square = false;

    for(int i=1; i<M.size(); i++){
        if(row_num != M[i].size()) throw invalid_argument("Not a valid matrix");
    };

    for(int i=1; i<M.size(); i++) order.push_back(i);

    for(int i=0; i<M.size(); i++){
        A.push_back({});
        for(int j=0; j<M[0].size(); j++){
            c = {M[i][j], 0.0};
            A[i].push_back(c);
        };
    };
};


vector<complex<double>>& matrix::operator[](int i){
    return A[i];
};


int matrix::get_dim(){
    return row_num;
};


void matrix::disp(){
    for(int i=0; i<col_num; i++){
        for(int j=0; j<row_num; j++){
            cout<<A[i][j]<<' ';
        };
        cout<<'\n';
    };
};


vector<complex<double>> matrix::solve(vector<complex<double>> b){
    if(square){ return lin_sys(b); }
    else{ return l_sqrs(b); };
};

vector<complex<double>> matrix::solve(vector<double> b){
    vector<complex<double>> v;
    complex<double> c;
    for(int i=0; i<b.size(); i++){
        c = {b[i], 0.0};
        v.push_back(c);
    };

    return solve(v);
};


vector<complex<double>> matrix::map(vector<complex<double>> v){
    vector<complex<double>> solution;
    for (int i=0; i< col_num; i++) solution.push_back(0.0);

    for(int i=0; i<col_num; i++){
        for(int j=0; j<row_num; j++) solution[i] += A[i][j] * v[j];
    };

    return solution;
};


void matrix::transform(){
    vector<vector<complex<double>>> M;
    int t;
    for(int i=0; i<row_num; i++){
        M.push_back({});
        for(int j=0; j<col_num; j++) M[i].push_back(A[j][i]);
    };

    t = col_num;
    col_num = row_num;
    row_num = t;
    A = M;
};


vector<complex<double>> matrix::max_eigen(){
    vector<complex<double>> p, x, y;
    double m=0.0;
    double c=1.0;
    for(int i=0; i<col_num; i++) x.push_back(1.0);
    for(int i=0; i<col_num; i++) p.push_back(1.0);

    while(c>1e-10){
        p=x;
        y = map(x);
        x = y;
        for(int i=0; i<col_num; i++) if(norm(y[i]) > m) m = norm(y[i]);
        for(int i=0; i<col_num; i++) x[i] = x[i]/m;
        
        c = 0.0;
        for(int i=0; i<col_num; i++) c+= norm(x[i] - p[i]);
    };

    return x;
};


vector<complex<double>> matrix::all_eigen_vals(){ // note this method only works for a square matrix

    if(col_num != row_num){
        //throw an error
    };
    bool c=true;
    complex<double> sig={0.0, 0.0};
    vector<complex<double>> x;
    vector<vector<complex<double>>> S=A;

    for(int i=0; i< col_num; i++) x.push_back({0.0, 0.0});

    while(c){
        sig=A[col_num-1][row_num-1];
        
        for(int i=0; i<col_num; i++) A[i][i] -= sig;
        qr_decomp();
        // do the matrix multiplication
        for(int i=0; i<row_num; i++){
            for(int j=0; j<row_num; j++){
                A[j][i] = 0.0;
                for(int k=0; k<row_num; k++){
                    A[j][i] += R[k][i]*Q[j][k];
                };
            };
        };

        c = false;
        for(int i=0; i<row_num-1; i++){
            A[i][i] += sig;
            if(norm(A[i][row_num-1]) > 1e-8) c=true;
        };
        A[row_num-1][row_num-1] += sig;

        cout<<endl;
        disp();
        cout<<endl;
    };

    for(int i=0; i< col_num; i++) x[i] = A[i][i];
    A=S;
    qr_current = false;

    return x;
};


vector<vector<complex<double>>> matrix::all_eigen_vec(){
    vector<complex<double>> b,  x = all_eigen_vals();
    for(int i=0; i<row_num; i++) b.push_back({0.0, 0.0});
    vector<vector<complex<double>>> R, S=A;

    for(int i=0; i<x.size(); i++){
        A=S;
        for(int j=0; j<row_num; j++) A[j][j] -= x[i];
        R.push_back(solve(b));
    };

    return R;
};


vector<complex<double>> matrix::lin_sys(vector<complex<double>> b){
    if(row_num != b.size()){
        throw invalid_argument("matrix and vector dimensions do not match");
    };

    vector<complex<double>> solution;
    for(int i=0; i<row_num; i++){
        solution.push_back(0.0);
    };

    if(!lu_current){ lu_decomp(); };

    for(int i=0; i<row_num; i++){               // solve lower triangular system (forward substitution)
        solution[order[i]] = b[order[i]];
        for(int j=i+1; j<row_num; j++){
            b[order[j]] -= LU[i][order[j]]*solution[order[i]];
        };
    };

    for(int i=row_num-1; i>=0; i--){            // solve upper triangular system (backward substitution)
        b[order[i]] = solution[order[i]]/LU[i][order[i]];
        for(int j=0; j<i; j++) solution[order[j]] -= LU[i][order[j]]*b[order[i]];
    };

    for(int i=0; i<row_num; i++) solution[i] = b[order[i]];
    return solution;
};


vector<complex<double>> matrix::l_sqrs(vector<complex<double>> b){
    vector<complex<double>> solution;
    complex<double> beta, gamma;

    if(row_num != b.size()){
        throw invalid_argument("matrix and vector dimensions do not match");
    };

    if(!qr_current) qr_decomp();

    for(int i=0; i<col_num; i++){   // transformations
        beta = {0.0, 0.0};
        solution.push_back({0.0, 0.0});
        for(int j=0; j<row_num; j++) solution[i] += Q[i][j]*b[j];
    };

    for(int i=col_num-1; i>=0; i--){   // back-substitution
        b[i] = solution[i]/R[i][i];
        for(int j=0; j<i; j++) solution[j] -= R[i][j]*b[i];
    };

    for(int i=0; i<col_num; i++) solution[i] = b[i];
    return solution;
};


void matrix::lu_decomp(){   // lu decomposition for solving linear systems
    LU = A;
    int t;
    int max; 
    vector<complex<double>> solution;

    for(int i=0; i<row_num; i++){
        order[i] = i;
    };

    for(int i=0; i<row_num; i++){ 
        max = i;  
        for(int j=i; j<row_num; j++){           // "pivot" i.e. adjust row ordering to ensure numerical stability   
            if (abs(LU[i][order[j]]) > abs(LU[i][order[max]])) max = j; 
        };
        t = order[max];
        order[max] = order[i];
        order[i] = t;

        for(int j=i+1; j<row_num; j++){         // defines the lower matrix in place
            LU[i][order[j]] /= LU[i][order[i]];
        };
        
        for(int j=i+1; j<row_num; j++){         // apply the transformation to the rest of the matrix
            for(int k=i+1; k<row_num; k++){
                LU[k][order[j]] -=  LU[i][order[j]]*LU[k][order[i]];
            };
        };
    };

    lu_current = true;
};


void matrix::qr_decomp(){
    vector<complex<double>> v;
    complex<double> c={0.0, 0.0};
    Q=A;
    R.clear();
    for(int i=0; i<col_num; i++) v.push_back(c);
    for(int i=0; i<col_num; i++) R.push_back(v);

    for(int i=0; i<col_num; i++){
        R[i][i] = {0.0, 0.0};
        for(int j=0; j<row_num; j++) R[i][i] += Q[i][j]*Q[i][j];
        R[i][i] = sqrt(R[i][i]);

        if(R[i][i] == 0.0){
            // throw an error
        };

        for(int j=0; j<row_num; j++) Q[i][j] /= R[i][i];

        for(int j=i+1; j<col_num; j++){
            R[j][i] = {0.0, 0.0};
            for(int k=0; k<row_num; k++) R[j][i] += Q[i][k]*Q[j][k];
            for(int k=0; k<row_num; k++) Q[j][k] -= Q[i][k]*R[j][i];
        };
    };

    qr_current = true;
};


matrix matrix::operator*(matrix B){

};


matrix matrix::operator*(complex<double> z){

};
