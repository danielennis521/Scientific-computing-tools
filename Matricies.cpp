#include<cmath>
#include<iostream>
#include<vector>
#include<complex>
#include "Matricies.h"

using namespace std;


matrix::matrix(vector<vector<complex<double>>> M){
    square = true;
    row_dim = M[0].size();
    col_dim = M.size();
    A = M;

    if (row_dim != col_dim)square = false;

    for(int i=1; i<M.size(); i++){
        if(row_dim != M[i].size()) throw invalid_argument("Not a valid matrix");
    };

    for(int i=1; i<M.size(); i++) order.push_back(i);
};


matrix::matrix(vector<vector<double>> M){
    square = true;
    row_dim = M[0].size();
    col_dim = M.size();
    complex<double> c;

    if (row_dim != col_dim)square = false;

    for(int i=1; i<M.size(); i++){
        if(row_dim != M[i].size()) throw invalid_argument("Not a valid matrix");
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
    return row_dim;
};


void matrix::disp(){
    for(int i=0; i<row_dim; i++){
        for(int j=0; j<row_dim; j++){
            cout<<A[j][i]<<' ';
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
    for (int i=0; i< col_dim; i++) solution.push_back(0.0);

    for(int i=0; i<col_dim; i++){
        for(int j=0; j<row_dim; j++) solution[i] += A[i][j] * v[j];
    };

    return solution;
};


void matrix::transform(){
    vector<vector<complex<double>>> M;
    int t;
    for(int i=0; i<row_dim; i++){
        M.push_back({});
        for(int j=0; j<col_dim; j++) M[i].push_back(A[j][i]);
    };

    t = col_dim;
    col_dim = row_dim;
    row_dim = t;
    A = M;
};


vector<complex<double>> matrix::max_eigen(){

};


vector<complex<double>> matrix::lin_sys(vector<complex<double>> b){
    if(row_dim != b.size()){
        throw invalid_argument("matrix and vector dimensions do not match");
    };

    vector<complex<double>> solution;
    for(int i=0; i<row_dim; i++){
        solution.push_back(0.0);
    };

    if(!decomp_current){ lu_decomp(); };

    for(int i=0; i<row_dim; i++){               // solve lower triangular system (forward substitution)
        solution[order[i]] = b[order[i]];
        for(int j=i+1; j<row_dim; j++){
            b[order[j]] -= LU[i][order[j]]*solution[order[i]];
        };
    };

    for(int i=row_dim-1; i>=0; i--){            // solve upper triangular system (backward substitution)
        b[order[i]] = solution[order[i]]/LU[i][order[i]];
        for(int j=0; j<i; j++){
            solution[order[j]] -= LU[i][order[j]]*b[order[i]];
        };
    };

    for(int i=0; i<row_dim; i++) solution[i] = b[order[i]];
    return solution;
};


vector<complex<double>> matrix::l_sqrs(vector<complex<double>> b){
    vector<complex<double>> solution;
    complex<double> beta, gamma;
    for(int i=0; i<row_dim; i++) solution.push_back(0.0);

    if(row_dim != b.size()){
        throw invalid_argument("matrix and vector dimensions do not match");
    };

    if(!decomp_current) qr_decomp();

    for(int i=0; i<col_dim; i++){              // apply the stransformations to the given values
        beta = LU[i][row_dim]*LU[i][row_dim];
        for(int j=0; j<row_dim; j++) beta += LU[i][j]*LU[i][j];

        gamma = LU[i][row_dim]*b[i];
        for(int k=i; k<row_dim; k++) gamma += LU[i][k]*LU[i][k];

        b[i] -= 2.0*gamma*LU[i][row_dim]/beta;
        for(int k=i+1; k<row_dim; k++) b[k] -= 2.0*gamma*LU[i][k]/beta;
    };

    for(int i=col_dim-1; i>=0; i--){            // backsubstitution to find the solution
        solution[i] = b[i]/A[i][i];
        for(int j=0; j<i; j++){
            b[j] -= A[i][j]*solution[i];
        };
    };

    return solution;
};


void matrix::lu_decomp(){   // lu decomposition for solving linear systems
    LU = A;
    int t;
    int max; 
    vector<complex<double>> solution;

    for(int i=0; i<row_dim; i++){
        order[i] = i;
    };

    for(int i=0; i<row_dim; i++){ 
        max = i;  
        for(int j=i; j<row_dim; j++){           // "pivot" i.e. adjust row ordering to ensure numerical stability   
            if (abs(LU[i][order[j]]) > abs(LU[i][order[max]])) max = j; 
        };
        t = order[max];
        order[max] = order[i];
        order[i] = t;

        for(int j=i+1; j<row_dim; j++){         // defines the lower matrix in place
            LU[i][order[j]] /= LU[i][order[i]];
        };
        
        for(int j=i+1; j<row_dim; j++){         // apply the transformation to the rest of the matrix
            for(int k=i+1; k<row_dim; k++){
                LU[k][order[j]] -=  LU[i][order[j]]*LU[k][order[i]];
            };
        };
    };

    decomp_current = true;
};

        
vector<complex<double>> matrix::qr_decomp(){     // qr decomposition for lest squares method
    LU = A;
    complex<double> alpha, beta, gamma;

    for(int i=0; i<row_dim; i++) LU[i].push_back(0.0);
    
    for(int i=0; i<col_dim; i++){
        alpha = 0.0;
        for(int j=0; j<row_dim; j++) alpha += LU[i][j]*LU[i][j];
        alpha = sqrt(alpha);

        LU[i][row_dim] = LU[i][i] - alpha;

        beta = LU[i][row_dim]*LU[i][row_dim];
        for(int j=0; j<row_dim; j++) beta += LU[i][j]*LU[i][j];

        for(int j=i; j<col_dim; j++){       // apply to remaining submatrix
            gamma = LU[i][row_dim]*LU[j][i];
            for(int k=i; k<row_dim; k++) gamma += LU[i][k]*LU[j][k];
            
            LU[j][i] -= 2.0*gamma*LU[i][row_dim]/beta;
            for(int k=i+1; k<row_dim; k++) LU[j][k] -= 2.0*gamma*LU[j][k]/beta;
        };
    };
    decomp_current = true;
};
