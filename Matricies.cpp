#include<cmath>
#include<iostream>
#include<vector>

// class defining matricies and useful methods
class matrix{

    public:

        matrix(std::vector<std::vector<double>> M, bool banded=false){
            this->banded = banded;
            square = true;
            dim = M[0].size();

            for(int i=1; i<M.size(); i++){
                if(dim != M[i].size()){
                    throw std::invalid_argument("Not a valid matrix");
                };
            };

            for(int i=1; i<M.size(); i++){
                order.push_back(i);
            };
            A = M;
        };

        // allow indexing
        std::vector<double>& operator[](int i){
            return A[i];
        };

        // return the dimension of the matrix
        int get_dim(){
            return dim;
        };

        // display the current values in the matrix
        void disp(){
            for(int i=0; i<dim; i++){
                for(int j=0; j<dim; j++){
                    std::cout<<A[j][order[i]]<<' ';
                };
                std::cout<<'\n';
            };
        };

        // will find either the exact solution or the least squares solution for given vector 
        std::vector<double> solve(std::vector<double> b){
            if(square){ return lin_sys(b); }
            else{ return l_sqrs(b); };
        };


    private:
        // needs to be reworked to store values for solving in the triangular matrix 
        std::vector<double> lin_sys(std::vector<double> b){
            //error handling
            if(dim != b.size()){
                throw std::invalid_argument("matrix and vector dimensions do not match");
            };

            std::vector<double> solution;
            for(int i=0; i<dim; i++){
                solution.push_back(0.0);
            };

            // obtain decomposition if needed
            if(banded){
                if(!decomp_current){ band_lu_decomp(); };
            }
            else{
                if(!decomp_current){ reg_lu_decomp(); };
            };

            // perform the forewards and backwards substitution to solve
            for(int i=0; i<dim; i++){
                solution[order[i]] = b[order[i]];
                for(int j=i+1; j<dim; j++){
                    b[order[j]] -= A[i][order[j]]*solution[order[i]];
                };
            };

            for(int i=dim-1; i>=0; i--){
                b[order[i]] = solution[order[i]]/A[i][order[i]];
                for(int j=0; j<i; j++){
                    solution[order[j]] -= A[i][order[j]]*b[order[i]];
                };
            };

            return b;
        };


        void reg_lu_decomp(){
            // needed variables for computation
            int t;
            int max; 
            std::vector<double> solution;

            for(int i=0; i<dim; i++){
                order[i] = i;
            };

            // perform decomposition to obtain a triangular matrix
            for(int i=0; i<dim; i++){
                // identify the largest element in the column
                max = i;
                for(int j=i; j<dim; j++){
                    if (std::abs(A[i][order[j]]) > std::abs(A[i][order[max]])){
                        max = j; 
                    };
                };

                // adjust the ordering
                t = order[max];
                order[max] = order[i];
                order[i] = t;

                // defines the lower matrix 
                for(int j=i+1; j<dim; j++){
                    A[i][order[j]] /= A[i][order[i]];
                };

                // apply row operations                
                for(int j=i+1; j<dim; j++){
                    for(int k=i+1; k<dim; k++){
                        A[k][order[j]] -=  A[i][order[j]]*A[k][order[i]];
                    };
                };
            };
        };


        void band_lu_decomp(){

        };


        // least squares method 
        std::vector<double> l_sqrs(std::vector<double> b){
            //error handling
            if(dim != b.size()){
                throw std::invalid_argument("matrix and vector dimensions do not match");
            };
            std::vector<double> solution;

            for(int i=0; i<dim; i++){
                solution.push_back(0.0);
            }; 
            return solution;
        };

        std::vector<std::vector<double>> A;  // the actual entries of the matrix
        std::vector<std::vector<double>> LU;  // the decomposition that will be used for linear systems
        std::vector<int> order;  // 
        int dim;  // size of the matrix
        bool square;  // true if the matrix is in fact square
        bool banded;  // true if the given matrix has multidiagonal form
        bool decomp_current; // true if no changes to matrix have been made since decomposition was last computed
        
};



int main(){

    std::vector<std::vector<double>> A = {
        {1.0, 2.0, 3.0},
        {2.0, 1.0, 2.0},
        {3.0, 2.0, 1.0},
        };

    matrix M(A);

    std::vector<double> b = {1.0, 1.0, 1.0};

    std::vector<double> sol = M.solve(b);
    M.disp();
    for(int i=0; i<M.get_dim(); i++){
        std::cout<<sol[i];
    };
}