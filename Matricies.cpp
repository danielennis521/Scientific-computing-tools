#include<cmath>
#include<iostream>
#include<vector>


class matrix{

    public:

        matrix(std::vector<std::vector<double>> M){
            square = true;
            row_dim = M[0].size();
            col_dim = M.size();
            A = M;

            if (row_dim != col_dim)square = false;

            for(int i=1; i<M.size(); i++){
                if(row_dim != M[i].size()) throw std::invalid_argument("Not a valid matrix");
            };

            for(int i=1; i<M.size(); i++) order.push_back(i);
        };


        std::vector<double>& operator[](int i){
            return A[i];
        };


        int get_dim(){
            return row_dim;
        };


        void disp(){
            for(int i=0; i<row_dim; i++){
                for(int j=0; j<row_dim; j++){
                    std::cout<<A[j][order[i]]<<' ';
                };
                std::cout<<'\n';
            };
        };


        std::vector<double> solve(std::vector<double> b){
            if(square){ return lin_sys(b); }
            else{ return l_sqrs(b); };
        };


        std::vector<double> transform(std::vector<double> v){
            std::vector<double> solution;
            for (int i=0; i< col_dim; i++) solution.push_back(0.0);

            for(int i=0; i<col_dim; i++){
                for(int j=0; j<row_dim; j++) solution[i] += A[i][j] * v[j];
            };

            return solution;
        };


        std::vector<double> max_eigen(){
            double lmbda = 0.0;
            double max = 1.0;
            std::vector<double> v;
            for (int i=0; i< col_dim; i++){
                v.push_back(std::rand()/RAND_MAX);
            };

            while (std::abs(lmbda - max) > 0.000001){
                lmbda = max;
                v = transform(v);
                max = std::abs(v[0]);
                for (int i=1; i< col_dim; i++){
                    if (std::abs(v[i]) > max) max = std::abs(v[i]);
                };

                for (int i=1; i< col_dim; i++) v[i] /= max;
            };

            lmbda = max;
            max = 0.0;
            for (int i=1; i< col_dim; i++) max += v[i]*v[i];
            max = std::sqrt(max);
            for (int i=1; i< col_dim; i++) v[i] /= max;
            return v;
        };


        matrix all_eigen(){
        };


    private:
        std::vector<double> lin_sys(std::vector<double> b){
            if(row_dim != b.size()){
                throw std::invalid_argument("matrix and vector dimensions do not match");
            };

            std::vector<double> solution;
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

            return b;
        };


        std::vector<double> l_sqrs(std::vector<double> b){
            std::vector<double> solution;
            double beta;
            double gamma;
            for(int i=0; i<row_dim; i++) solution.push_back(0.0);

            if(row_dim != b.size()){
                throw std::invalid_argument("matrix and vector dimensions do not match");
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


        // lu decomposition for solving linear systems
        void lu_decomp(){
            LU = A;
            int t;
            int max; 
            std::vector<double> solution;

            for(int i=0; i<row_dim; i++){
                order[i] = i;
            };

            for(int i=0; i<row_dim; i++){ 
                max = i;  
                for(int j=i; j<row_dim; j++){           // "pivot" i.e. adjust row ordering to ensure numerical stability   
                    if (std::abs(LU[i][order[j]]) > std::abs(LU[i][order[max]])) max = j; 
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


        // qr decomposition for lest squares method
        std::vector<double> qr_decomp(){
            LU = A;
            double sgn;
            double alpha;
            double beta;
            double gamma;

            for(int i=0; i<row_dim; i++) LU[i].push_back(0.0);
            
            for(int i=0; i<col_dim; i++){
                if(LU[i][i] < 0.0){sgn = -1.0;}
                else{sgn = 1.0;}
                alpha = 0.0;
                for(int j=0; j<row_dim; j++) alpha += LU[i][j]*LU[i][j];
                alpha = std::sqrt(alpha);

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



        std::vector<std::vector<double>> A;         // the actual entries of the matrix
        std::vector<std::vector<double>> LU;        // the triangular decomposition of the matrix
        std::vector<int> order;                     // track order of rows rather than actually interchange
        int row_dim;
        int col_dim;  
        bool square;  
        bool decomp_current;
        
};



int main(){

    std::vector<std::vector<double>> A = {
        {2.0, 1.0, 0.0, 0.0, 0.0},
        {1.0, 2.0, 1.0, 0.0, 0.0},
        {0.0, 1.0, 2.0, 1.0, 0.0},
        {0.0, 0.0, 1.0, 2.0, 1.0},
        {0.0, 0.0, 0.0, 1.0, 2.0}
        };

    matrix M(A);

    std::vector<double> b = {1.0, 1.0, 1.0, 1.0, 1.0};

    std::vector<double> sol = M.solve(b);
    M.disp();
    for(int i=0; i<M.get_dim(); i++){
        std::cout<<sol[i]<<' ';
    };
}
