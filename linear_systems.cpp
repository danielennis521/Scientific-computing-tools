#include<cmath>
#include<iostream>
#include<vector>

class lin_sys{

    public:
        lin_sys(std::vector<std::vector<float>> A, std::vector<float> b, int system_type = 0){
            dim = b.size();
            for (int i=0; i<dim; i++){
                if(A[i].size() != dim){
                    throw std::invalid_argument("system dimensions are not square");
                };
            };
            this->A = A;
            this->b = b;
            this->system_type = system_type;
            for(int i=0; i<dim; i++){
                order.push_back(i);
                solution.push_back(0.0);
            };
            //solve();
        };

        // calls the appropriate solution method for current system
        void solve(bool perturbation = false){

            general_solve();
        };

        // method for general changes to system, new solution automatically generated 
        void update(){

        };


        // for updates by a rank one matrix, in this scenario the new solution can be generated faster
        // than can be done for a general update 
        void perturbe(std::vector<float> u, std::vector<float> v){
            //u is the direction of the perturbation and v is the scaling factors

            for(int i=0; i<dim; i++){
                for(int j=0; j<dim; j++){
                    A[i][j] += u[j]*v[i];
                };
            };
            solve(true);
        };


        void display(){
            for (int i=0; i<dim; i++){
                for(int j=0; j<dim; j++){
                    std::cout<<A[j][i]<<' ';
                };
                std::cout<<' ';
                if(i == dim/2){
                    std::cout<<'=';
                }else{
                    std::cout<<' ';
                };
                std::cout<<' ';
                std::cout<<b[i]<<'\n';

            };

            for(int i=0; i<dim; i++){
                std::cout<<solution[i]<<' ';
            };
            std::cout<<'\n';

                        for(int i=0; i<dim; i++){
                std::cout<<order[i]<<' ';
            };
            std::cout<<'\n';
        };

    private: 

        void general_solve(){
            int t;
            int max;
            float sum; 

            // perform decomposition to obtain a triangular matrix
            for(int i=0; i<dim; i++){
                order[i] =  i;
            };

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

                // apply row operations                
                for(int j=i+1; j<dim; j++){

                    for(int k=i+1; k<dim; k++){
                        A[k][order[j]] -=  A[i][order[j]]*A[k][order[i]]/A[i][order[i]];
                    };
                    b[order[j]] -=  A[i][order[j]]*b[order[i]]/A[i][order[i]];
                };
            };

            // perform the back substitution
            for(int i=dim-1; i>=0; i--){
                solution[i] = b[order[i]]/A[i][order[i]];
                for(int j=0; j<i; j++){
                    b[order[j]] -= A[i][order[j]]*solution[i];

                };
            };
        };


        void banded_solve(){

        };


        void triangular_solve(){

        };


        // update the constant vector, the right hand side of Ax = b
        void update_constants(){

        };

        // update the coefficient matrix, the left hand side of Ax = b
        void update_coefficients(){
            
        };


        // specifies if the system is a general one or one of a special for for which
        // a specialized solution mathod may be used
        int system_type;
        int dim;

        std::vector<std::vector<float>> A;
        std::vector<float> b;
        std::vector<float> solution;

        std::vector<int> order;

};


int main(){

    std::vector<std::vector<float>> A = {
        {1.0, 2.0, 3.0},
        {2.0, 1.0, 2.0},
        {3.0, 2.0, 1.0}
        };

    std::vector<float> b = {1.0, 1.0, 1.0};


    lin_sys test(A, b);
    test.display();
    test.solve();
    test.display();
}