#include "poisson.h"
#include <math.h>

using namespace std;


Vector dyn_func(Vector x){
    Vector y = MatrixXd::Zero(x.rows(), 1);
    for(int i = 0; i < x.rows(); i++){
        y(i) = (3*x(i) + pow(x(i),2))*exp(x(i));
    }
    return y;
}
/*
int main(int, char *[]){
    VectorXi dims = MatrixXi::Zero(1, 1);
    dims(0) = 5;
    cout << dims.cols() << endl;
    double max_error = pow(-12, 10);
    MatrixXd domain = MatrixXd::Zero(2, 1);
    domain(0,0) = 0.0;
    domain(1,0) = 1.0;
    Poisson poisson = Poisson(max_error, dims, domain, dyn_func);
    poisson.show();
    VectorXd solution = poisson.solve();
    cout << "Solution: " << endl;
    cout << solution << endl;
    }*/



