#include "../PDE/include/poisson.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Core"

#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"

using namespace Eigen;

using namespace std;


/*Vector dyn_func(Vector x){
    Vector y = MatrixXd::Zero(x.rows(), 1);
    for(int i = 0; i < x.rows(); i++){
        y(i) = (3*x(i) + pow(x(i),2))*exp(x(i));
    }
    return y;
}*/

int main(int, char *[]){
	int m = 4;
	VectorXd u(m), b(m);
	SparseMatrix<double> B(m, m);
	// fill A and b
	B.insert(0, 0) = 1;
	B.insert(1, 1) = 1;
	B.insert(2, 2) = 1;
	B.insert(3, 3) = 2;
	b << 1, 2, 3, 4;
	ConjugateGradient<SparseMatrix<double>, Upper> cg;
	cg.compute(B);
	u = cg.solve(b);
	cout << u;
	/*
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
    cout << solution << endl;*/
    }



