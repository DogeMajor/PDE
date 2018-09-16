#ifndef POISSON_H
#define POISSON_H
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"
#include <iostream>
#include <math.h>

using namespace Eigen;

typedef Eigen::Matrix<double, Dynamic, 1> Vector;
typedef Vector (* VectorFunction)(Vector x);
typedef double (* Function)(Vector x);
typedef unsigned int uint;
typedef Eigen::Triplet<double> T; //For filling sparse matrices

class Poisson
{
    public:
        Poisson(double m_error, VectorXi dims, MatrixXd dom, VectorFunction fn, Function scalar_fn);
        void set_matrix();
        int get_shift_number(int max_dim);
        int get_under_dim(int dim_number);
        int to_index(VectorXi coords);
        VectorXi to_coords(int index);
        double eval_func(VectorXi coords);
        VectorXd vectorize_scalar_func();
        SparseMatrix<double> get_diff_matrix() const;
        Vector derivative(Vector u) const;
        Vector solve();
        void show() const;
    protected:
        double max_error;
        VectorXi dimensions;
        MatrixXd domain;
        SparseMatrix<double> A;
        VectorXd h;
        VectorFunction func;
        Function scalar_func;

};

#endif
