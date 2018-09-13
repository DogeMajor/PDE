#ifndef POISSON_H
#define POISSON_H
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"
#include <iostream>

using namespace Eigen;

typedef Eigen::Matrix<double, Dynamic, 1> Vector;
//typedef Eigen::Matrix<int, Dynamic, 1> IntVector; //Jag e så finlandsvensk...
typedef Vector (* Function)(Vector x);
typedef unsigned int uint;
typedef Eigen::Triplet<double> T; //For filling sparse matrices

class Poisson
{
    public:
        Poisson(double m_error, VectorXi dims, MatrixXd dom, Function fn);
        MatrixXd get_diff_matrix() const;
        Vector derivative(Vector u) const;
        Vector solve();
        void show() const;
    protected:
        double max_error;
        VectorXi dimensions;
        MatrixXd domain;
        SparseMatrix<double> A; //is equal to Matrix<double, Dynamic, Dynamic>
        //MatrixXf A;
        VectorXd h;
        Function func;

};

#endif
