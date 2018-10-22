#ifndef SOLVER_H
#define SOLVER_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Core"
#include "node.h"
#include "element.h"
#include "mesh.h"
#include <math.h>

#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"
#include <iostream>
#include "Function.h"
#include <vector>

using namespace std;
using namespace Eigen;

typedef double(*Function)(VectorXd x);
typedef Eigen::Triplet<double> Tri;


template <int Dim, typename T>
class Solver {

public:
	Solver();
	Solver(PDE<Element<Dim, Dim + 1, T> > &p, Mesh<Dim, Dim + 1, T> &m) :pde(p) mesh(m) {}
	~Solver() {}
	void build_mesh();
	void refine();
	SparseMatrix<double> get_stiffness_matrix() const;
	VectorXd solve();
	void show() const;

private:
	PDE<Element<Dim, Dim + 1, T> > pde;
	Mesh<Dim, Dim + 1, T> mesh;

};

template <int Dim, typename T>
Solver<Dim, N, T>::Solver() {
}

template <int Dim, typename T>
void Solver<Dim, N, T>::show() const {
	mesh.show();
}
#endif
