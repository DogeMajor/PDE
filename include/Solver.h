#ifndef SOLVER_H
#define SOLVER_H
#include <iostream>
#include "PDE.h"
#include "mesh.h"
#include <math.h>

#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"


using namespace std;
using namespace Eigen;

//typedef double(*Function)(VectorXd x);
typedef Eigen::Triplet<double> Tri;


template <int Dim, typename T>
class Solver {

public:
	Solver();
	Solver(PDE<Dim, T> p, Mesh<Dim, Dim + 1, T> *m);
	~Solver() {}
	void set_pde(PDE<Dim, T> p) { pde = p; }
	void set_mesh(Mesh<Dim, Dim + 1, T> *m) { mesh = m; }
	//void build_mesh();
	//void refine();
	MatrixXd get_stiffness_matrix(int unique_nodes) const;
	//SparseMatrix<double> get_stiffness_matrix() const;
	//VectorXd solve();
	void show() const;

private:
	PDE<Dim, T> pde;
	Mesh<Dim, Dim + 1, T>* mesh;

};

template <int Dim, typename T>
Solver<Dim, T>::Solver() {
}

template <int Dim, typename T>
Solver<Dim, T>::Solver(PDE<Dim, T> p, Mesh<Dim, Dim + 1, T> *m) {
	pde = p;
	mesh = m;
}

template <int Dim, typename T>
MatrixXd Solver<Dim, T>::get_stiffness_matrix(int unique_nodes) const {
	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	MatrixXd stiffness = MatrixXd::Zero(unique_nodes, unique_nodes);
	int I, J;
	SimplexFunction<T> fn_i, fn_j;
	while (iter != nullptr) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			fn_i = iter->data.get_function(i);
			for (int j = i + 1; j < Dim + 1; j++) {
				fn_j = iter->data.get_function(j);
				J = iter->data[j].get_index();
				cout << pde.A(iter->data, fn_i, fn_j) << endl;
				stiffness(I, J) += pde.A(iter->data, fn_i, fn_j);
			}
		}
		iter = iter->next;
	}
	return stiffness;
}

/*template <int Dim, typename T>
SparseMatrix<double> Solver<Dim, T>::get_stiffness_matrix() const {
	MeshNode <Element<Dim, Dim+1, T> >* iter = mesh.get_top_mesh_node();
	std::vector<Tri> tripletList;
	tripletList.reserve((2 * dimensions.rows() + 1)*all_dims);
	while (iter->next != nullptr) {

		
		tripletList.push_back(Tri(n, n, 2 * h_factor));
	}
	SparseMatrix<double> stiffness(all_dims, all_dims);
	stiffness.setFromTriplets(tripletList.begin(), tripletList.end());
	return stiffness;
}*/

template <int Dim, typename T>
void Solver<Dim, T>::show() const {
	mesh->show();
}

#endif