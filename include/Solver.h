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
	Solver(PDE<Dim, T> p, Mesh<Dim, Dim + 1, T> *m, BoundaryConditions<T> b);
	~Solver() {}
	void set_pde(PDE<Dim, T> p) { pde = p; }
	void set_mesh(Mesh<Dim, Dim + 1, T> *m) { mesh = m; mesh->reset_indices();}
	MatrixXd get_stiffness_matrix(int n) const;
	MatrixXd get_inner_stiffness_matrix(MatrixXd stiffness);
	MatrixXd get_boundary_matrix(MatrixXd stiffness);
	VectorXd get_f_vec(int n) const;
	VectorXd get_boundary_coeffs();
	void refine() { mesh->refine(); mesh->reset_indices(boundaries); }
	VectorXd solve();
	MatrixXd get_solution_values(VectorXd solution);//In the same order as indexing of nodes
	Mesh<Dim, Dim + 1, T> get_mesh() { return *mesh; }
	void show() const;

private:
	PDE<Dim, T> pde;
	Mesh<Dim, Dim + 1, T>* mesh;
	BoundaryConditions<T> boundaries;

};

template <int Dim, typename T>
Solver<Dim, T>::Solver() {
}

template <int Dim, typename T>
Solver<Dim, T>::Solver(PDE<Dim, T> p, Mesh<Dim, Dim + 1, T> *m, BoundaryConditions<T> b) {
	pde = p;
	mesh = m;
	boundaries = b;
	mesh->reset_indices(boundaries);
}

template <int Dim, typename T>
MatrixXd Solver<Dim, T>::get_stiffness_matrix(int n) const {
	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	MatrixXd stiffness = MatrixXd::Zero(n+1, n+1);
	int I, J;
	SimplexFunction<T> fn_i, fn_j;
	while (iter != nullptr) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			fn_i = iter->data.get_function(i);
			for (int j = i; j < Dim + 1; j++) {
				fn_j = iter->data.get_function(j);
				J = iter->data[j].get_index();
				stiffness(I, J) += pde.A(iter->data, fn_i, fn_j);
			}
		}
		iter = iter->next;
	}
	return stiffness;
}

template <int Dim, typename T>
MatrixXd Solver<Dim, T>::get_inner_stiffness_matrix(MatrixXd stiffness) {
	int sz = mesh->get_max_inner_index() + 1;
	return stiffness.block(0, 0, sz, sz);
}

template <int Dim, typename T>
MatrixXd Solver<Dim, T>::get_boundary_matrix(MatrixXd stiffness) {
	int rows = mesh->get_max_inner_index() + 1;
	int cols = mesh->get_max_outer_index() - mesh->get_max_inner_index();
	return stiffness.block(0, rows, rows, cols);
}

template <int Dim, typename T>
VectorXd Solver<Dim, T>::get_f_vec(int n) const {
	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	VectorXd f_vec = VectorXd::Zero(n + 1);
	int I;
	int max_index = mesh->get_max_inner_index();
	SimplexFunction<T> fn_i;
	while (iter != nullptr) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			fn_i = iter->data.get_function(i);
			if (I <= max_index) {f_vec(I) += pde.f(iter->data, fn_i);}
		}
		iter = iter->next;
	}
	return f_vec;
}

template <int Dim, typename T>
VectorXd Solver<Dim, T>::get_boundary_coeffs() {
	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	int sz = mesh->get_max_outer_index() - mesh->get_max_inner_index();
	int min_I = mesh->get_max_inner_index() + 1;
	VectorXd coeffs(sz);
	int I;
	T loc;
	while (iter != nullptr) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			loc = iter->data[i].get_location();
			if (boundaries.cond(loc)) {
				coeffs(I - min_I) = boundaries.val(loc);
			}
		}
		iter = iter->next;
	}
	return coeffs;
}

template <int Dim, typename T>
VectorXd Solver<Dim, T>::solve() {
	int sol_size = mesh->get_max_outer_index() + 1;
	int inner_nodes = mesh->get_max_inner_index() + 1;
	MatrixXd total_stiffness = get_stiffness_matrix(mesh->get_max_outer_index());
	MatrixXd stiffness = get_inner_stiffness_matrix(total_stiffness);
	VectorXd f_vec = get_f_vec(mesh->get_max_inner_index());

	MatrixXd bound_mat = get_boundary_matrix(total_stiffness);
	VectorXd bound_coeffs = get_boundary_coeffs();
	VectorXd b_vec = bound_mat * bound_coeffs;

	VectorXd inner_coeffs = stiffness.inverse()*(f_vec-b_vec);
	VectorXd sol(sol_size);
	sol.head(inner_nodes) = inner_coeffs;
	sol.tail(sol_size - inner_nodes) = bound_coeffs;
	return sol;
}

template <int Dim, typename T>
MatrixXd Solver<Dim, T>::get_solution_values(VectorXd solution) {
	MatrixXd values = MatrixXd::Zero(solution.size(), Dim + 1);
	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	int I;
	T loc;
	while (iter != nullptr) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			loc = iter->data[i].get_location();
			for (int J = 0; J < Dim; J++) {
				values(I, J) = loc[J];
			}
			values(I, Dim) = solution(I);
		}
		iter = iter->next;
	}
	return values;
}

template <int Dim, typename T>
void Solver<Dim, T>::show() const {
	mesh->show();
}

#endif
