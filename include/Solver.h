#ifndef SOLVER_H
#define SOLVER_H
#include <iostream>
#include "PDE.h"
#include "Mesh.h"
#include <math.h>

#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"


using namespace std;
using namespace Eigen;

//typedef double(*Function)(VectorXd x);
//typedef Eigen::Triplet<double> Tri;


template <int Dim, typename T>
class Solver {

public:
	Solver();
	Solver(PDE<Dim, T> p, Mesh<Dim, Dim + 1, T> *m, BoundaryConditions<T> b);
	~Solver() {}
	void refine();
	map<array<int, 2>, double> get_sparse_stiffness_map() const;
	SparseMatrix<double> get_sparse_stiffness_matrix(int n) const;
	MatrixXd get_stiffness_matrix(int n) const;
	SparseMatrix<double> get_sparse_inner_stiffness_matrix(map<array<int, 2>, double> stiffness_map);
	SparseMatrix<double> get_sparse_boundary_matrix(map<array<int, 2>, double> stiffness_map);
	VectorXd get_f_vec(int n) const;
	VectorXd get_boundary_coeffs();
	VectorXd conj_gradient_solve(const SparseMatrix<double> &A, VectorXd b);

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

template<int Dim, typename T>
void Solver<Dim, T>::refine() {
	mesh->refine(); mesh->reset_indices(boundaries);
}


template <int Dim, typename T>
map<array<int, 2>, double> Solver<Dim, T>::get_sparse_stiffness_map() const{
	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	map<array<int, 2>, double> stiffness_map;
	int I, J;
	SimplexFunction<T> fn_i, fn_j;
	while (iter != nullptr) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			fn_i = iter->data.get_function(i);
			for (int j = i; j < Dim + 1; j++) {
				fn_j = iter->data.get_function(j);
				J = iter->data[j].get_index();
				stiffness_map[{I, J}] += pde.A(iter->data, fn_i, fn_j);
				if (I != J) {
					stiffness_map[{J, I}] += pde.A(iter->data, fn_i, fn_j);
				}
			}
		}
		iter = iter->next;
	}
	return stiffness_map;
}

template <int Dim, typename T>
SparseMatrix<double> Solver<Dim, T>::get_sparse_stiffness_matrix(int n) const{
	
	map<array<int, 2>, double> stiffness_map = get_sparse_stiffness_map();
	SparseMatrix<double> stiffness(n + 1, n + 1);
	stiffness.reserve(VectorXi::Constant(n+1, factorial(n)));
	typedef map< array<int, 2>, double>::iterator Map_iter;
	for (Map_iter map_iter = stiffness_map.begin(); map_iter != stiffness_map.end(); map_iter++) {
		stiffness.coeffRef(map_iter->first[0], map_iter->first[1]) = map_iter->second;
	}
	//stiffness.makeCompressed();
	return stiffness;
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
				if (I != J) {
					stiffness(J, I) += pde.A(iter->data, fn_i, fn_j);
				}
			}
		}
		iter = iter->next;
	}
	return stiffness;
}

template <int Dim, typename T>
SparseMatrix<double> Solver<Dim, T>::get_sparse_inner_stiffness_matrix(map<array<int,2>, double> stiffness_map) {
	int sz = mesh->get_max_inner_index() + 1;
	SparseMatrix<double> inner_stiffness(sz, sz);
	inner_stiffness.reserve(VectorXi::Constant(sz, factorial(sz)));
	typedef map< array<int, 2>, double>::iterator Map_iter;
	for (Map_iter map_iter = stiffness_map.begin(); map_iter != stiffness_map.end(); map_iter++) {
		if ((map_iter->first[0] < sz) && (map_iter->first[1] < sz)) {
			inner_stiffness.coeffRef(map_iter->first[0], map_iter->first[1]) = map_iter->second;
		}
	}
	return inner_stiffness;
}

template <int Dim, typename T>
SparseMatrix<double> Solver<Dim, T>::get_sparse_boundary_matrix(map<array<int, 2>, double> stiffness_map) {
	int max_inner = mesh->get_max_inner_index() + 1;
	int max_outer = mesh->get_max_outer_index() + 1;
	SparseMatrix<double> bound_mat(max_inner, max_outer - max_inner);
	bound_mat.reserve(max_inner* factorial(Dim+1));
	typedef map< array<int, 2>, double>::iterator Map_iter;
	for (Map_iter map_iter = stiffness_map.begin(); map_iter != stiffness_map.end(); map_iter++) {
		if ((map_iter->first[1] >= max_inner) && (map_iter->first[0] < max_inner)) {
			bound_mat.coeffRef(map_iter->first[0], map_iter->first[1] - max_inner) = map_iter->second;
		}
	}
	return bound_mat;
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
			//if (I <= max_index) {f_vec(I) += pde.f(iter->data, fn_i);}
			if (I <= max_index) { f_vec(I) = f_vec(I) + pde.f_monte_carlo(iter->data, fn_i, i, 50); }
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

template <int Dim, typename T>//A.x = b, solve x
VectorXd Solver<Dim, T>::conj_gradient_solve(const SparseMatrix<double> &A, VectorXd b) {
	VectorXd solution(A.cols());
	ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;
	cg.compute(A);
	solution = cg.solve(b);
	b = A * solution;
	return cg.solve(b);
}

template <int Dim, typename T>
VectorXd Solver<Dim, T>::solve() {
	int inner_size = mesh->get_max_inner_index() + 1;
	int outer_size = mesh->get_max_outer_index() + 1;
	map<array<int, 2>, double> stiffness_map = get_sparse_stiffness_map();
	SparseMatrix<double> stiffness = get_sparse_inner_stiffness_matrix(stiffness_map);
	VectorXd f_vec = get_f_vec(inner_size - 1);

	VectorXd bound_coeffs = get_boundary_coeffs();
	SparseMatrix<double> bound_mat = get_sparse_boundary_matrix(stiffness_map);
	VectorXd b_vec(inner_size);
	b_vec = bound_mat * bound_coeffs;

	VectorXd vec = f_vec - b_vec;
	VectorXd inner_coeffs(inner_size);
	inner_coeffs = conj_gradient_solve(stiffness, vec);

	VectorXd coeffs(outer_size);
	coeffs.head(inner_size) = inner_coeffs;
	coeffs.tail(outer_size - inner_size) = bound_coeffs;
	return coeffs;
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
