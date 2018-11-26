#ifndef SOLVER_H
#define SOLVER_H
#include <iostream>
#include "PDE.h"
#include "Mesh.h"
#include "MeshFiller.h"
#include <math.h>
#include <atomic>         // std::atomic, std::atomic_flag, ATOMIC_FLAG_INIT
#include <thread>         // std::thread, std::this_thread::yield
#include <vector>
#include <future>


#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"


using namespace std;
using namespace Eigen;


template <int Dim, typename T>
class Solver {

public:
	Solver();
	Solver(PDE<Dim, T> p, Mesh<Dim, Dim + 1, T> *m, BoundaryConditions<T> b);
	~Solver() {}
	void fill_mesh_covering_box(T mid_point, VectorXd lenghts);
	void refine();
	pair<SparseMatrix<double>, VectorXd> get_sparse_stiffness_matrix_and_f_vec(int sz);
	SparseMatrix<double> get_sparse_stiffness_matrix(int sz) const;
	VectorXd get_f_vec(int n);
	VectorXd get_f_vec_async(int parts);
	VectorXd get_f_vec_part(int start_ind, int stop_ind);
	VectorXd get_boundary_coeffs();
	VectorXd conj_gradient_solve(const SparseMatrix<double> &A, VectorXd b);

	VectorXd solve(int parts = 4);//4 cores!!
	MatrixXd get_solution_values(VectorXd solution);//In the same order as indexing of nodes
	MatrixXd get_grid_values();
	void save_grid(string file_name);
	Mesh<Dim, Dim + 1, T> get_mesh() { return *mesh; }
	void show() const;

private:
	PDE<Dim, T> pde;
	Mesh<Dim, Dim + 1, T>* mesh;
	BoundaryConditions<T> boundaries;
	MeshFiller<Dim, Dim + 1, T> mesh_filler;

};

template <int Dim, typename T>
Solver<Dim, T>::Solver() {
}

template <int Dim, typename T>
Solver<Dim, T>::Solver(PDE<Dim, T> p, Mesh<Dim, Dim + 1, T> *m, BoundaryConditions<T> b) {
	pde = p;
	mesh = m;
	boundaries = b;
	mesh->set_element_divider(b);
	mesh->reset_indices(boundaries);
}

template<int Dim, typename T>
void Solver<Dim, T>::refine() {
	mesh->refine();
	mesh->reset_indices(boundaries);
}

template<int Dim, typename T>
void Solver<Dim, T>::fill_mesh_covering_box(T mid_point, VectorXd lengths) {
	vector <Vertex<Dim, T> *> box_vertices = mesh_filler.build_box_vertices(mid_point, lengths);
	mesh_filler.build_mesh(box_vertices, mesh, boundaries);
}

template <int Dim, typename T>
pair<SparseMatrix<double>, VectorXd> Solver<Dim, T>::get_sparse_stiffness_matrix_and_f_vec(int sz) {

	SparseMatrix<double> stiffness(sz, sz);
	int vec_sz = mesh->get_max_inner_index() + 1;
	VectorXd f_vec = VectorXd::Zero(vec_sz);
	stiffness.reserve(sz*factorial(Dim + 1));

	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	int max_index = mesh->get_max_inner_index();
	int I, J;
	SimplexFunction<T> fn_i, fn_j;
	while (iter != nullptr) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			fn_i = iter->data.get_function(i);
			if (I <= max_index) { f_vec(I) = f_vec(I) + pde.f_monte_carlo(iter->data, fn_i, i, 50); }
			for (int j = i; j < Dim + 1; j++) {
				fn_j = iter->data.get_function(j);
				J = iter->data[j].get_index();
				stiffness.coeffRef(I, J) = stiffness.coeffRef(I, J) + pde.A(iter->data, fn_i, fn_j);
				if (I != J) {
					stiffness.coeffRef(J, I) = stiffness.coeffRef(J, I) + pde.A(iter->data, fn_i, fn_j);
				}
			}
		}
		iter = iter->next;
	}
	stiffness.makeCompressed();
	return pair<SparseMatrix<double>, VectorXd>(stiffness, f_vec);
}


template <int Dim, typename T>
SparseMatrix<double> Solver<Dim, T>::get_sparse_stiffness_matrix(int sz) const{
	
	SparseMatrix<double> stiffness(sz, sz);
	stiffness.reserve(sz*factorial(Dim + 1));

	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	int I, J;
	SimplexFunction<T> fn_i, fn_j;
	while (iter != nullptr) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			fn_i = iter->data.get_function(i);
			for (int j = i; j < Dim + 1; j++) {
				fn_j = iter->data.get_function(j);
				J = iter->data[j].get_index();
				stiffness.coeffRef(I, J) = stiffness.coeffRef(I, J) + pde.A(iter->data, fn_i, fn_j);
				if (I != J) {
					stiffness.coeffRef(J, I) = stiffness.coeffRef(J, I) + pde.A(iter->data, fn_i, fn_j);
				}
			}
		}
		iter = iter->next;
	}
	stiffness.makeCompressed();
	return stiffness;
}


template <int Dim, typename T>
VectorXd Solver<Dim, T>::get_f_vec(int n) {
	MeshNode<Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
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

template <int Dim, typename T>//To test if parallel computing can be used
VectorXd Solver<Dim, T>::get_f_vec_async(int parts) {
	int max_nodes = mesh->how_many_nodes() + 1;
	int n_k = max_nodes / parts;
	int max_Index = mesh->get_max_inner_index();
	VectorXd f_vec = VectorXd::Zero(max_Index + 1);
	vector<std::future<VectorXd> > futures;
	
	int start_ind = 0;
	int stop_ind = 0;
	for (int k = 0; k < parts + 1; k++) {
		start_ind = k * n_k;
		stop_ind = (k + 1)*n_k;
		futures.push_back(std::async(&Solver::get_f_vec_part, this, start_ind, stop_ind));
	}

	for (int k = 0; k < parts + 1; k++) {
		f_vec = f_vec + futures[k].get();
	}

	return f_vec;
}

template <int Dim, typename T>
VectorXd Solver<Dim, T>::get_f_vec_part(int start_ind, int stop_ind) {
	if (stop_ind > mesh->how_many_nodes()) { stop_ind = mesh->how_many_nodes(); }
	int max_Index = mesh->get_max_inner_index();
	VectorXd f_vec = VectorXd::Zero(max_Index+1);
	int I;
	SimplexFunction<T> fn_i;
	int index = start_ind;
	MeshNode<Element<Dim, Dim + 1, T> >* iter = mesh->get_node(index);
	while (index < stop_ind) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			fn_i = iter->data.get_function(i);
				//if (I <= max_index) {f_vec(I) += pde.f(iter->data, fn_i);}
			if (I <= max_Index) { f_vec(I) = f_vec(I) + pde.f_monte_carlo(iter->data, fn_i, i, 50); }
		}
		index++;
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
				if (I > min_I - 1) {
					loc = iter->data[i].get_location();
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
VectorXd Solver<Dim, T>::solve(int parts) {
	Timer timer = Timer();
	int inner_size = mesh->get_max_inner_index() + 1;
	int outer_size = mesh->get_max_outer_index() + 1;

	SparseMatrix<double> total_stiffness = get_sparse_stiffness_matrix(outer_size);
	//VectorXd f_vec = get_f_vec(inner_size - 1);
	VectorXd f_vec = get_f_vec_async(parts);
	SparseMatrix<double> stiffness = total_stiffness.block(0, 0, inner_size, inner_size);
	
	SparseMatrix<double> bound_mat = total_stiffness.block(0, inner_size, inner_size, outer_size - inner_size);
	VectorXd bound_coeffs = get_boundary_coeffs();
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
	SolverDAO<Dim, Dim + 1, T> solver_dao;
	return solver_dao.get_solution_values(mesh, solution);
}

template <int Dim, typename T>
MatrixXd Solver<Dim, T>::get_grid_values() {
	MeshDAO<Dim, Dim + 1, T> mesh_dao;
	return mesh_dao.get_grid_values(mesh);
}

template <int Dim, typename T>
void Solver<Dim, T>::save_grid(string file_name) {
	MeshDAO<Dim, Dim + 1, T> mesh_dao;
	mesh_dao.save_grid(file_name, mesh);
}

template <int Dim, typename T>
void Solver<Dim, T>::show() const {
	mesh->show();
}

#endif
