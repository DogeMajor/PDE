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
	bool is_at_boundary(VectorXd coords) { return boundaries.cond(coords); }
	double value_at_boundary(VectorXd coords) { return boundaries.val(coords); }
	//void build_mesh();
	//void refine();
	MatrixXd get_stiffness_matrix(int n) const;
	MatrixXd get_inner_stiffness_matrix(int n) const;
	VectorXd get_vector_part(int n) const;
	VectorXd set_boundary_vals(VectorXd &sol);
	//SparseMatrix<double> get_stiffness_matrix() const;
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
	//cout << "max_index" <<max_index << endl;
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
				//cout << pde.A(iter->data, fn_i, fn_j) << endl;
				stiffness(I, J) += pde.A(iter->data, fn_i, fn_j);
				/*if (I == J && I == 0) {
					iter->data[i].show();
					cout << "volume: " << iter->data.get_volume() << endl;
					cout << fn_i.coeff << endl;
					cout << fn_j.coeff << endl;
					cout << "Integration result" << endl;
					cout << pde.A(iter->data, fn_i, fn_j) << endl;
					cout << "Diag. el. for index " << I << ": " << stiffness(I, J) << endl; 
				}*/
			}
		}
		iter = iter->next;
	}
	return stiffness;
}

template <int Dim, typename T>
VectorXd Solver<Dim, T>::get_vector_part(int n) const {
	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	//cout << "max_index" <<max_index << endl;
	VectorXd f_vec = VectorXd::Zero(n + 1);
	int I;
	SimplexFunction<T> fn_i;
	while (iter != nullptr) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			fn_i = iter->data.get_function(i);
			f_vec(I) += pde.f(iter->data, fn_i);
		}
		iter = iter->next;
	}
	return f_vec;
}

template <int Dim, typename T>
VectorXd Solver<Dim, T>::set_boundary_vals(VectorXd &sol) {
	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	int I;
	T loc;
	while (iter != nullptr) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			loc = iter->data[i].get_location();
			if (boundaries.cond(loc)) {
				sol(I) = boundaries.val(loc);
			}
		}
		iter = iter->next;
	}
	return sol;
}

template <int Dim, typename T>
MatrixXd Solver<Dim, T>::get_inner_stiffness_matrix(int n) const {
	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	//cout << "max_index" <<max_index << endl;
	MatrixXd stiffness = MatrixXd::Zero(n + 1, n + 1);
	int I, J;
	T loc_i, loc_j;
	SimplexFunction<T> fn_i, fn_j;
	int max_index = mesh->get_max_inner_index();
	while (iter != nullptr) {
		for (int i = 0; i < Dim + 1; i++) {
			I = iter->data[i].get_index();
			fn_i = iter->data.get_function(i);
			loc_i = iter->data[i].get_location();
			for (int j = i; j < Dim + 1; j++) {
				
				fn_j = iter->data.get_function(j);
				J = iter->data[j].get_index();
				loc_j = iter->data[j].get_location();
				//cout << pde.A(iter->data, fn_i, fn_j) << endl;
				if ((I <= max_index) && (J <= max_index)) {
					cout << "index added to inner_stiffness matrix at " << I << "," << J << endl;
					stiffness(I, J) += pde.A(iter->data, fn_i, fn_j);
				}
				
				/*if (I == J && I == 0) {
					iter->data[i].show();
					cout << "volume: " << iter->data.get_volume() << endl;
					cout << fn_i.coeff << endl;
					cout << fn_j.coeff << endl;
					cout << "Integration result" << endl;
					cout << pde.A(iter->data, fn_i, fn_j) << endl;
					cout << "Diag. el. for index " << I << ": " << stiffness(I, J) << endl;
				}*/
			}
		}
		iter = iter->next;
	}
	return stiffness;
}

template <int Dim, typename T>
VectorXd Solver<Dim, T>::solve() {
	MatrixXd stiffness = get_stiffness_matrix(mesh->get_max_outer_index());
	VectorXd f_vec = get_vector_part(mesh->get_max_outer_index());
	cout << "Showing stiffness matrix!" << endl;
	cout << stiffness << endl;
	cout << f_vec << endl;
	return stiffness.inverse()*f_vec;
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
