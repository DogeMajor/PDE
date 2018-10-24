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
	void set_mesh(Mesh<Dim, Dim + 1, T> *m) { mesh = m; max_index = mesh->reset_indices();}
	//void build_mesh();
	//void refine();
	MatrixXd get_stiffness_matrix(int unique_nodes) const;
	VectorXd get_vector_part() const;
	//SparseMatrix<double> get_stiffness_matrix() const;
	void refine() { mesh->refine(); max_index = mesh->reset_indices(-1); }
	VectorXd solve();
	MatrixXd get_solution_values(VectorXd solution);//In the same order as indexing of nodes
	void show() const;

private:
	PDE<Dim, T> pde;
	Mesh<Dim, Dim + 1, T>* mesh;
	int max_index;//Hoe many unique Node<Dim,  T>s exist in the Mesh

};

template <int Dim, typename T>
Solver<Dim, T>::Solver() {
}

template <int Dim, typename T>
Solver<Dim, T>::Solver(PDE<Dim, T> p, Mesh<Dim, Dim + 1, T> *m) {
	pde = p;
	mesh = m;
	max_index = mesh->reset_indices();
}

template <int Dim, typename T>
MatrixXd Solver<Dim, T>::get_stiffness_matrix(int unique_nodes) const {
	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	//cout << "max_index" <<max_index << endl;
	MatrixXd stiffness = MatrixXd::Zero(max_index+1, max_index+1);
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
			}
		}
		iter = iter->next;
	}
	return stiffness;
}

template <int Dim, typename T>
VectorXd Solver<Dim, T>::get_vector_part() const {
	MeshNode <Element<Dim, Dim + 1, T> >* iter = mesh->get_top_mesh_node();
	//cout << "max_index" <<max_index << endl;
	VectorXd f_vec = VectorXd::Zero(max_index + 1);
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
VectorXd Solver<Dim, T>::solve() {
	MatrixXd stiffness = get_stiffness_matrix(max_index);
	VectorXd f_vec = get_vector_part();
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
