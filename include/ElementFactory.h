#ifndef ELEMENTFACTORY_H
#define ELEMENTFACTORY_H

#include "element.h"


template <int Dim, int N, typename T>
class ElementFactory {

public:
	ElementFactory() {}
	~ElementFactory() {}
	vector <Node <Dim, T> > build_nodes(vector <T> &locations);
	Element<Dim, N, T> build(vector <T> locations);
	Element<Dim, N, T> build(vector <Node <Dim, T>* > nodes_vec);
	MatrixXd get_inv_matrix(vector <Node <Dim, T>* > nodes_vec);
	SimplexFunction<T> build_function(MatrixXd M, int node_no);
	vector <SimplexFunction <T> > build_functions(vector <Node <Dim, T>* > nodes_vec);
};

template <int Dim, int N, typename T>
vector <Node <Dim, T> > ElementFactory<Dim, N, T>::build_nodes(vector <T> &locations) {
	vector <Node <Dim, T> > nodes_vec;
	NodeFactory<Dim, T> node_factory;
	for (int i = 0; i < locations.size(); i++) {
		nodes_vec.push_back(node_factory.build(locations[i]));
	}
	return nodes_vec;
}

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementFactory<Dim, N, T>::build(vector <T> locations) {
	vector <Node <Dim, T>* > nodes_vec(N, nullptr);
	for (int i = 0; i < locations.size(); i++) {
		nodes_vec[i] = new Node<Dim, T>(locations[i]);
	}
	vector <SimplexFunction <T> > funcs = build_functions(nodes_vec);
	Element<Dim, N, T> el(nodes_vec, funcs);
	return el;
}

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementFactory<Dim, N, T>::build(vector <Node <Dim, T>* > nodes_vec) {
	vector <SimplexFunction <T> > funcs = build_functions(nodes_vec);
	Element<Dim, N, T> el(nodes_vec, funcs);
	return el;
}

template <int Dim, int N, typename T>
MatrixXd ElementFactory<Dim, N, T>::get_inv_matrix(vector <Node <Dim, T>* > nodes_vec) {
	MatrixXd M(Dim + 1, Dim + 1);
	for (int node = 0; node < Dim + 1; node++) {
		for (int col = 0; col < Dim; col++) {
			M(node, col) = nodes_vec[node]->get_location()[col];
		}
		M(node, Dim) = 1;
	}
	return M.inverse();
}

template <int Dim, int N, typename T>
SimplexFunction<T> ElementFactory<Dim, N, T>::build_function(MatrixXd M, int node_no) {
	VectorXd unit_vec = VectorXd::Zero(N);
	unit_vec(node_no) = 1;
	VectorXd coeffs = M * unit_vec;
	SimplexFunction<T> fn_obj;
	fn_obj.coeff = coeffs;
	return fn_obj;
}

template <int Dim, int N, typename T>
vector <SimplexFunction <T> > ElementFactory<Dim, N, T>::build_functions(vector <Node <Dim, T>* > nodes_vec) {
	MatrixXd M_inv = get_inv_matrix(nodes_vec);
	vector <SimplexFunction <T> > functions;
	for (int i = 0; i < Dim + 1; i++) {
		functions.push_back(build_function(M_inv, i));
	}
	return functions;
}

#endif
