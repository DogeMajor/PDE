#ifndef ELEMENTFACTORY_H
#define ELEMENTFACTORY_H

#include "Element.h"


template <int Dim, int N, typename T>
class ElementFactory {

public:
	ElementFactory() {}
	~ElementFactory() {}
	vector <Vertex <Dim, T> > build_vertices(vector <T> &locations);
	Element<Dim, N, T> build(vector <T> locations);
	Element<Dim, N, T> build(vector <Vertex <Dim, T>* > vertices_vec);
	MatrixXd get_inv_matrix(vector <Vertex <Dim, T>* > vertices_vec);
	SimplexFunction<T> build_function(MatrixXd M, int vertex_no);
	vector <SimplexFunction <T> > build_functions(vector <Vertex <Dim, T>* > vertices_vec);
};

template <int Dim, int N, typename T>
vector <Vertex <Dim, T> > ElementFactory<Dim, N, T>::build_vertices(vector <T> &locations) {
	vector <Vertex <Dim, T> > vertices_vec;
	VertexFactory<Dim, T> vertex_factory;
	for (int i = 0; i < locations.size(); i++) {
		vertices_vec.push_back(vertex_factory.build(locations[i]));
	}
	return vertices_vec;
}

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementFactory<Dim, N, T>::build(vector <T> locations) {
	vector <Vertex <Dim, T>* > vertices_vec(N, nullptr);
	for (int i = 0; i < locations.size(); i++) {
		vertices_vec[i] = new Vertex<Dim, T>(locations[i]);
	}
	vector <SimplexFunction <T> > funcs = build_functions(vertices_vec);
	Element<Dim, N, T> el(vertices_vec, funcs);
	return el;
}

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementFactory<Dim, N, T>::build(vector <Vertex <Dim, T>* > vertices_vec) {
	vector <SimplexFunction <T> > funcs = build_functions(vertices_vec);
	Element<Dim, N, T> el(vertices_vec, funcs);
	return el;
}

template <int Dim, int N, typename T>
MatrixXd ElementFactory<Dim, N, T>::get_inv_matrix(vector <Vertex <Dim, T>* > vertices_vec) {
	MatrixXd M(Dim + 1, Dim + 1);
	for (int vertex = 0; vertex < Dim + 1; vertex++) {
		for (int col = 0; col < Dim; col++) {
			M(vertex, col) = vertices_vec[vertex]->get_location()[col];
		}
		M(vertex, Dim) = 1;
	}
	return M.inverse();
}

template <int Dim, int N, typename T>
SimplexFunction<T> ElementFactory<Dim, N, T>::build_function(MatrixXd M, int vertex_no) {
	VectorXd unit_vec = VectorXd::Zero(N);
	unit_vec(vertex_no) = 1;
	VectorXd coeffs = M * unit_vec;
	SimplexFunction<T> fn_obj;
	fn_obj.coeff = coeffs;
	return fn_obj;
}

template <int Dim, int N, typename T>
vector <SimplexFunction <T> > ElementFactory<Dim, N, T>::build_functions(vector <Vertex <Dim, T>* > vertices_vec) {
	MatrixXd M_inv = get_inv_matrix(vertices_vec);
	vector <SimplexFunction <T> > functions;
	for (int i = 0; i < Dim + 1; i++) {
		functions.push_back(build_function(M_inv, i));
	}
	return functions;
}

#endif
