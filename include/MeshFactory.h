#ifndef MESHFACTORY_H
#define MESHFACTORY_H
#include <iostream>
#include <map>
#include <array>
#include "Vertex.h"
#include "Element.h"
#include "Mesh.h"

using namespace Eigen;
using namespace std;

template<int Dim>
class SimplexDivisions {
public:
	SimplexDivisions();
	~SimplexDivisions() {}
	array<int, Dim + 1> simplex_indices(int n) { return to_simplex[n]; }

private:
	map<int, array<int, Dim+1> > to_simplex;

};


template<>
class SimplexDivisions<2> {
public:
	SimplexDivisions();
	~SimplexDivisions() {}
	array<int, 3> simplex_indices(int n) { return to_simplex[n]; }

private:
	map<int, array<int, 3> > to_simplex;

};

SimplexDivisions<2>::SimplexDivisions() {
	to_simplex[0] = array<int, 3>({ 0,1,3 });
	to_simplex[1] = array<int, 3>({ 1,2,3 });
}


template<>
class SimplexDivisions<3> {
public:
	SimplexDivisions();
	~SimplexDivisions() {}
	array<int, 4> simplex_indices(int n) { return to_simplex[n]; }

private:
	map<int, array<int, 4> > to_simplex;

};

SimplexDivisions<3>::SimplexDivisions() {
	to_simplex[0] = array<int, 4>({ 1,2,0,4 });
	to_simplex[1] = array<int, 4>({ 1,2,4,5 });
	to_simplex[2] = array<int, 4>({ 1,2,3,5 });
	to_simplex[3] = array<int, 4>({ 5,6,2,4 });
	to_simplex[4] = array<int, 4>({ 5,6,3,2 });
	to_simplex[5] = array<int, 4>({ 5,6,3,7 });
}


template<int Dim>
MatrixXd get_box_coordinates(VectorXd midpoint, VectorXd lengths) {};

template<>
MatrixXd get_box_coordinates<2>(VectorXd midpoint, VectorXd lengths) {
	MatrixXd box_coords(4, 2);
	VectorXd corner(2);
	VectorXd delta(2);
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			delta[0] = -0.5*pow(-1, i)*lengths[0];
			delta[1] = -0.5*pow(-1, j)*lengths[1];
			corner = midpoint + delta;
			box_coords.row(i + j * 2) = corner;
		}
	}
	return box_coords;
}

template<>
MatrixXd get_box_coordinates<3>(VectorXd midpoint, VectorXd lengths) {
	MatrixXd box_coords(8, 3);
	VectorXd corner(3);
	VectorXd delta(3);
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				delta[0] = -0.5*pow(-1, i)*lengths[0];
				delta[1] = -0.5*pow(-1, j)*lengths[1];
				delta[2] = -0.5*pow(-1, k)*lengths[2];
				corner = midpoint + delta;
				box_coords.row(i + 2*j + 4*k) = corner;
			}
		}
	}
	return box_coords;
}


template<int Dim, int N, typename T>
class MeshFactory {

private:
	SimplexDivisions<Dim> simplex_divisions;
	VertexFactory<Dim, T> vertex_factory;
	ElementFactory<Dim, N, T> element_factory;

public:
	MeshFactory() {
		simplex_divisions = SimplexDivisions<Dim>();
		vertex_factory = VertexFactory<Dim, T>();
		element_factory = ElementFactory<Dim, N, T>();
	}
	~MeshFactory() {}
	T to_location(VectorXd coords);
	vector <Vertex<Dim, T>* > build_box_vertices(T mid_point, VectorXd lengths);
	vector <Vertex<Dim, T>* > get_simplex_vertices(int i, vector <Vertex<Dim, T>* > box_vertices);
	vector < Element<Dim, N, T> > cover_box_with_simplices(T mid_point, VectorXd lengths);
	Element<Dim, N, T> build_simplex(int i, vector <Vertex<Dim, T>* > box_vertices);
	Mesh<Dim, N, T> build_mesh(vector <Vertex<Dim, T>* > box_vertices);
};

template<int Dim, int N, typename T>//Add non-const access method [] in Point class if needed
T MeshFactory<Dim, N, T>::to_location(VectorXd coords) {
	T loc(Dim);
	for (int i = 0; i < coords.size(); i++) {
		loc[i] = coords[i];
	}
	return loc;
}

template<int Dim, int N, typename T>
vector <Vertex<Dim, T>* > MeshFactory<Dim, N, T>::build_box_vertices(T mid_point, VectorXd lengths) {
	vector <Vertex<Dim, T>* > box_vertices;
	VectorXd coord_vec(Dim);
	T loc(Dim);
	int vert_no = pow(2, Dim);

	MatrixXd box_coords = get_box_coordinates<Dim>(mid_point, lengths);
	for (int i = 0; i < vert_no; i++) {
		for (int j = 0; j < Dim; j++) {
			coord_vec[j] = box_coords(i, j);
		}
		loc = to_location(coord_vec);
		box_vertices.push_back(new Vertex<Dim, T>(loc));
	}
	return box_vertices;
}

template<int Dim, int N, typename T>
vector <Vertex<Dim, T>* > MeshFactory<Dim, N, T>::get_simplex_vertices(int i, vector <Vertex<Dim, T>* > box_vertices) {
	vector <Vertex<Dim, T>* > simplex_vertices(Dim + 1, nullptr);
	array<int, Dim + 1> indices = simplex_divisions.simplex_indices(i);
	for (int j = 0; j < Dim + 1; j++) {
		simplex_vertices[j] = box_vertices[indices[j]];
	}
	return simplex_vertices;
}

template<int Dim, int N, typename T>
vector <Element<Dim, N, T> > MeshFactory<Dim, N, T>::cover_box_with_simplices(T mid_point, VectorXd lengths) {
	vector <Element<Dim, N, T> > simplices(Dim * (Dim + 1) / 2);
	vector <Vertex<Dim, T>* > box_vertices = build_box_vertices(mid_point, lengths);
	vector<Vertex<Dim, T>* > simplex_verts;
	Element<Dim, N, T>  el;
	int max_ind = Dim * (Dim + 1) / 2;
	for (int i = 0; i < max_ind; i++) {
		simplex_verts = get_simplex_vertices(i, box_vertices);
		el = element_factory.build(simplex_verts);
		simplices[i] = element_factory.build(simplex_verts);
	}
	return simplices;
}

template<int Dim, int N, typename T>
Element<Dim, N, T> MeshFactory<Dim, N, T>::build_simplex(int i, vector <Vertex<Dim, T>* > box_vertices) {
	vector<Vertex<Dim, T>* > simplex_verts = get_simplex_vertices(i, box_vertices);
	return element_factory.build(simplex_verts);
}

template<int Dim, int N, typename T>
Mesh<Dim, N, T> MeshFactory<Dim, N, T>::build_mesh(vector <Vertex<Dim, T>* > box_vertices) {
	Mesh<Dim, N, T> mesh;
	int max = Dim * (Dim + 1) / 2;
	for(int i = 0; i < max; i++) {
		mesh.push(build_simplex(i, box_vertices));
		cout << "Element no " << i <<endl;
		mesh.get_top().show();
	}

	return mesh;
}

#endif
