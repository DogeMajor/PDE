#ifndef BOUNDARYELEMENTADDER_H
#define BOUNDARYELEMENTADDER_H
#include "Element.h"
#include "ElementFactory.h"
#include "Function.h"

template <int Dim, int N, typename T>
class BoundaryElementAdder {

public:
	BoundaryElementAdder(BoundaryConditions<T> b);
	~BoundaryElementAdder() {}
	map<array<int, 2>, Vertex<Dim, T>* > get_mid_vertices_map(Element<Dim, N, T> &el, map<array<int, 2>, Vertex<Dim, T>* > &all_boundary_vertices);
	vector<Element <Dim, N, T>* > add(Element <Dim, N, T>& el, map< array<int, 2>, Vertex<Dim, T>* > &commons);
	Element<Dim, N, T> get_corner_element(int I, vector <Vertex <Dim, T>* >  midpoint_vertices, Element <Dim, N, T>& el);
	Element<Dim, N, T> get_inner_element(int I, vector <Vertex <Dim, T>* >  midpoint_vertices);

private:
	ElementFactory <Dim, N, T> factory;
	BoundaryConditions<T> boundaries;

};

template <int Dim, int N, typename T>
BoundaryElementAdder<Dim, N, T>::BoundaryElementAdder(BoundaryConditions<T> b) {
	factory = ElementFactory <Dim, N, T>();
	boundaries = b;
}

template <int Dim, int N, typename T>
map<array<int, 2>, Vertex<Dim, T>* > BoundaryElementAdder<Dim, N, T>::get_mid_vertices_map(Element<Dim, N, T> &el, map<array<int, 2>, Vertex<Dim, T>* > &all_boundary_vertices) {

}

template <int Dim, int N, typename T>
vector<Element <Dim, N, T>* > BoundaryElementAdder<Dim, N, T>::add(Element <Dim, N, T>& el, map< array<int, 2>, Vertex<Dim, T>* > &commons) {

}

template <int Dim, int N, typename T>
Element<Dim, N, T> BoundaryElementAdder<Dim, N, T>::get_corner_element(int I, vector <Vertex <Dim, T>* >  midpoint_vertices, Element <Dim, N, T>& el) {

}

template <int Dim, int N, typename T>
Element<Dim, N, T> BoundaryElementAdder<Dim, N, T>::get_inner_element(int I, vector <Vertex <Dim, T>* >  midpoint_vertices) {

}


#endif