#ifndef MESH_H
#define MESH_H
#include <iostream>
#include "node.h"
#include "element.h"
#include "ElementFactory.h"
#include "ElementDivider.h"
#include "Function.h"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include <math.h>

using namespace std;
using namespace Eigen;


template <typename Nodetype>
struct MeshNode {
	Nodetype data;
	MeshNode<Nodetype>* next;
};


template <int Dim, int N, typename T>
class Mesh : Counter<Mesh<Dim, N, T> > {

public:
	Mesh();
	Mesh(Element<Dim, N, T> &t);
	~Mesh();
	bool push(Element<Dim, N, T> &t);//To the top
	bool push(MeshNode<Element<Dim, N, T> > *previous, Element<Dim, N, T> &t);//Adds efter the previous element!
	bool pop();//From the top
	bool pop(MeshNode<Element<Dim, N, T> > *previous);//Deletes the next element!!!
	int how_many() const;// { return objects_alive; }
	int how_many_nodes() const;
	int get_max_inner_index() const { return max_inner_index; }
	int get_max_outer_index() const { return max_outer_index; }
	MeshNode<Element<Dim, N, T> > * get_top_mesh_node() { return top; }
	Element<Dim, N, T> get_top() { return top->data; }
	Element<Dim, N, T> get_last();
	Element<Dim, N, T> get_element(int item_no);
	void refine();
	int set_inner_and_init_outer_indices(int index, BoundaryConditions<T> boundaries);
	int set_outer_indices(int index, BoundaryConditions<T> boundaries);
	void reset_indices(BoundaryConditions<T> boundaries);
	bool operator==(const Mesh<Dim, N, T> & m) const;
	bool operator!=(const Mesh<Dim, N, T> & m) const;
	void show() const;

private:
	MeshNode<Element<Dim, N, T> > *top;
	ElementDivider<Dim, N, T> divider;
	int node_counter;
	int max_inner_index;
	int max_outer_index;

};

template <int Dim, int N, typename T>
Mesh<Dim, N, T>::Mesh() {
	top = nullptr;
	node_counter = 0;
	max_inner_index = -1;
	max_outer_index = -1;
}

template <int Dim, int N, typename T>
Mesh<Dim, N, T>::Mesh(Element<Dim, N, T> &t) {
	top = new MeshNode <Element<Dim, N, T> >{ t, nullptr };
	node_counter = 1;
	max_inner_index = -1;
	max_outer_index = -1;
}

template <int Dim, int N, typename T>
Mesh<Dim, N, T>::~Mesh() {
	while (pop() != false) {}
	cout << "Mesh destroyed!" << endl;
	cout << "Nodes left in the mesh: " << node_counter << endl;
}

template <int Dim, int N, typename T>
bool Mesh<Dim, N, T>::push(Element<Dim, N, T> &t) {
	top = new MeshNode <Element<Dim, N, T> >{ t, top };
	node_counter++;
	return true;
}

template <int Dim, int N, typename T>
bool Mesh<Dim, N, T>::push(MeshNode<Element<Dim, N, T> > *previous, Element<Dim, N, T> &t) {
	if (previous == nullptr) {return push(t);}
	MeshNode <Element<Dim, N, T> >* temp = previous->next;
	previous->next = new MeshNode <Element<Dim, N, T> >{ t, temp };
	node_counter++;
	return true;
}

template <int Dim, int N, typename T>
bool Mesh<Dim, N, T>::pop() {
	if (top != nullptr) {
		MeshNode<Element<Dim, N, T> >* old_top = top;
		top = old_top->next;
		delete old_top;
		node_counter--;
		return true;
	}
	return false;
}

template <int Dim, int N, typename T>
bool Mesh<Dim, N, T>::pop(MeshNode<Element<Dim, N, T> > *previous) {
	if (previous == nullptr) {return pop();}
	if (previous->next != nullptr) {
		MeshNode<Element<Dim, N, T> >* old = previous->next;
		previous->next = old->next;
		delete old;
		node_counter--;
		return true;
	}
	return false;
}

template <int Dim, int N, typename T>
int Mesh<Dim, N, T>::how_many() const {
	return objects_alive;
}

template <int Dim, int N, typename T>
int Mesh<Dim, N, T>::how_many_nodes() const {
	return node_counter;
}

template <int Dim, int N, typename T>
bool Mesh<Dim, N, T>::operator==(const Mesh<Dim, N, T> & m) const {
	bool same_top = (m.top == top);
	return same_top;
}

template <int Dim, int N, typename T>
bool Mesh<Dim, N, T>::operator!=(const Mesh<Dim, N, T> & m) const {
	return !(*this == m);
}

template <int Dim, int N, typename T>
Element<Dim,N,T> Mesh<Dim, N,T>::get_last() {
	if (top->next == nullptr) { return top->data; }
	MeshNode <Element<Dim, N, T> >* iter = top;
	while (iter->next != nullptr) { iter = iter->next; }
	return iter->data;
}

template <int Dim, int N, typename T>
Element<Dim, N, T> Mesh<Dim, N, T>::get_element(int item_no) {
	MeshNode <Element<Dim, N, T> >* iter = top;
	for (int i = 0; i < item_no; i++) { iter = iter->next; }
	return iter->data;
}

template <int Dim, int N, typename T>
void Mesh<Dim, N, T>::refine() {
	MeshNode<Element<Dim, N, T> >* original_node = top;
	MeshNode<Element<Dim, N, T> >* previous = top;
	MeshNode<Element<Dim, N, T> >* before_original_node = nullptr;
	vector<Element<Dim, N, T>* > new_els;
	map< array<int, 2>, Node<Dim, T>* > commons;

	while(original_node != nullptr){
		new_els = divider.divide(original_node->data, commons);
		for (int i = 0; i < new_els.size(); i++) {
			push(previous, *new_els[i]);
			previous = previous->next;
		}
		pop(before_original_node);		
		before_original_node = previous;
		original_node = previous->next;
		previous = previous->next;
	}
	map< array<int, 2>, Node<Dim, T>* >::iterator map_iter = commons.begin();
	commons.erase(map_iter, commons.end());
}

template <int Dim, int N, typename T>
int Mesh<Dim, N, T>::set_inner_and_init_outer_indices(int index, BoundaryConditions<T> boundaries) {
	MeshNode<Element<Dim, N, T> >* iter = top;
	while (iter != nullptr) {
		index = iter->data.set_inner_node_indices(index, boundaries);
		iter->data.set_outer_node_indices_to(-1, boundaries);
		iter = iter->next;
	}
	return index;
}

template <int Dim, int N, typename T>//Only works if run after set_inner_and_init_outer_indices!!
int Mesh<Dim, N, T>::set_outer_indices(int index, BoundaryConditions<T> boundaries) {
	MeshNode<Element<Dim, N, T> >* iter = top;
	while (iter != nullptr) {
		index = iter->data.set_indices(index);
		iter = iter->next;
	}
	return index;
}

template <int Dim, int N, typename T>
void Mesh<Dim, N, T>::reset_indices(BoundaryConditions<T> boundaries) {
	max_inner_index = set_inner_and_init_outer_indices(max_inner_index, boundaries);
	max_outer_index = set_outer_indices(max_inner_index, boundaries);
}

template <int Dim, int N, typename T>
void Mesh<Dim, N, T>::show() const {
	cout << "Amount of nodes: " << node_counter << endl;
	MeshNode<Element<Dim, N, T> >* iter = top;
	int count = 0;
	while (iter != nullptr) {
		cout << "Showing element no " << count << endl;
		iter->data.show();
		iter = iter->next;
		count++;
	}
	cout << endl;
}


#endif
