#ifndef MESH_H
#define MESH_H
#include <iostream>
#include "node.h"
#include "element.h"
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
	bool push(Element<Dim, N, T> &t);//To the top  OK
	bool pop();//From the top  OK
	bool push(MeshNode<Element<Dim, N, T> > *previous, Element<Dim, N, T> &t);//Adds efter the previous lement!
	bool pop(MeshNode<Element<Dim, N, T> > *previous);//Deletes the next element!!!
	int how_many() const;// { return objects_alive; }
	int how_many_nodes() const;
	//vector <T> divide_element(T &el);
	MeshNode<Element<Dim, N, T> > * get_top_mesh_node() { return top; }
	const Element<Dim, N, T> get_top() { return top->data; }
	const Element<Dim, N, T> get_last();
	void refine();
	bool operator==(const Mesh<Dim, N, T> & m) const;
	bool operator!=(const Mesh<Dim, N, T> & m) const;
	void show() const;

private:
	MeshNode<Element<Dim, N, T> > *top;
	ElementDivider<Dim, N, T> divider;
	int node_counter;

};

template <int Dim, int N, typename T>
Mesh<Dim, N, T>::Mesh() {
	top = nullptr;
	node_counter = 0;
}

template <int Dim, int N, typename T>
Mesh<Dim, N, T>::Mesh(Element<Dim, N, T> &t) {
	top = new MeshNode <Element<Dim, N, T> >{ t, nullptr };
	node_counter = 1;
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
bool Mesh<Dim, N, T>::pop() {
	if (top != nullptr) {
		MeshNode<Element<Dim, N, T> > temp = *top;
		delete top;
		top = temp.next;
		node_counter--;
		return true;
	}
	return false;
}

template <int Dim, int N, typename T>
bool Mesh<Dim, N, T>::pop(MeshNode<Element<Dim, N, T> > *previous) {
	if (previous->next != nullptr && previous != nullptr) {
		MeshNode<Element<Dim, N, T> >* temp = previous->next;
		previous->next = temp->next;
		delete temp;
		node_counter--;
		return true;
	}
	return false;
}

template <int Dim, int N, typename T>
bool Mesh<Dim, N, T>::push(MeshNode<Element<Dim, N, T> > *previous, Element<Dim, N, T> &t) {
	MeshNode <Element<Dim, N, T> >* temp = previous->next;
	previous->next = new MeshNode <Element<Dim, N, T> >{ t, temp };
	node_counter++;
	return true;
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
const Element<Dim,N,T> Mesh<Dim, N,T>::get_last() {
	if (top->next == nullptr) { return top->data; }
	MeshNode <Element<Dim, N, T> >* iter = top;
	while (iter->next != nullptr) { iter = iter->next; }
	return iter->data;
}

template <int Dim, int N, typename T>
void Mesh<Dim, N, T>::refine() {
	MeshNode<Element<Dim, N, T> >* previous = top;//For future needs...

	MeshNode<Element<Dim, N, T> >* iter = top->next;;//For future needs...
	vector<Element<Dim, N, T>* > new_els = divider.divide(top->data);
	pop();//Delete the first node - the next old node has to be deleted with some other method
	push(*new_els[0]);
	previous = top;
	for (int i = 1; i < new_els.size(); i++) {
		push(*new_els[i]);
	}
	//while (iter != nullptr) {
		//new_els = divider.divide(iter->data);
	//}

}

template <int Dim, int N, typename T>
void Mesh<Dim, N, T>::show() const {
	//BaseMesh<int, N T>* next_BaseMesh = get_next();
	cout << "Amount of nodes: " << node_counter << endl;
	MeshNode<Element<Dim, N, T> >* iter = top;
	int count = 0;
	while (iter != nullptr) {
		cout << "Showing elemnt no " << count << endl;
		iter->data.show();
		iter = iter->next;
		count++;
	}
}



#endif


