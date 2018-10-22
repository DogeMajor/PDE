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
	const Element<Dim, N, T> get_element(int item_no);
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
		MeshNode<Element<Dim, N, T> >* old_top = top;
		top = old_top->next;
		delete old_top;
		
		cout << "heh" << endl;
		//cout << "new top: " << endl;
		//if (top != nullptr) { top->data.show(); }
		node_counter--;
		return true;
	}
	return false;
}

template <int Dim, int N, typename T>//NOT OK!!!!
bool Mesh<Dim, N, T>::pop(MeshNode<Element<Dim, N, T> > *previous) {
	if (previous == nullptr) {
		//cout << "Deleting element" << endl;
		//top->data.show();
		return pop();
	}
	if (previous->next != nullptr) {
		MeshNode<Element<Dim, N, T> >* old = previous->next;
		//temp->data.show();
		//delete previous->next;
		previous->next = old->next;
		cout << "Deleting element" << endl;
		old->data.show();
		delete old;
		node_counter--;
		return true;
	}
	return false;
}

template <int Dim, int N, typename T>
bool Mesh<Dim, N, T>::push(MeshNode<Element<Dim, N, T> > *previous, Element<Dim, N, T> &t) {
	if (previous == nullptr) {
		return push(t);
	}
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
const Element<Dim, N, T> Mesh<Dim, N, T>::get_element(int item_no) {
	MeshNode <Element<Dim, N, T> >* iter = top;
	for (int i = 0; i < item_no; i++) { iter = iter->next; }
	return iter->data;
}

template <int Dim, int N, typename T>
void Mesh<Dim, N, T>::refine() {
	MeshNode<Element<Dim, N, T> >* original_node = top;
	MeshNode<Element<Dim, N, T> >* previous = nullptr;
	MeshNode<Element<Dim, N, T> >* before_original_node = nullptr;
	//original_node->data.show();
	vector<Element<Dim, N, T>* > new_els;
	//do{
	for (int j = 0; j < 2; j++) {
		new_els = divider.divide(original_node->data);
		//if (before_original_node != nullptr) { before_original_node->next->data.show(); }
		cout << "Deleting element succeeded: " << pop(before_original_node);
		for (int i = 0; i < new_els.size(); i++) {
			push(previous, *new_els[i]);
			previous = previous->next;
		}
		cout << "heh" << endl;
		show();
		cout << "heh" << endl;
		before_original_node = previous;
		previous->next->data.show();
		//pop(previous);
		original_node = previous->next;
		//if (before_original_node != nullptr) { pop(before_original_node); }
		
		cout << "Original node is nullptr:" << bool(original_node == nullptr) << endl;
		//original_node->data.show();
	}
	//}
	//while (original_node != nullptr);

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
	cout << endl;
}



#endif


