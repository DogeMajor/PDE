#ifndef LINKEDMESH_H
#define LINKEDMESH_H
#include <iostream>
//#include "node.h"
//#include "element.h"
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


template <typename T>
class LinkedMesh : Counter<LinkedMesh<T> > {

public:
	LinkedMesh();
	LinkedMesh(T &t);
	LinkedMesh(const LinkedMesh<T> &mesh);
	~LinkedMesh();
	bool push(T &t);//To the top  OK
	bool pop();//From the top  OK
	LinkedMesh <T>& operator=(const LinkedMesh<T> &m);
	bool operator==(const LinkedMesh<T> &m) const;
	bool operator!=(const LinkedMesh<T> &m) const;
	int how_many() const;// { return objects_alive; }
	int how_many_nodes() const;
	//vector <T> divide_element(T &el);
	const T get_top() { return top->data; }
	const T get_last();
	void refine();
	void show() const;

private:
	MeshNode<T> *top;
	int node_counter;

};

template <typename T>
LinkedMesh<T>::LinkedMesh() {
	top = nullptr;
	node_counter = 0;
}

template <typename T>
LinkedMesh<T>::LinkedMesh(T &t) {
	top = new MeshNode<T>{ t, nullptr };
	node_counter = 1;
}

template <typename T>
LinkedMesh<T>::~LinkedMesh() {
	while (pop() != false) {}
	cout << "LinkedMesh destroyed!" << endl;
	cout << "Nodes left in the mesh: " << node_counter << endl;
}

template <typename T>
bool LinkedMesh<T>::push(T &t) {
	top = new MeshNode<T>{ t, top };
	node_counter++;
	return true;
}

template <typename T>
bool LinkedMesh<T>::pop() {
	if (top != nullptr) {
		MeshNode<T> temp = *top;
		delete top;
		top = temp.next;
		node_counter--;
		return true;
	}
	return false;
}

template <typename T>
LinkedMesh<T>::LinkedMesh(const LinkedMesh<T> &mesh) {//If the speed of copying matters you should use doubly linked list!
	//while (pop() != false) { cout << "Element popped out of the top" << endl; }
	//pop();
	cout << "COPY:  top destroyed!!" << endl;
	//top = mesh.top;
	push(mesh.top->data);
}

template <typename T>
int LinkedMesh<T>::how_many() const {
	return objects_alive;
}

template <typename T>
int LinkedMesh<T>::how_many_nodes() const {
	return node_counter;
}

template <typename T>
LinkedMesh <T>& LinkedMesh<T>::operator=(const LinkedMesh<T> & m) {
	cout << "ASSIGNMENT OP WAS CALLED!!!!" << endl;
	if (*this != m) {
		top = m.get_element();
		if (next != nullptr) {
			if (m.next != nullptr) { next = m.next; }
			else {
				delete next;
				next = nullptr;
			}
		}
		else {
			if (m.next != nullptr) { next = new LinkedMesh(m.get_next->get_element()); }
		}

	}
	return *this;
}

template <typename T>
bool LinkedMesh<T>::operator==(const LinkedMesh<T> & m) const {
	bool same_top = (m.top == top);
	return same_top;
}

template <typename T>
bool LinkedMesh<T>::operator!=(const LinkedMesh<T> & m) const {
	return !(*this == m);
}


template <typename T>
const T LinkedMesh<T>::get_last(){
	if (top->next == nullptr) { return top->data; }
	MeshNode <T>* iter = top;
	while (iter->next != nullptr) {iter = iter->next;}
	return iter->data;
}

template <typename T>
void LinkedMesh<T>::refine() {

}

template <typename T>
void LinkedMesh<T>::show() const {
	//BaseMesh<T>* next_BaseMesh = get_next();
	cout << "To be fixed..." << endl;
	cout << "Amount of nodes: " << node_counter << endl;
	//top.show();
	/*do{
		top.show();
		this = next;
	}
	while(this != nullptr);*/
}

#endif