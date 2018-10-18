#ifndef BASEMESH_H
#define BASEMESH_H
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
	bool push(T &t);//To the top
	bool pop();//From the top
	//LinkedMesh <T>& operator=(const LinkedMesh<T> &m);
	bool operator==(const LinkedMesh<T> &m) const;
	bool operator!=(const LinkedMesh<T> &m) const;
	int how_many() const;// { return objects_alive; }
	//vector <T> divide_element(T &el);
	const MeshNode<T> get_top() { return *top; }
	const MeshNode<T> get_last();
	void show() const;

private:
	MeshNode<T> *top;

};

template <typename T>
LinkedMesh<T>::LinkedMesh() {
	top = nullptr;
}

template <typename T>
LinkedMesh<T>::LinkedMesh(T &t) {
	top = new MeshNode<T>{ t, nullptr };
}

template <typename T>
LinkedMesh<T>::~LinkedMesh() {
	while (pop() != false) {}
	cout << "LinkedMesh destroyed!" << endl;
}

template <typename T>
bool LinkedMesh<T>::push(T &t) {
	top = new MeshNode<T>{ t, top };
	return true;
}

template <typename T>
bool LinkedMesh<T>::pop() {
	if (top != nullptr) {
		MeshNode<T> temp = *top;
		delete top;
		top = temp.next;
		cout << "top element destroyed!" << endl;
		return true;
	}
	return false;
}

template <typename T>
LinkedMesh<T>::LinkedMesh(const LinkedMesh<T> &mesh) {
	if (*this != mesh) {
		top = mesh.top;
		next = mesh.next;
	}
	return *this;
}

template <typename T>
int LinkedMesh<T>::how_many() const {
	return objects_alive;
}
/*template <typename T>
LinkedMesh <T>& LinkedMesh<T>::operator=(const LinkedMesh<T> & m) {
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
}*/

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
const MeshNode<T> LinkedMesh<T>::get_last(){
	MeshNode <T>* iter = top;
	while (iter->next != nullptr) {iter = next;}
	return *iter;
}

template <typename T>
void LinkedMesh<T>::show() const {
	//BaseMesh<T>* next_BaseMesh = get_next();
	cout << "To be fixed..." << endl;
	//top.show();
	/*do{
		top.show();
		this = next;
	}
	while(this != nullptr);*/
}

#endif