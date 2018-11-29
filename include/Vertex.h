#ifndef VERTEX_H
#define VERTEX_H
#include <iostream>
#include "Counter.h"
using namespace std;


template <int Dim, typename T>//It is assumed that class T has a function int size() const!
class Vertex: Counter<Vertex<Dim, T> >{

public:
    Vertex();
    Vertex(const T &loc);
    Vertex(const Vertex &a);
    ~Vertex();
    void set_index(int ind);
    void set_shared_elements(int shared_els);
	const T& get_location() const { return location; }
	int get_index() const { return index; }
	int get_shared_elements() const { return shared_elements; }
	int how_many() const;
    Vertex<Dim,T>& operator=(const Vertex &a);
    bool operator==(const Vertex<Dim, T> &a) const;
    bool operator!=(const Vertex<Dim, T> &a) const;
    void show() const;

private:
    T location;
    int index;
    int shared_elements;
};


template <int Dim, typename T>
Vertex<Dim, T>::Vertex(){
    shared_elements = 0;
    index = -1;
}

template <int Dim, typename T>
Vertex<Dim, T>::Vertex(const T &loc){
    location = loc;
    shared_elements = 0;
    index = -1;
}

template <int Dim, typename T>
Vertex<Dim,T>::Vertex(const Vertex &a){
    location = a.location;
    index = a.index;
    shared_elements = a.shared_elements;
}

template <int Dim, typename T>
Vertex<Dim, T>::~Vertex(){
}

template <int Dim, typename T>
void Vertex<Dim,T>::set_index(int ind){
    index = ind;
}

template <int Dim, typename T>
void Vertex<Dim,T>::set_shared_elements(int shared_els){
    shared_elements = shared_els;
}

template <int Dim, typename T>
int Vertex<Dim, T>::how_many() const {
	return objects_alive;
}

template <int Dim, typename T>
Vertex<Dim,T>& Vertex<Dim,T>::operator=(const Vertex &a){
    if(*this!=a){
        location = a.location;
        index = a.index;
        shared_elements = a.shared_elements;
    }
    return *this;
}

template <int Dim, typename T>
bool Vertex<Dim,T>::operator==(const Vertex<Dim, T> &a) const{
    bool same_location = (location == a.location);
    bool same_index = (index == a.index);
    bool same_shared_elements = (shared_elements == a.shared_elements);
    return same_location && same_index && same_shared_elements;
}

template <int Dim, typename T>
bool Vertex<Dim,T>::operator!=(const Vertex<Dim, T> &a) const{
    return !(*this==a);
}

template <int Dim, typename T>
void Vertex<Dim,T>::show() const{
    cout <<"index: " << index << endl;
    cout <<"location: " << endl;
    T loc = get_location();
	if(loc.size()==Dim){
		for (int i = 0; i < Dim; i++) { cout << loc[i] << ", "; }
	}
	cout << endl;
    cout <<"Amount of shared elements: " << shared_elements << endl;
}


template <int Dim, typename T>
class VertexFactory{

public:
	VertexFactory() {}
	~VertexFactory() {}
	Vertex<Dim, T> build(T loc);

};

template <int Dim, typename T>
Vertex<Dim, T>  VertexFactory<Dim, T>::build(T loc) {
	Vertex<Dim, T> Vertex(loc);
	return Vertex;
}


#endif
