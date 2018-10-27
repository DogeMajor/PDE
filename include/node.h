#ifndef NODE_H
#define NODE_H
#include <iostream>
#include "Counter.h"
using namespace std;


template <int Dim, typename T>
class Node: Counter<Node<Dim, T> >{

public:
    Node();
    Node(const T &loc);
    Node(const Node &a);//copy constructor
    ~Node();
    void set_index(int ind);
    void set_shared_elements(int shared_els);
    const T& get_location() const;
    int get_index() const;
    int get_shared_elements() const;
	int how_many() const;
    Node<Dim,T>& operator=(const Node &a);
    bool operator== (const Node<Dim, T> &a) const;
    bool operator!=(const Node<Dim, T> &a) const;
    void show() const;

private:
    T location;
    int index;
    int shared_elements;
};


template <int Dim, typename T>
Node<Dim, T>::Node(){
    shared_elements = 0;
    index = -1;
}

template <int Dim, typename T>
Node<Dim, T>::Node(const T &loc){
	//for (int i = 0; i < Dim; i++) {location[i] = loc[i];}
    location = loc;
    shared_elements = 0;
    index = -1;
}

template <int Dim, typename T>
Node<Dim,T>::Node(const Node &a){
    location = a.location;
    index = a.index;
    shared_elements = a.shared_elements;
}

template <int Dim, typename T>
Node<Dim, T>::~Node(){
}

template <int Dim, typename T>
void Node<Dim,T>::set_index(int ind){
    index = ind;
}

template <int Dim, typename T>
void Node<Dim,T>::set_shared_elements(int shared_els){
    shared_elements = shared_els;
}

template <int Dim, typename T>
const T& Node<Dim,T>::get_location() const{
    return location;
}

template <int Dim, typename T>
int Node<Dim,T>::get_index() const{
    return index;
}

template <int Dim, typename T>
int Node<Dim,T>::get_shared_elements() const{
    return shared_elements;
}

template <int Dim, typename T>//Counting alive instances only
int Node<Dim, T>::how_many() const {
	return objects_alive;
}

template <int Dim, typename T>
Node<Dim,T>& Node<Dim,T>::operator=(const Node &a){
    if(*this!=a){
        location = a.location;
        index = a.index;
        shared_elements = a.shared_elements;
    }
    return *this;
}

template <int Dim, typename T>
bool Node<Dim,T>::operator==(const Node<Dim, T> &a) const{
    bool same_location = (location == a.location);
    bool same_index = (index == a.index);
    bool same_shared_elements = (shared_elements == a.shared_elements);
    return same_location && same_index && same_shared_elements;
}

template <int Dim, typename T>
bool Node<Dim,T>::operator!=(const Node<Dim, T> &a) const{
    return !(*this==a);
}

template <int Dim, typename T>//Does not show location when T = VectorXd!!!
void Node<Dim,T>::show() const{
    cout <<"index: " << index << endl;
    cout <<"location: " << endl;
    if(index != -1){
        T loc = get_location();
        for(int i=0; i<Dim; i++){cout << location[i] << ", ";}
		cout << endl;
    }
	
    //cout <<"Amount of shared elements: " << shared_elements << endl;
    //cout <<"Amount of all nodes: " << how_many() << endl;
}

#endif
