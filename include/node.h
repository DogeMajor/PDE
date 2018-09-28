#ifndef NODE_H
#define NODE_H
#include <iostream>
#include "baseNode.h"
using namespace std;



template <int Dim, typename T>
class Node: public BaseNode{

public:
    Node();
    Node(T &loc);
    Node(const Node &a);//copy constructor
    ~Node();
    void set_index(int ind);
    void set_shared_elements(int shared_els);
    T get_location() const;
    int get_index() const;
    int get_shared_elements() const;
    int get_node_amount() const;
    Node<Dim,T>& operator=(const Node &a);
    bool operator== (const Node &a) const;
    bool operator!=(const Node &a) const;
    virtual void show() const;

private:
    T location;
    int index;
    int shared_elements;
    static int node_amount;
};

template <int Dim, typename T>
int Node<Dim, T>::node_amount=0;

template <int Dim, typename T>
Node<Dim, T>::Node() : BaseNode(){
    shared_elements = 0;
    index = 0;
    node_amount++;
}

template <int Dim, typename T>
Node<Dim, T>::Node(T &loc) : BaseNode(){
    location = loc;
    shared_elements = 0;
    index = 0;
    node_amount++;
}

template <int Dim, typename T>
Node<Dim,T>::Node(const Node &a) : BaseNode(){
    location = a.location;
    index = a.index;
    shared_elements = a.shared_elements;
    node_amount++;
}

template <int Dim, typename T>
Node<Dim, T>::~Node(){
    node_amount--;
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
T Node<Dim,T>::get_location() const{
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

template <int Dim, typename T>
int Node<Dim,T>::get_node_amount() const{
    return node_amount;
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
bool Node<Dim,T>::operator==(const Node &a) const{
    bool same_location = (location == a.location);
    bool same_index = (index == a.index);
    bool same_shared_elements = (shared_elements == a.shared_elements);
    return same_location && same_index && same_shared_elements;

}

template <int Dim, typename T>
bool Node<Dim,T>::operator!=(const Node &a) const{
    return !(*this==a);
}


template <int Dim, typename T>
void Node<Dim,T>::show() const{
    cout <<"index: " << index << endl;
    cout <<"location: " << endl;
    cout << location << endl;
    cout <<"Amount of shared elements: " << shared_elements << endl;
    cout <<"Amount of all nodes: " << node_amount << endl;
}

#endif
