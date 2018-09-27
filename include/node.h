#ifndef NODE_H
#define NODE_H
#include <iostream>
using namespace std;


template <int Dim, typename T>
class Node{

public:
    Node();
    Node(T &loc);
    Node(const Node &a);//copy constructor
    void set_index(int ind);
    void set_shared_elements(int shared_els);
    T get_location() const;
    int get_index() const;
    int get_shared_elements() const;
    Node<Dim,T>& operator=(const Node &a);
    bool operator== (const Node &a) const;
    bool operator!=(const Node &a) const;
    void show() const;

private:
    T location;
    int index;
    int shared_elements;
};

template <int Dim, typename T>
Node<Dim, T>::Node(){
    shared_elements = 0;
    index = 0;
}

template <int Dim, typename T>
Node<Dim, T>::Node(T &loc){
    location = loc;
    shared_elements = 0;
    index = 0;
}

template <int Dim, typename T>
Node<Dim,T>::Node(const Node &a){
    location = a.location;
    index = a.index;
    shared_elements = a.shared_elements;
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
    bool same_location = bool(this->get_location() == a.get_location());
    bool same_index = bool(this->get_index() == a.get_index());
    bool same_shared_elements = (this->get_shared_elements() == a.get_shared_elements());
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
}

#endif
