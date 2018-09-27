#ifndef NODE_H
#define NODE_H
#include <iostream>
using namespace std;


template <int Dim, typename T>
class Node{

public:
    Node(T &loc);
    Node(const Node &a);//copy constructor
    void set_index(int ind);
    void set_neighbour_amount(int neighbours);
    T get_location() const;
    int get_index() const;
    int get_neighbour_amount() const;
    Node<Dim,T>& operator=(const Node &a);
    bool operator== (const Node &a) const;
    bool operator!=(const Node &a) const;
    void show() const;

private:
    T location;
    int index;
    int neighbour_amount;
    int dimension;
};

template <int Dim, typename T>
Node<Dim, T>::Node(T &loc){
    location = loc;
    neighbour_amount = 0;
    index = 0;
    dimension = Dim;
}

template <int Dim, typename T>
Node<Dim,T>::Node(const Node &a){
    location = a.location;
    index = a.index;
    neighbour_amount = a.neighbour_amount;
}

template <int Dim, typename T>
void Node<Dim,T>::set_index(int ind){
    index = ind;
}

template <int Dim, typename T>
void Node<Dim,T>::set_neighbour_amount(int neighbours){
    neighbour_amount = neighbours;
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
int Node<Dim,T>::get_neighbour_amount() const{
    return neighbour_amount;
}

template <int Dim, typename T>
Node<Dim,T>& Node<Dim,T>::operator=(const Node &a){
    if(*this!=a){
        location = a.location;
        index = a.index;
        neighbour_amount = a.neighbour_amount;
    }
    return *this;
}

template <int Dim, typename T>
bool Node<Dim,T>::operator==(const Node &a) const{
    bool same_location = bool(this->get_location() == a.get_location());
    bool same_index = bool(this->get_index() == a.get_index());
    bool same_neighbour_amount = (this->get_neighbour_amount() == a.get_neighbour_amount());
    return same_location && same_index && same_neighbour_amount;

}

template <int Dim, typename T>
bool Node<Dim,T>::operator!=(const Node &a) const{
    return !(*this==a);
}


template <int Dim, typename T>
void Node<Dim,T>::show() const{
    cout <<"index: " << index << endl;
    cout <<"dimension: " << dimension << endl;
    cout <<"location: " << endl;
    cout << location << endl;
    cout <<"Amount of neighbours: " << neighbour_amount << endl;
}

#endif
