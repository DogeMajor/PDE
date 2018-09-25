#ifndef NODE_H
#define NODE_H
#include <iostream>
using namespace std;


template <typename T>
class Node{

public:
    Node(T &loc);
    void set_index(int ind);
    void set_neighbour_amount(int neighbours);
    T get_location() const;
    int get_index() const;
    int get_neighbour_amount() const;
    bool operator== (const Node &a) const;
    bool operator!=(const Node &a) const;
    void show() const;

private:
    T location;
    int index;
    int neighbour_amount;
};

template <typename T>
Node<T>::Node(T &loc){
    location = loc;
    neighbour_amount = 0;
    index = 0;
}

template <typename T>
void Node<T>::set_index(int ind){
    index = ind;
}

template <typename T>
void Node<T>::set_neighbour_amount(int neighbours){
    neighbour_amount = neighbours;
}

template <typename T>
T Node<T>::get_location() const{
    return location;
}

template <typename T>
int Node<T>::get_index() const{
    return index;
}

template <typename T>
int Node<T>::get_neighbour_amount() const{
    return neighbour_amount;
}


template <typename T>
bool Node<T>::operator==(const Node &a) const{
    bool same_location = bool(this->get_location() == a.get_location());
    bool same_index = bool(this->get_index() == a.get_index());
    bool same_neighbour_amount = (this->get_neighbour_amount() == a.get_neighbour_amount());
    return same_location && same_index && same_neighbour_amount;

}

template <typename T>
bool Node<T>::operator!=(const Node &a) const{
    return !(*this==a);
}


template <typename T>
void Node<T>::show() const{
    cout <<"index: " << endl;
    cout << index << endl;
    cout <<"location: " << location << endl;
    cout <<"Amount of neighbours: " << neighbour_amount << endl;
}

#endif
