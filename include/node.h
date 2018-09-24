#ifndef NODE_H
#define NODE_H
#include <iostream>
using namespace std;


template <int N, typename T>
class Node{

public:
    Node(T &loc);
    void set_index(int ind);
    void set_neighbour_amount(int neighbours);
    T get_location() const;
    int get_index() const;
    int get_neighbour_amount() const;
    void show() const;

private:
    T location;
    int index;
    int neighbour_amount;
};

template <int N,typename T>
Node<N,T>::Node(T &loc){
    location = loc;
    neighbour_amount = 0;
    index = 0;
}

template <int N,typename T>
void Node<N,T>::set_index(int ind){
    index = ind;
}

template <int N,typename T>
void Node<N,T>::set_neighbour_amount(int neighbours){
    neighbour_amount = neighbours;
}

template <int N,typename T>
T Node<N,T>::get_location() const{
    return location;
}

template <int N,typename T>
int Node<N,T>::get_index() const{
    return index;
}

template <int N,typename T>
int Node<N,T>::get_neighbour_amount() const{
    return neighbour_amount;
}


template <int N,typename T>
void Node<N,T>::show() const{
    cout <<"index: " << endl;
    cout << index << endl;
    cout <<"location: " << location << endl;
    cout <<"Amount of neighbours: " << neighbour_amount << endl;
}

#endif
