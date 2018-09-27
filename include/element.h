#ifndef ELEMENT_H
#define ELEMENT_H
#include <iostream>
#include "node.h"
using namespace std;


//We simply want to use the already existing nodes and don't need to worry about garbage collection.
template <int Dim, int N, typename T>
class Element{

public:
    Element();
    Element(Node<Dim,T>* nod[N]);
    Element(Element &el);
    ~Element();
    void set_neighbours();
    void set_indices();
    Node<Dim,T> operator[](int i);
    Element<Dim,N,T>& operator=(const Element &el);
    bool operator==(const Element &el) const;
    bool operator!=(const Element &el) const;
    void show() const;

private:
    Node<Dim,T>* nodes[N];
    int dimension;

};

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(){
    for(int i=0; i<N; i++){
        nodes[i] = new Node<Dim,T>;
    }
    dimension = Dim;
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(Node<Dim,T>* nod[N]){
    for(int i=0; i<N; i++){
        nodes[i] = nod[i];
    }
    dimension = Dim;
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(Element &el){
    for(int i=0; i<N; i++){
        nodes[i] = el.nodes[i];
    }
    dimension = Dim;
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::~Element(){
    int neighbours = 0;
    for(int i=0; i<N; i++){
        neighbours = nodes[i]->get_neighbour_amount();
        nodes[i]->set_neighbour_amount(neighbours-1);
        if(neighbours-1 <= 0){
            delete nodes[i];
        }
    }
    cout << "Element destroyed!" << endl;
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::set_neighbours(){
    int neighbours = 0;
    for(int i=0; i<N; i++){
        neighbours = nodes[i]->get_neighbour_amount();
        nodes[i]->set_neighbour_amount(neighbours+N-1);
    }
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::set_indices(){
    for(int i=0; i<N; i++){
        nodes[i]->set_index(i);
    }
}

template <int Dim, int N, typename T>
Node<Dim,T> Element<Dim,N,T>::operator[](int i){
    return *(nodes[i]);
}

template <int Dim, int N, typename T>
Element<Dim,N,T>& Element<Dim,N,T>::operator=(const Element &el){
    if(*this != el){
        for(int i=0; i<N; i++){
            nodes[i] = el.nodes[i];
        }
    }
    dimension = el.dimension;
    return *this;
}


template <int Dim, int N, typename T>
bool Element<Dim,N,T>::operator==(const Element &el) const{
    bool result = true;

    for(int i=0; i<N; i++){
        result = result && bool(*(nodes[i]) == *(el.nodes[i]));
        }
    return result;
}

template <int Dim, int N, typename T>
bool Element<Dim,N,T>::operator!=(const Element &el) const{
    return !(*this == el);
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::show() const{
    int n = N;
    cout <<"#elements: " << n << endl;
    for(int i=0; i<N; i++){
        nodes[i]->show();
    }
}

#endif

