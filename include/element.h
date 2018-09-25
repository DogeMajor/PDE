#ifndef ELEMENT_H
#define ELEMENT_H
#include <iostream>
#include "node.h"
using namespace std;


//We simply want to use the already existing nodes and don't need to worry about garbage collection.
template <int N, typename T>
class Element{

public:
    Element();
    Element(Node<T>* nod[N]);
    Element(Element &el);
    ~Element();
    Node<T> operator[](int i) const;
    Element<N,T>& operator=(const Element &el);
    bool operator==(const Element &el) const;
    bool operator!=(const Element &el) const;
    void show() const;

private:
    Node<T>* nodes[N];

};

template <int N, typename T>
Element<N,T>::Element(){
    for(int i=0; i<N; i++){
        nodes[i] = new Node<T>;
    }
}

template <int N, typename T>
Element<N,T>::Element(Node<T>* nod[N]){
    for(int i=0; i<N; i++){
        nodes[i] = nod[i];
    }
}

template <int N, typename T>
Element<N,T>::Element(Element &el){
    for(int i=0; i<N; i++){
        nodes[i] = el.nodes[i];
    }
}

template <int N, typename T>
Element<N,T>::~Element(){
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

template <int N, typename T>
Node<T> Element<N,T>::operator[](int i) const{
    return *(nodes[i]);
}

template <int N, typename T>
Element<N,T>& Element<N,T>::operator=(const Element &el){
    if(*this != el){
        for(int i=0; i<N; i++){
            nodes[i] = el.nodes[i];
        }
    }
    return *this;
}


template <int N, typename T>
bool Element<N,T>::operator==(const Element &el) const{
    bool result = true;

    for(int i=0; i<N; i++){
        result = result && bool(*(nodes[i]) == *(el.nodes[i]));
        }
    return result;
}

template <int N, typename T>
bool Element<N,T>::operator!=(const Element &el) const{
    return !(*this == el);
}

template <int N, typename T>
void Element<N,T>::show() const{
    int n = N;
    cout <<"#elements: " << n << endl;
    for(int i=0; i<N; i++){
        nodes[i]->show();
    }
}

#endif

