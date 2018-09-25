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
    ~Element();
    Node<T> operator[](int i) const;
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
        //cout << "index: " << i << endl;
        nodes[i] = nod[i];
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
    //nodes = NULL;
    cout << "Element destroyed!" << endl;
}

template <int N, typename T>
Node<T> Element<N,T>::operator[](int i) const{
    return *(nodes[i]);
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

