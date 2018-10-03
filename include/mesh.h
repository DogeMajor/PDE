#ifndef MESH_H
#define MESH_H
#include <iostream>
//#include "node.h"
#include "element.h"
#include <math.h>

using namespace std;
using namespace Eigen;


template <int Dim, int N, typename T>
class Mesh{

public:
    Mesh();
    Mesh(Element<Dim,N,T> &t);
    Mesh(Element<Dim,N,T> t, Mesh<Dim,N,T> *n);
    ~Mesh();
    const Mesh <Dim,N,T>& operator=(const Mesh<Dim,N,T> &m);
    bool operator==(const Mesh<Dim,N,T> &m) const;
    bool operator!=(const Mesh<Dim,N,T> &m) const;
    Element<Dim,N,T> get_element() const;
    const Mesh<Dim,N,T>* get_next();
    double weak_form_element(int i, int j) const;//Calculating the stiffness matrix
    void show() const;

private:
    Element <Dim,N,T> top;
    Mesh<Dim,N,T> *next;

};

template <int Dim, int N, typename T>
Mesh<Dim,N,T>::Mesh(){
    //top = new Element<Dim,N,T>;
    //top = new Element<Dim,N,T>;
    next = nullptr;
}

template <int Dim, int N, typename T>
Mesh<Dim,N,T>::Mesh(Element<Dim,N,T> &t){
    top = t;
    next = nullptr;
}


template <int Dim, int N, typename T>
Mesh<Dim,N,T>::Mesh(Element<Dim,N,T> t, Mesh<Dim,N,T> *n){
    top = t;
    next = n;
}

template <int Dim, int N, typename T>
Mesh<Dim,N,T>::~Mesh(){
    //delete top;
    //delete next;
    delete next;
    cout << "Mesh destroyed!" << endl;
}

template <int Dim, int N, typename T>
const Mesh <Dim,N,T>& Mesh<Dim,N,T>::operator=(const Mesh<Dim,N,T> & m){
    if(*this!=m){
        top = m.get_element();
        next = m.get_next();
    }
    return *this;
}

template <int Dim, int N, typename T>
bool Mesh<Dim,N,T>::operator==(const Mesh<Dim,N,T> & m) const{
    bool same_top = (m.top == top);
    bool same_next = (m.next == next);
    return same_top && same_next;
}

template <int Dim, int N, typename T>
bool Mesh<Dim,N,T>::operator!=(const Mesh<Dim,N,T> & m) const{
    return !(*this == m);
}

template <int Dim, int N, typename T>
Element<Dim,N,T> Mesh<Dim,N,T>::get_element() const{
    return top;
}

template <int Dim, int N, typename T>
const Mesh<Dim,N,T>* Mesh<Dim,N,T>::get_next(){
    return next;
}

template <int Dim, int N, typename T>
double Mesh<Dim,N,T>::weak_form_element(int i, int j) const{
    return 0.0;
}

template <int Dim, int N, typename T>
void Mesh<Dim,N,T>::show() const{
    //Mesh<Dim,N,T>* next_mesh = get_next();
    cout << "To be fixed..." << endl;
    top.show();
    /*do{
        top.show();
        this = next;
    }
    while(this != nullptr);*/
}

#endif


