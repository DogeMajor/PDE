#ifndef SIMPLEMESH_H
#define SIMPLEMESH_H
#include <iostream>
//#include "node.h"
//#include "element.h"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include <math.h>

using namespace std;
using namespace Eigen;


template <typename T>
class SimpleMesh{

public:
    SimpleMesh();
    SimpleMesh(T &t);
    SimpleMesh(T &t, SimpleMesh<T> *n);
    ~SimpleMesh();
    const SimpleMesh <T>& operator=(const SimpleMesh<T> &m);
    bool operator==(const SimpleMesh<T> &m) const;
    bool operator!=(const SimpleMesh<T> &m) const;
    T get_element() const;
    const SimpleMesh<T>* get_next();
    void show() const;

private:
    T top;
    SimpleMesh *next;

};

template <typename T>
SimpleMesh<T>::SimpleMesh(){
    //top = new Element<T>;
    //top = new Element<T>;
    next = nullptr;
}

template <typename T>
SimpleMesh<T>::SimpleMesh(T &t){
    top = t;
    next = nullptr;
}


template <typename T>
SimpleMesh<T>::SimpleMesh(T &t, SimpleMesh<T> *n){
    top = t;
    next = n;
}

template <typename T>
SimpleMesh<T>::~SimpleMesh(){
    //delete top;
    //delete next;
    delete next;
    cout << "SimpleMesh destroyed!" << endl;
}

template <typename T>
const SimpleMesh <T>& SimpleMesh<T>::operator=(const SimpleMesh<T> & m){
    if(*this!=m){
        top = m.get_element();
        if(next!=nullptr){
            if(m.next != nullptr){
                next = m.next;
            }
            else{
                delete next;
                next = nullptr;
            }
        }

        else{
            if(m.next!=nullptr){
                next = new SimpleMesh(m.get_next->get_element());
            }
        }

    }
    return *this;
}

template <typename T>
bool SimpleMesh<T>::operator==(const SimpleMesh<T> & m) const{
    bool same_top = (m.top == top);
    bool same_next = (m.next == next);
    return same_top && same_next;
}

template <typename T>
bool SimpleMesh<T>::operator!=(const SimpleMesh<T> & m) const{
    return !(*this == m);
}

template <typename T>
T SimpleMesh<T>::get_element() const{
    return top;
}

template <typename T>
const SimpleMesh<T>* SimpleMesh<T>::get_next(){
    return next;
}

template <typename T>
void SimpleMesh<T>::show() const{
    //SimpleMesh<T>* next_SimpleMesh = get_next();
    cout << "To be fixed..." << endl;
    //top.show();
    /*do{
        top.show();
        this = next;
    }
    while(this != nullptr);*/
}

#endif



