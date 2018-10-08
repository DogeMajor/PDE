#ifndef BASEMESH_H
#define BASEMESH_H
#include <iostream>
//#include "node.h"
//#include "element.h"
#include "FunctionHandler.h"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include <math.h>

using namespace std;
using namespace Eigen;


template <typename T>
class BaseMesh{

public:
    BaseMesh();
    BaseMesh(T t);
    BaseMesh(T &t, BaseMesh<T> *n);
    ~BaseMesh();
    void set_top(T t);
    const BaseMesh <T>& operator=(const BaseMesh<T> &m);
    bool operator==(const BaseMesh<T> &m) const;
    bool operator!=(const BaseMesh<T> &m) const;
    T get_element() const;
    const BaseMesh<T>* get_next();
    void show() const;

private:
    T top;
    BaseMesh *next;
    vector <SimplexFunction <T> > functions;

};

template <typename T>
BaseMesh<T>::BaseMesh(){
    //top = new Element<T>;
    //top = new Element<T>;
    next = nullptr;
}

template <typename T>
BaseMesh<T>::BaseMesh(T t){
    top = t;
    next = nullptr;
}


template <typename T>
BaseMesh<T>::BaseMesh(T &t, BaseMesh<T> *n){
    top = t;
    next = n;
}

template <typename T>
BaseMesh<T>::~BaseMesh(){
    //delete top;
    //delete next;
    if(next!=nullptr){delete next;}
    cout << "BaseMesh destroyed!" << endl;
}

template <typename T>
void BaseMesh<T>::set_top(T t){
    cout << "Setting to top to value" << endl;
    //t.show();
    top = t;
}


template <typename T>
const BaseMesh <T>& BaseMesh<T>::operator=(const BaseMesh<T> & m){
    if(*this!=m){
        top = m.get_element();
        functions = m.functions;
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
                next = new BaseMesh(m.get_next->get_element());
            }
        }

    }
    return *this;
}

template <typename T>
bool BaseMesh<T>::operator==(const BaseMesh<T> & m) const{
    bool same_top = (m.top == top);
    bool same_next = (m.next == next);
    bool same_fns = (m.functions == functions);
    return same_top && same_next;
}

template <typename T>
bool BaseMesh<T>::operator!=(const BaseMesh<T> & m) const{
    return !(*this == m);
}

template <typename T>
T BaseMesh<T>::get_element() const{
    return top;
}

template <typename T>
const BaseMesh<T>* BaseMesh<T>::get_next(){
    return next;
}

template <typename T>
void BaseMesh<T>::show() const{
    //BaseMesh<T>* next_BaseMesh = get_next();
    cout << "To be fixed..." << endl;
    //top.show();
    /*do{
        top.show();
        this = next;
    }
    while(this != nullptr);*/
}

#endif
