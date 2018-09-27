#ifndef MESH_H
#define MESH_H
#include <iostream>
#include "node.h"
#include "element.h"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include <math.h>

using namespace std;
using namespace Eigen;


template <int Dim, int N, typename T>
class Mesh{

public:
    Mesh();
    Mesh(Element<Dim,N,T> *top);
    ~Mesh();
    Element<Dim,N,T>& get_next();
    Matrix<double, Dim, Dim> get_simplex_matrix(Element<Dim,N,T> &el) const;//For calculating the volume
    double get_element_volume(Element<Dim,N,T> &el) const;
    double weak_form_element(int i, int j) const;
    void show() const;

private:
    Element<Dim,N,T> *top;
    //Element<N,T> *next;

};

template <int Dim, int N, typename T>
Mesh<Dim,N,T>::Mesh(){
    top = new Element<Dim,N,T>;
    //next = new Element<N,T>;
}


template <int Dim, int N, typename T>
Mesh<Dim,N,T>::Mesh(Element<Dim,N,T> *top){
    top = top;
}

template <int Dim, int N, typename T>
Mesh<Dim,N,T>::~Mesh(){
    delete top;
    //delete next;
    cout << "Mesh destroyed!" << endl;
}


template <int Dim, int N, typename T>
Element<Dim,N,T>& Mesh<Dim,N,T>::get_next(){
    return *top;
}


template <int Dim, int N, typename T>
Matrix<double, Dim, Dim> Mesh<Dim,N,T>::get_simplex_matrix(Element<Dim,N,T> &el) const{
    //For a simplex, N == Dim+1
    MatrixXd simplex_mat = MatrixXd::Zero(Dim,Dim);
    for(int col=0; col<Dim; col++){
        for(int row=0; row<Dim; row++){
            simplex_mat(row, col) = el[row+1].get_location()(col)-el[row].get_location()(col);
        }
    }
    return simplex_mat;
}

template <int Dim, int N, typename T>
double Mesh<Dim,N,T>::get_element_volume(Element<Dim,N,T> &el) const{
    MatrixXd simplex_mat = get_simplex_matrix(el);
    return simplex_mat.determinant();
}

template <int Dim, int N, typename T>
double Mesh<Dim,N,T>::weak_form_element(int i, int j) const{
    return 0.0;
}

template <int Dim, int N, typename T>
void Mesh<Dim,N,T>::show() const{
    top->show();
}

#endif


