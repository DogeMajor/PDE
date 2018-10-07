#ifndef ELEMENT_H
#define ELEMENT_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include "node.h"
#include <vector>
#include <functional>
//#include "Function.h"
#include "FunctionHandler.h"

using namespace std;
using namespace Eigen;

class BaseElement{
public:
    virtual ~BaseElement(){}
    virtual void show() const = 0;

};

int factorial(int n){
    return (n != 0)? n*factorial(n-1) : 1;
}



//We simply want to use the already existing nodes and don't need to worry about garbage collection.
template <int Dim, int N, typename T>
class Element{

public:
    Element();
    Element(Node<Dim,T> *nod[N]);
    Element(Node<Dim,T> *nod[N], vector <SimplexFunction> funcs);
    Element(const Element &el);//copy constructor
    ~Element();
    void increase_shared_elements();
    void set_indices();
    SimplexFunction get_function(int node_no);
    Node<Dim,T> operator[](int i);
    Element<Dim,N,T>& operator=(Element &el);
    bool operator==(const Element &el) const;
    bool operator!=(const Element &el) const;
    Matrix<double, Dim, Dim> get_simplex_matrix(Element &el) const;
    double get_volume() const;
    void show() const;

private:
    Node<Dim,T>* nodes[N];
    vector <SimplexFunction> functions;

};

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(){
    for(int i=0; i<N; i++){
        nodes[i] = new Node<Dim,T>;
    }
    increase_shared_elements();
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(Node<Dim,T> *nod[N]){
    for(int i=0; i<N; i++){//If node has no shared_elements it must be a new one!
        if(nod[i]->get_shared_elements() <= 0) {nodes[i] = new Node<Dim,T>(*nod[i]);}
        else {nodes[i] = nod[i];}
    }
    increase_shared_elements();
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(Node<Dim,T> *nod[N], vector <SimplexFunction> funcs){
    for(int i=0; i<N; i++){//If node has no shared_elements it must be a new one!
        if(nod[i]->get_shared_elements() <= 0) {nodes[i] = new Node<Dim,T>(*nod[i]);}
        else {nodes[i] = nod[i];}
    }
    increase_shared_elements();
    functions = funcs;
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(const Element &el){
    for(int i=0; i<N; i++){
        nodes[i] = el.nodes[i];
    }
    increase_shared_elements();
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::~Element(){
    int shared_elements = 0;
    for(int i=0; i<N; i++){
        shared_elements = nodes[i]->get_shared_elements();
        nodes[i]->set_shared_elements(shared_elements-1);
        if(shared_elements-1 <= 0){
            delete nodes[i];
        }
    }
    cout << "Element destroyed!" << endl;
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::increase_shared_elements(){
    int shared_elements = 0;
    for(int i=0; i<N; i++){
        shared_elements = nodes[i]->get_shared_elements();
        nodes[i]->set_shared_elements(shared_elements+1);
    }
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::set_indices(){
    for(int i=0; i<N; i++){
        nodes[i]->set_index(i);
    }
}

template <int Dim, int N, typename T>
SimplexFunction Element<Dim,N,T>::get_function(int node_no){
     return functions[node_no];
}

template <int Dim, int N, typename T>
Node<Dim,T> Element<Dim,N,T>::operator[](int i){
    return *(nodes[i]);
}

template <int Dim, int N, typename T>
Element<Dim,N,T>& Element<Dim,N,T>::operator=(Element &el){
    if(*this != el){
        for(int i=0; i<N; i++){
            if(nodes[i]->get_shared_elements() <= 0) {delete nodes[i];}
        }
        for(int i=0; i<N; i++){
            nodes[i] = el.nodes[i];
        }
    }
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
Matrix<double, Dim, Dim> Element<Dim,N,T>::get_simplex_matrix(Element &el) const{
    //For a simplex, N == Dim+1
    MatrixXd simplex_mat = MatrixXd::Zero(Dim,Dim);
    for(int col=0; col<Dim; col++){
        for(int row=0; row<Dim; row++){
            simplex_mat(row, col) = el[row+1].get_location()[col]-el[row].get_location()[col];
        }
    }
    return simplex_mat;
}

template <int Dim, int N, typename T>
double Element<Dim,N,T>::get_volume() const{
    Element<Dim,N,T> temp = *this;
    MatrixXd simplex_mat = (Dim == N-1)? get_simplex_matrix(temp): MatrixXd::Zero(Dim,Dim);
    return simplex_mat.determinant()/factorial(Dim);
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::show() const{
    cout <<"#elements: " << N << endl;
    for(int i=0; i<N; i++){
        nodes[i]->show();
    }
}


template <int Dim, int N, typename T>
class ElementFactory{

public:
    ElementFactory(){}
    ~ElementFactory(){}
    Element<Dim, N, T> build(Node<Dim,T> *nod[N]);
    MatrixXd get_inv_matrix(Node<Dim,T>* nodes[]);
    SimplexFunction build_function(MatrixXd M, int node_no);
    vector <SimplexFunction> build_functions(Node<Dim,T>* nodes[]);
};

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementFactory<Dim,N,T>::build(Node<Dim,T> *nod[N]){
    vector <SimplexFunction> funcs = build_functions(nod);
    Element<Dim, N, T> el(nod, funcs);
    return el;
}

template <int Dim, int N, typename T>
MatrixXd ElementFactory<Dim,N,T>::get_inv_matrix(Node<Dim,T>* nodes[]){
    MatrixXd M(Dim+1, Dim+1);
    for(int node=0; node<Dim+1; node++){
        for(int col=0; col<Dim; col++){
            M(node,col) = nodes[node]->get_location()[col];
        }
        M(node,Dim) = 1;
    }
    return M.inverse();
}

template <int Dim, int N, typename T>
SimplexFunction ElementFactory<Dim,N,T>::build_function(MatrixXd M, int node_no){
    VectorXd unit_vec = VectorXd::Zero(Dim+1);
    unit_vec(node_no) = 1;
    VectorXd coeffs = M*unit_vec;
    SimplexFunction fn_obj;
    fn_obj.coeff = coeffs;
    return fn_obj;
}

template <int Dim, int N, typename T>
vector <SimplexFunction> ElementFactory<Dim,N,T>::build_functions(Node<Dim,T>* nodes[]){
    MatrixXd M_inv = get_inv_matrix(nodes);
    vector <SimplexFunction> functions;
    for(int i=0; i<Dim+1; i++){
        functions.push_back(build_function(M_inv, i));
    }
    return functions;
}

#endif
