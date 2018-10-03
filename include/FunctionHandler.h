#ifndef FUNCTIONHANDLER_H
#define FUNCTIONHANDLER_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include "node.h"
#include "element.h"
#include <math.h>
#include <functional>

using namespace std;
using namespace Eigen;

template <typename Type>
struct Function{
    Type coeff;
    virtual double operator(Type coords) = 0;
};

struct SimplexFunction{
    VectorXd coeff;
    double operator()(VectorXd coords){return coeff.head(coeff.rows()).dot(coords) + coeff.tail(1);}
};

template <int Dim, int N, typename T>//Same as for the Element template
class FunctionHandler{
/*This class creates (pyramid) functions on demand, it does not store anything*/
public:
    FunctionHandler(){}
    ~FunctionHandler(){}
    SimplexFunction build_function(Element<int Dim, int N, typename T> &el, int node_no);
    std::function<double(T)> wrap_function(SimplexFunction f);

};


template <int Dim, int N, typename T>
SimplexFunction FunctionHandler<Dim, N, T>::build_function(Element<int Dim, int N, typename T> &el, int node_no){
    VectorXd coeffs(Dim+1);
    Node<Dim, T> node =  el[int node_no];

}
