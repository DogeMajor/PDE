#ifndef FUNCTIONHANDLER_H
#define FUNCTIONHANDLER_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include "node.h"
#include "element.h"
#include "Function.h"
#include <math.h>
#include <functional>
#include <vector>

using namespace std;
using namespace Eigen;
typedef double (* Function)(VectorXd x);


/*
template <int Dim, int N, typename T>//Same as for the Element template hmm... this starts resembling PDE...
class FunctionAnalyzer{
    public:
    FunctionAnalyzer(BilinearFunction &bl_fn, Function fn){A = bl_fn; f = fn;}
    FunctionAnalyzer(){}
    ~FunctionAnalyzer(){}
    double bilinear_form(VectorXd a, VectorXd b);
    double sobolev_dot_product(Element<Dim, N, T> &el, SimplexFunction a, SimplexFunction b);
    double sobolev_f(Element<Dim, N, T> &el, SimplexFunction a);
    private:
    BilinearFunction A; //A(u,v) = f(v) for all v in V, dim(V) = #nodes
    Function f;

};

template <int Dim, int N, typename T>
double FunctionAnalyzer<Dim, N, T>::bilinear_form(VectorXd a, VectorXd b){
    return A(a, b);
}

template <int Dim, int N, typename T>
double FunctionAnalyzer<Dim, N, T>::sobolev_dot_product(Element<Dim, N, T> &el, SimplexFunction a, SimplexFunction b){
    VectorXd coords = VectorXd::Zero(Dim);
    return A(a.gradient(coords), b.gradient(coords))*el.get_volume();
}

template <int Dim, int N, typename T>//Chooses one point in the middle of simplex, returns f(P_mid)*a(P_mid)*el.volume()
double FunctionAnalyzer<Dim, N, T>::sobolev_f(Element<Dim, N, T> &el, SimplexFunction a){
    VectorXd avg_location = VectorXd::Zero(Dim);
    for(int i=0; i<N; i++){
        avg_location += el[i].get_location();
    }
    avg_location = (1.0/double(N))*avg_location;
    return f(avg_location)*a(avg_location)*el.get_volume();
}

*/
#endif
