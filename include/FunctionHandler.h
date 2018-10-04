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
typedef double (* Function)(VectorXd x);

struct SimplexFunction{
    VectorXd coeff;
    double operator()(VectorXd coords){
        double result = dot_product(coeff,coords);
        return result + double(coeff.tail(1)[0]);
        }
    double dot_product(VectorXd a, VectorXd b){
        int min_size = min(a.size(), b.size());
        double result = 0.0;
        for(int i=0; i<min_size; i++){result += a[i]*b[i];}
        return result;
        }
    VectorXd gradient(VectorXd coords){
        int coeff_size = coeff.size();
        return coeff.head(coeff_size-1);
    }
};


struct BilinearFunction{
    MatrixXd mat;
    double operator()(VectorXd x, VectorXd y){
        double result = (x.transpose())*mat*y;
        return result;
        }
};

template <int Dim, int N, typename T>//Same as for the Element template
class FunctionGenerator{
    public:
    FunctionGenerator(){}
    ~FunctionGenerator(){}
    MatrixXd get_inv_matrix(Element<Dim, N, T> &el);
    SimplexFunction build_function(MatrixXd M, int node_no);
    //std::function<double(&SimplexFunction)> wrap_function(SimplexFunction f);

};

template <int Dim, int N, typename T>
MatrixXd FunctionGenerator<Dim, N, T>::get_inv_matrix(Element<Dim, N, T> &el){
    MatrixXd M(Dim+1, Dim+1);
    for(int node=0; node<N; node++){
        for(int col=0; col<Dim; col++){
            M(node,col) = el[node].get_location()[col];
        }
        M(node,Dim) = 1;
    }
    return M.inverse();
}

template <int Dim, int N, typename T>
SimplexFunction FunctionGenerator<Dim, N, T>::build_function(MatrixXd M, int node_no){
    VectorXd unit_vec = VectorXd::Zero(Dim+1);
    unit_vec(node_no) = 1;
    VectorXd coeffs = M*unit_vec;
    SimplexFunction fn_obj;
    fn_obj.coeff = coeffs;
    return fn_obj;
}

/*template <int Dim, int N, typename T>
std::function<double(&SimplexFunction)>  FunctionGenerator<Dim, N, T>::wrap_function(SimplexFunction f){
    //MyValue&)> fifth = &MyValue::fifth;
    std::function<double(&SimplexFunction)> func = &SimplexFunction::();
    return func;
}*/

template <int Dim, int N, typename T>//Same as for the Element template
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
    //std::function<double(&SimplexFunction)> wrap_function(SimplexFunction f);

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

#endif
