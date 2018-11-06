#ifndef PDE_H
#define PDE_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Core"
#include "mesh.h"

#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"
//#include "Function.h"

using namespace std;
using namespace Eigen;

typedef double (* Function)(VectorXd x);
typedef Eigen::Triplet<double> Tri;


template <int Dim, typename T>
class PDE{

public:
	PDE();
    PDE(BilinearFunction bl_fn, Function fn){A_kernel = bl_fn; f_kernel = fn;}
	const BilinearFunction & get_bilinear_func() const { return A_kernel; }
    const double A(Element<Dim, Dim+1,T> el, SimplexFunction<T> a, SimplexFunction<T> b) const;//Integrates A_kernel*a*b over Element simplex
    double f(Element<Dim, Dim+1, T> &el, SimplexFunction<T> a) const;//Integrates f_kernel*a over Element simplex
	double b(Element<Dim, Dim + 1, T> &el, SimplexFunction<T> a, SimplexFunction<T> b, BoundaryConditions<T> boundaries) const;//Surface integral of phi_i * grad(boundary_fn) on element's edge
	VectorXd to_VectorXd(T &location) const;
	PDE<Dim, T>& operator=(const PDE<Dim, T> &p);
	//void show() const;

 private:
     BilinearFunction A_kernel;
     Function f_kernel;
     //vector <double> boundary_conds;

};

template <int Dim, typename T>
PDE<Dim, T>::PDE() {}

template <int Dim, typename T>
const double PDE<Dim, T>::A(Element<Dim, Dim+1, T> el, SimplexFunction<T> a, SimplexFunction<T> b) const{
    VectorXd coords = VectorXd::Zero(Dim);
	return A_kernel(a.gradient(), b.gradient())*el.get_volume();
}

template <int Dim, typename T>//Chooses one point in the middle of simplex, returns f(P_mid)*a(P_mid)*el.volume()
VectorXd PDE<Dim, T>::to_VectorXd(T &location) const {
	VectorXd x = VectorXd::Zero(Dim);
	for (int i = 0; i < Dim; i++) {
		x(i) = location[i];
	}
	return x;
}

template <int Dim, typename T>//Chooses one point in the middle of simplex, returns f(P_mid)*a(P_mid)*el.volume()
double PDE<Dim, T>::f(Element<Dim, Dim+1, T> &el, SimplexFunction<T> a) const{
    T avg_location = el[0].get_location();
    for(int i=1; i<Dim+1; i++){
        avg_location = avg_location + el[i].get_location();
    }
    avg_location = (1.0/double(Dim+1))*avg_location;
    return f_kernel(to_VectorXd(avg_location))*a(avg_location)*el.get_volume();
}

template <int Dim, typename T>//Assumption: Element has nodes at boundary!!
double PDE<Dim, T>::b(Element<Dim, Dim + 1, T> &el, SimplexFunction<T> a, SimplexFunction<T> b, BoundaryConditions<T> boundaries) const {
	cout << "To be coded" << endl;
	vector <Node<Dim, T>* > old_nodes = el.get_nodes();
	vector <Node<Dim, T>* > surface_nodes;
	for (int i = 0; i < old_nodes.size(); i++) {
		if (boundaries.cond(old_nodes[i]->get_location())) {
			surface_nodes.push_back(old_nodes[i]);
		}
	}

	return 0;
}

template <int Dim, typename T>
PDE<Dim, T>& PDE<Dim, T>::operator=(const PDE<Dim, T> &p) {
	A_kernel = p.A_kernel;
	f_kernel = p.f_kernel;
	return *this;
}



#endif
