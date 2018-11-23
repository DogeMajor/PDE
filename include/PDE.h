#ifndef PDE_H
#define PDE_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Core"
#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"
#include "Mesh.h"
#include "Randomizer.h"

using namespace std;
using namespace Eigen;

typedef double (* Function)(VectorXd x);
typedef Eigen::Triplet<double, int> Tri;


double sum(vector<double> x) {
	double sum = 0.0;
	for (int i = 0; i < x.size(); i++) { sum += x[i]; }
	return sum;
}

template <int Dim, typename T>
class PDE{

public:
	PDE();
    PDE(BilinearFunction bl_fn, Function fn){A_kernel = bl_fn; f_kernel = fn;}
	void set_timer(Timer s) { timer = s; }
	const BilinearFunction & get_bilinear_func() const { return A_kernel; }
    const double A(Element<Dim, Dim+1,T> el, SimplexFunction<T> a, SimplexFunction<T> b) const;//Integrates A_kernel*a*b over Element simplex
    double f(Element<Dim, Dim+1, T> &el, SimplexFunction<T> a) const;//Integrates f_kernel*a over Element simplex
	T get_random_location(Element<Dim, Dim + 1, T> &el, Randomizer &randomizer) const;
	double f_monte_carlo(Element<Dim, Dim + 1, T> &el, SimplexFunction<T> a, int fn_index, int n=20) const;//Integrates f_kernel*a over Element simplex

	double b(Element<Dim, Dim + 1, T> &el, SimplexFunction<T> a, SimplexFunction<T> b, BoundaryConditions<T> boundaries) const;//Surface integral of phi_i * grad(boundary_fn) on element's edge
	VectorXd to_VectorXd(T &location) const;
	PDE<Dim, T>& operator=(const PDE<Dim, T> &p);

 private:
     BilinearFunction A_kernel;
     Function f_kernel;
	 VolumeCalculator<Dim, T> volume_calculator;
	 Timer timer;

};

template <int Dim, typename T>
PDE<Dim, T>::PDE() {
	volume_calculator = VolumeCalculator<Dim, T>();
	timer = Timer();
}

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

template <int Dim, typename T>
double PDE<Dim, T>::f_monte_carlo(Element<Dim, Dim + 1, T> &el, SimplexFunction<T> a, int fn_index, int n) const {
	T loc;
	double sum = 0;
	double var = 0;
	T avg = el.get_avg_location();//OBS! If randomizer works poorly or n is small then the actual avg location of the generated locs is different
	Randomizer randomizer(timer.get_nanoseconds());
	for (int i = 0; i <n; i++) {
		loc = get_random_location(el, randomizer);
		var += dist_squared<Dim, T>(avg,loc);
		sum = sum + f_kernel(to_VectorXd(loc))*a(loc);
	}
	var = (1 / double(n))*var;
	el.set_f_variation(fn_index, var);
	return sum*el.get_volume()*(1/double(n));
}

template <int Dim, typename T>//Chooses one point in the middle of simplex, returns f(P_mid)*a(P_mid)*el.volume()
T PDE<Dim, T>::get_random_location(Element<Dim, Dim + 1, T> &el, Randomizer &randomizer) const{
	T loc;
	int sz = el.vertices_size();
	vector<double> prob(sz);
	prob = randomizer.get_convex_coeffs(sz);
	loc = prob[0] * el[0].get_location();
	for (int i = 1; i < el.vertices_size(); i++) {
		loc = loc + prob[i] * el[i].get_location();
	}
	return loc;
}

template <int Dim, typename T>//Assumption: Element has nodes at the boundary!!
double PDE<Dim, T>::b(Element<Dim, Dim + 1, T> &el, SimplexFunction<T> a, SimplexFunction<T> b, BoundaryConditions<T> boundaries) const {
	cout << "To be coded" << endl;
	vector <Node<Dim, T>* > old_nodes = el.get_nodes();
	vector <Node<Dim, T>* > surface_nodes;
	T avg_loc;
	for (int i = 0; i < old_nodes.size(); i++) {
		if (boundaries.cond(old_nodes[i]->get_location())) {
			surface_nodes.push_back(old_nodes[i]);
		}
	}

	double vol = volume_calculator.get_dist_volume(old_nodes);
}

template <int Dim, typename T>
PDE<Dim, T>& PDE<Dim, T>::operator=(const PDE<Dim, T> &p) {
	A_kernel = p.A_kernel;
	f_kernel = p.f_kernel;
	return *this;
}


#endif
