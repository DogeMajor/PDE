#ifndef PDE_H
#define PDE_H
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <random>
#include <chrono>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Core"
#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"
#include "mesh.h"

using namespace std;
using namespace Eigen;

typedef double (* Function)(VectorXd x);
typedef Eigen::Triplet<double, int> Tri;

typedef std::mt19937 Engine;
typedef std::uniform_real_distribution<> Distribution;

class Generator {

private:
	Engine eng;
	Distribution dist;

public:
	Generator() {}
	Generator(double min, double max, int seed) : dist(min, max){
		eng.seed(seed);
	}
	Distribution::result_type gen() { return dist(eng); }
};



class Randomizer {
public:
	Randomizer(int seed) { generator = Generator(0, 1, seed); }
	~Randomizer() {}
	double prob() { return generator.gen(); }
	int random_int(int max) { return int(max*generator.gen()); }
	vector<double> randomize_items(vector<double> x);
	vector<double> get_convex_coeffs(int dim);

private:
	Generator generator;
};

vector<double> Randomizer::get_convex_coeffs(int dim) {
	vector <double> coeffs(dim);
	double sum = 0.0;
	for (int i = 0; i < dim - 1; i++) {
		coeffs[i] = (1 - sum)*prob();
		sum += coeffs[i];
	}
	coeffs[dim - 1] = 1 - sum;
	return randomize_items(coeffs);
}

vector<double> Randomizer::randomize_items(vector<double> x){
	vector<double> items;
	int index;
	int sz = x.size();
	for (int i = 0; i < sz; i++) {
		index = random_int(sz - i);
		items.push_back(x[index]);
		x.erase(x.begin() + index);
	}
	return items;
}


class Seeder {
private:
	chrono::high_resolution_clock::time_point start_time;
public:
	Seeder() {
		start_time = chrono::high_resolution_clock::now();
	}

	int get_nanoseconds() const{
		chrono::high_resolution_clock::time_point current_time =
			chrono::high_resolution_clock::now();
		return chrono::duration_cast<chrono::nanoseconds>(current_time - start_time).count();
	}
	int get_milliseconds() const {
		chrono::high_resolution_clock::time_point current_time =
			chrono::high_resolution_clock::now();
		return chrono::duration_cast<chrono::milliseconds>(current_time - start_time).count();
	}
	int get_seconds() const {
		chrono::high_resolution_clock::time_point current_time =
			chrono::high_resolution_clock::now();
		return chrono::duration_cast<chrono::seconds>(current_time - start_time).count();
	}
};


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
	void set_seeder(Seeder s) { seeder = s; }
	const BilinearFunction & get_bilinear_func() const { return A_kernel; }
    const double A(Element<Dim, Dim+1,T> el, SimplexFunction<T> a, SimplexFunction<T> b) const;//Integrates A_kernel*a*b over Element simplex
    double f(Element<Dim, Dim+1, T> &el, SimplexFunction<T> a) const;//Integrates f_kernel*a over Element simplex
	T get_random_location(Element<Dim, Dim + 1, T> &el, Randomizer &randomizer) const;
	double f_monte_carlo(Element<Dim, Dim + 1, T> &el, SimplexFunction<T> a, int n=20) const;//Integrates f_kernel*a over Element simplex

	double b(Element<Dim, Dim + 1, T> &el, SimplexFunction<T> a, SimplexFunction<T> b, BoundaryConditions<T> boundaries) const;//Surface integral of phi_i * grad(boundary_fn) on element's edge
	VectorXd to_VectorXd(T &location) const;
	PDE<Dim, T>& operator=(const PDE<Dim, T> &p);

 private:
     BilinearFunction A_kernel;
     Function f_kernel;
	 VolumeCalculator<Dim, T> volume_calculator;
	 Seeder seeder;

};

template <int Dim, typename T>
PDE<Dim, T>::PDE() {
	volume_calculator = VolumeCalculator<Dim, T>();
	seeder = Seeder();
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
double PDE<Dim, T>::f_monte_carlo(Element<Dim, Dim + 1, T> &el, SimplexFunction<T> a, int n) const {
	T loc;
	double sum = 0;
	double var = 0;
	T avg = el.get_avg_location();//OBS! If randomizer works poorly or n is small then the actual avg location of the generated locs is different
	Randomizer randomizer(seeder.get_nanoseconds());
	for (int i = 0; i <n; i++) {
		loc = get_random_location(el, randomizer);
		var += dist_squared<Dim, T>(avg,loc);
		sum = sum + f_kernel(to_VectorXd(loc))*a(loc);
	}
	var = (1 / double(n))*var;
	return sum*el.get_volume()*(1/double(n));
}

template <int Dim, typename T>//Chooses one point in the middle of simplex, returns f(P_mid)*a(P_mid)*el.volume()
T PDE<Dim, T>::get_random_location(Element<Dim, Dim + 1, T> &el, Randomizer &randomizer) const{
	T loc;
	int sz = el.nodes_size();
	vector<double> prob(sz);
	prob = randomizer.get_convex_coeffs(sz);
	loc = prob[0] * el[0].get_location();
	for (int i = 1; i < el.nodes_size(); i++) {
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
