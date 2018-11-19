#include "../include/Point.h"
#include "../include/Vertex.h"
#include "../include/Element.h"
#include "../include/Randomizer.h"
#include "../include/PDE.h"
#include "../include/HelpfulTools.h"
#include "../include/Function.h"
#include <math.h>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Core"
#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"

using namespace std;
using namespace Eigen;

double PI = 3.14159;
double f_kern(VectorXd coords){
    return coords.transpose()*coords;
}

double f_kern_sin(VectorXd coords) {//Particle in a N-Dim box...
	double result = sin(coords[0] * PI);
	for (int i = 1; i < coords.rows(); i++) {
		result *= sin(coords[i] * PI);
	}
	return result;
}

//N-dim box's boundary
bool bound_cond(VectorXd coords) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) {return true;}
	}
	return false;
}

bool bound_is_inside(VectorXd coords) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] <= 0.0) || (coords[i] >= 1.0)) { return false; }
	}
	return true;
}

double bound_val(VectorXd coords) {
	if (coords[1] == 1.0) { return 1; }
	return 0;
}

VectorXd bound_normal(VectorXd coords) {
	int sz = coords.size();
	VectorXd result = VectorXd::Zero(sz);
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) {
			result(i) = (coords[i] == 0.0) ? -1 : 1;
			return result;
		}
	}
	return result;
}

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test PDE" ) {
	double PRODUCT_F_PHI2 = 0.0759909;
	VectorXd location(2);
	location << 0.0, 0.0;
	Vertex<2, VectorXd> node1(location);
	location << 1.0, 0.0;
	Vertex<2, VectorXd> node2(location);
	location << 1.0, 1.0;
	Vertex<2, VectorXd> node3(location);
	vector<Vertex<2, VectorXd> *> vertices(3, nullptr);
	location << 0.0, 0.0;
	vertices[0] = new Vertex<2, VectorXd>(location);
	location << 1.0, 0.0;
	vertices[1] = new Vertex<2, VectorXd>(location);
	location << 1.0, 1.0;
	vertices[2] = new Vertex<2, VectorXd>(location);
	vector<SimplexFunction <VectorXd> > funcs(3);
	VectorXd coeffs(3);
	coeffs << -1, 0, 1;
	funcs[0].coeff = coeffs;
	coeffs << 1, -1, 0;
	funcs[1].coeff = coeffs;
	coeffs << 0, 1, 0;
	funcs[2].coeff = coeffs;

	Element <2, 3, VectorXd> element(vertices, funcs);

    BilinearFunction bl_fn;
    bl_fn.mat = MatrixXd::Identity(2,2);
    PDE<2, VectorXd>  pde(bl_fn, f_kern_sin);

	BoundaryConditions<VectorXd> boundaries = {bound_cond, bound_is_inside, bound_val, bound_normal };

	chrono::high_resolution_clock::time_point start_t;
	Seeder seeder = Seeder();

	Randomizer randomizer = Randomizer(seeder.get_nanoseconds());
	
	SECTION("Random_prob funtions should work") {
		Generator rng_gen(0, 1, seeder.get_nanoseconds());
		double result;
		double avg = 0;
		int max_iter = 30;
		for (int i = 0; i < max_iter; i++) {
			result = randomizer.prob();
			avg += result;
			cout  << result << endl;
		}
		avg = avg * (1 / double(max_iter));
		cout << "Avg: " << avg << endl;
		REQUIRE(abs(avg - 0.5) < 0.10);
		vector<double> coeffs = randomizer.get_convex_coeffs(5);
		cout << "sum" << sum(coeffs) << endl;
		REQUIRE(limit_decimals(sum(coeffs), 3) == 1.000);
		show_vector<vector<double> >(coeffs);
		vector<double> items = randomizer.randomize_items(coeffs);
		REQUIRE(abs(sum(items)-1) < 0.001);
		show_vector<vector<double> >(items);
	}

	SECTION("BoundaryConditions for 2-D box should work") {
		PDE<2, VectorXd>  new_pde(BilinearFunction bl_fn, Function f_kernel);
	}

	SECTION( "Test constructing PDE" ){
		VectorXd temp(2);
		temp << 0.5, 0.5;
		REQUIRE(boundaries.cond(temp)==false);
		REQUIRE(boundaries.val(temp) == 0.0);
		temp << 0.5, 1;
		REQUIRE(boundaries.cond(temp) == true);
		REQUIRE(boundaries.val(temp) == 1.0);
    }

	SECTION( "Test inner product A(.,.)" ){
        REQUIRE( pde.A(element, funcs[0], funcs[0]) == 0.5 );
        REQUIRE( pde.A(element, funcs[0], funcs[1]) == -0.5 );
        REQUIRE( pde.A(element, funcs[0], funcs[2]) == 0.0 );
        REQUIRE( pde.A(element, funcs[1], funcs[1]) == 1.0 );
        REQUIRE( pde.A(element, funcs[1], funcs[2]) == -0.5 );
        REQUIRE( pde.A(element, funcs[2], funcs[2]) == 0.5 );
    }

	SECTION("Test f(.,.)") {

		cout <<"result f" << pde.f(element, funcs[2]) << endl;
		REQUIRE(abs(PRODUCT_F_PHI2 - pde.f(element, funcs[2])) < 0.1);

	}

	SECTION("Getting random location should succeed") {
		vector<double> probs = randomizer.get_convex_coeffs(2);
		cout << probs[0] <<", " << probs[1] << endl;
		VectorXd rand_loc = pde.get_random_location(element, randomizer);
		cout << rand_loc << endl;
		for (int i = 0; i > rand_loc.size(); i++) {
			REQUIRE(rand_loc[i] <= 1);
			REQUIRE(rand_loc[i] >= 0);
		}
	}

	SECTION("Test f_monte_carlo(.,.)") {//More accurate integration!!

		cout << funcs[2].coeff << endl;
		cout << "Monte carlo intergals" << endl;
		cout << pde.f_monte_carlo(element, funcs[2], 2, 10) << endl;
		REQUIRE(pde.f_monte_carlo(element, funcs[2], 2, 20) < 0.09);
		cout << pde.f_monte_carlo(element, funcs[2], 2, 50) << endl;
		REQUIRE(abs(pde.f_monte_carlo(element, funcs[2], 2, 100) - PRODUCT_F_PHI2) < 0.02);
		cout << pde.f_monte_carlo(element, funcs[2], 2, 1000) << endl;
		REQUIRE(element.get_avg_f_variation() < 0.10);
		cout << "var" << element.get_avg_f_variation() << endl;

	}

    SECTION( "Test inner product with f" ){
        REQUIRE( limit_decimals(pde.f(element, funcs[0]),3) == 0.125 );
    }
}
