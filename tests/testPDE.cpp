#include "../include/point.h"
#include "../include/node.h"
#include "../include/element.h"
#include "../include/PDE.h"
#include "../include/HelpfulTools.h"
#include <math.h>

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

	VectorXd location(2);
	location << 0.0, 0.0;
	Node <2, VectorXd> node1(location);
	location << 1.0, 0.0;
	Node <2, VectorXd> node2(location);
	location << 1.0, 1.0;
	Node <2, VectorXd> node3(location);
	vector<Node<2, VectorXd> *> nodes(3, nullptr);
	location << 0.0, 0.0;
	nodes[0] = new Node<2, VectorXd>(location);
	location << 1.0, 0.0;
	nodes[1] = new Node<2, VectorXd>(location);
	location << 1.0, 1.0;
	nodes[2] = new Node<2, VectorXd>(location);
	vector<SimplexFunction <VectorXd> > funcs(3);
	VectorXd coeffs(3);
	coeffs << -1, 0, 1;
	funcs[0].coeff = coeffs;
	coeffs << 1, -1, 0;
	funcs[1].coeff = coeffs;
	coeffs << 0, 1, 0;
	funcs[2].coeff = coeffs;

	Element <2, 3, VectorXd> element(nodes, funcs);

    BilinearFunction bl_fn;
    bl_fn.mat = MatrixXd::Identity(2,2);
    PDE<2, VectorXd>  pde(bl_fn, f_kern_sin);

	BoundaryConditions<VectorXd> boundaries = {bound_cond, bound_val, bound_normal };
	
	chrono::high_resolution_clock::time_point start_t;
	Seeder seeder = { start_t };
	seeder.start_clock();
	RNGType new_gen(seeder.get_nanoseconds());
	//new_gen.seed(seeder.get_nanoseconds());

	//srand(time(NULL));
	//RNGType new_gen(time(0));
	Randomizer randomizer = Randomizer(new_gen);
	
	SECTION("Random_prob funtions should work") {
		
		double result;
		double avg = 0;
		for (int i = 0; i < 50; i++) {
			//nanoseconds_time();
			//new_gen.seed(seeder.get_nanoseconds());
			//result = random01(new_gen);
			result = randomizer.random01();
			avg += result;
			cout  << result << endl;
		}
		avg = avg * (1 / double(50));
		cout << "Avg: " << avg << endl;
		REQUIRE(abs(avg - 0.5) < 0.05);
		vector<double> coeffs = get_convex_coeffs(new_gen, 5);
		cout << "sum" << sum(coeffs) << endl;
		show_vector<vector<double> >(coeffs);
		vector<double> items = randomize_items(new_gen, coeffs);
		cout << "Rand" << endl;
		show_vector<vector<double> >(items);
	}

	/*SECTION("BoundaryConditions for 2-D box should work") {
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
		//cout << pde.f(element, funcs[0]) << endl;
		//cout << pde.f(element, funcs[1]) << endl;

		cout << pde.f(element, funcs[2]) << endl;

	}

	SECTION("Getting random location should succeed") {
		VectorXd rand_loc = pde.get_random_location(element);
		cout << rand_loc << endl;
		for (int i = 0; i > rand_loc.size(); i++) {
			REQUIRE(rand_loc[i] <= 1);
			REQUIRE(rand_loc[i] >= 0);
		}
		
	}

	SECTION("Test f_monte_carlo(.,.)") {//More accurate integration!!
		//cout << pde.f(element, funcs[0]) << endl;
		//cout << pde.f(element, funcs[1]) << endl;
		element.show();
		cout << funcs[2].coeff << endl;
		cout << "Monte carlo intergals" << endl;
		cout << pde.f_monte_carlo(element, funcs[2], 10) << endl;
		cout << pde.f_monte_carlo(element, funcs[2], 20) << endl;
		REQUIRE(pde.f_monte_carlo(element, funcs[2], 20) < 0.08);
		//cout << pde.f_monte_carlo(element, funcs[2], 30) << endl;
		cout << pde.f_monte_carlo(element, funcs[2], 50) << endl;
		cout << pde.f_monte_carlo(element, funcs[2], 100) << endl;
		cout << pde.f_monte_carlo(element, funcs[2], 600) << endl;
		//cout << pde.f_monte_carlo(element, funcs[2], 4000) << endl;
		//cout << pde.f_monte_carlo(element, funcs[2], 5000) << endl;

	}


    SECTION( "Test inner product with f" ){
        REQUIRE( limit_decimals(pde.f(element, funcs[0]),7) == 0.0925925 );
    }*/
}
