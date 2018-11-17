#include "../include/Point.h"
#include "../include/Vertex.h"
#include "../include/Element.h"
#include "../include/BoundaryElementAdder.h"
#include "../include/HelpfulTools.h"
#include <math.h>

using namespace std;
using namespace Eigen;

double PI = 3.14159;
double f_kern(VectorXd coords) {
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
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) { return true; }
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


TEST_CASE("Test PDE") {
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


	BoundaryConditions<VectorXd> boundaries = { bound_cond, bound_is_inside, bound_val, bound_normal };


	SECTION("Test constructing ElementAdder") {
		BoundaryElementAdder <2, 3, VectorXd> el_adder2(boundaries);
		
	}

	
}
