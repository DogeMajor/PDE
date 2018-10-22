#include "../include/point.h"
#include "../include/node.h"
#include "../include/element.h"
#include "../include/PDE.h"
#include "../include/Solver.h"
#include "../include/HelpfulTools.h"
#include <math.h>

using namespace std;
using namespace Eigen;

double f_kern(VectorXd coords) {
	return coords.transpose()*coords;
}

#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"


TEST_CASE("Test Solver") {

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

	vector<Node<2, VectorXd> *> nodes2;
	nodes2.push_back(nodes[0]);
	nodes2.push_back(nodes[2]);
	location << 0, 1;
	nodes2.push_back(new Node<2, VectorXd>(location));
	ElementFactory <2, 3, VectorXd> factory;
	Element <2, 3, VectorXd> element2 = factory.build(nodes2);//Too lazy to find out whta the funcs would be...

	
	Mesh <2, 3, VectorXd> mesh(element);
	mesh.push(element2);

	BilinearFunction bl_fn;
	bl_fn.mat = MatrixXd::Identity(2, 2);
	PDE<2, VectorXd>  pde(bl_fn, f_kern);

	SECTION("Test default constructor") {
		Solver<2, VectorXd>  new_solver;
	}

	SECTION("Test inner product A(.,.)") {
		REQUIRE(pde.A(element, funcs[0], funcs[0]) == 0.5);
		REQUIRE(pde.A(element, funcs[0], funcs[1]) == -0.5);
		REQUIRE(pde.A(element, funcs[0], funcs[2]) == 0.0);
		REQUIRE(pde.A(element, funcs[1], funcs[1]) == 1.0);
		REQUIRE(pde.A(element, funcs[1], funcs[2]) == -0.5);
		REQUIRE(pde.A(element, funcs[2], funcs[2]) == 0.5);
	}


	SECTION("Test inner product with f") {
		REQUIRE(limit_decimals(pde.f(element, funcs[0]), 7) == 0.0925925);
	}
}
