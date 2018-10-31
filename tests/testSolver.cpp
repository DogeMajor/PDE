#include "../include/point.h"
#include "../include/Function.h"
#include "../include/node.h"
#include "../include/element.h"
#include "../include/PDE.h"
#include "../include/Solver.h"
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
		result *= sin(coords[i]*PI);
	}
	return result;
}



#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"


/*TEST_CASE("Test Solver") {

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
	Element <2, 3, VectorXd> element2 = factory.build(nodes2);//Too lazy to find out what the funcs would be...

	
	Mesh <2, 3, VectorXd> mesh(element);
	mesh.push(element2);
	Mesh <2, 3, VectorXd>* mesh_ptr;
	mesh_ptr = &mesh;
	
	BilinearFunction bl_fn;
	bl_fn.mat = MatrixXd::Identity(2, 2);
	PDE<2, VectorXd> pde(bl_fn, f_kern);
	//Solver<2, VectorXd> solver;
	//solver.set_pde(pde);
	//solver.set_mesh(&mesh);
	Solver<2, VectorXd> solver(pde, mesh_ptr);

	SECTION("Test default constructor") {
		Solver<2, VectorXd>  new_solver;
	}

	SECTION("Test get_stiffness matrix(MatrixXd)") {
		MatrixXd stiffness_mat = solver.get_stiffness_matrix(6);
		cout << stiffness_mat << endl;
	}

}*/


TEST_CASE("Test Solver with Point -based Mesh") {

	vector <double> vec1 = { 0.0, 0.0 };
	vector <double> vec2 = { 1.0, 0.0 };
	vector <double> vec3 = { 1.0, 1.0 };
	vector <double> vec4 = { 0.0, 1.0 };
	vector <double> vec5 = { 2.0, 0.0 };
	vector <double> vec6 = { 2.0, 1.0 };
	Point <2, double> point1(vec1);
	Point <2, double> point2(vec2);
	Point <2, double> point3(vec3);
	Point <2, double> point4(vec4);
	Point <2, double> point5(vec5);
	Point <2, double> point6(vec6);
	Node <2, Point <2, double> > n_1(point1);
	Node <2, Point <2, double> > n_2(point2);
	Node <2, Point <2, double> > n_3(point3);
	vector<Node <2, Point <2, double>  > * > node_vec(3, nullptr);
	node_vec[0] = new Node<2, Point <2, double> >(point1);
	node_vec[1] = new Node<2, Point <2, double> >(point2);
	node_vec[2] = new Node<2, Point <2, double> >(point3);
	ElementFactory<2, 3, Point <2, double> > factory;
	Element<2, 3, Point <2, double> > el1 = factory.build(node_vec);
	vector<Node <2, Point <2, double> > *> node_vec2;//(3, nullptr);
	node_vec2.push_back(node_vec[0]);
	node_vec2.push_back(node_vec[2]);
	node_vec2.push_back(new Node<2, Point <2, double> >(point4));
	Element<2, 3, Point <2, double> > el2 = factory.build(node_vec2);
	Mesh<2, 3, Point <2, double> > mesh(el1);
	mesh.push(el2);
	//mesh.get_top().show();
	//Element<2, 3, Point <2, double> > element2 = factory.build(node_vec);//Too lazy to find out what the funcs would be...

	Mesh<2, 3, Point <2, double> >* mesh_ptr;
	mesh_ptr = &mesh;

	BilinearFunction bl_fn;
	bl_fn.mat = MatrixXd::Identity(2, 2);
	cout << bl_fn.mat << endl;
	PDE<2, Point <2, double> > pde(bl_fn, f_kern_sin);
	cout << pde.A(el1, el1.get_function(0), el1.get_function(1));
	//Solver<2, VectorXd> solver;
	//solver.set_pde(pde);
	//solver.set_mesh(&mesh);
	Solver<2, Point <2, double> > solver(pde, mesh_ptr);
	//solver.show();
	
	SECTION("Test default constructor") {
		Solver<2, Point <2, double>>  new_solver;
	}

	SECTION("Test get_stiffness matrix(MatrixXd)") {
		MatrixXd stiffness_mat = solver.get_stiffness_matrix(10);
		MatrixXd mat_should_be(4,4);
		mat_should_be << 1, 0, -.5, -.5, 0, 1, -.5, 0, 0, 0, 1, 0, 0, -.5, 0, 1;
		REQUIRE(stiffness_mat == mat_should_be);
	}

	SECTION("Calculating vector f should succeed") {
		VectorXd f_vec = solver.get_vector_part();
		VectorXd vec_should_be(4);
		vec_should_be << 0.25, 0.25, 0.125, 0.125;
		for (int i = 0; i < 4; i++) { REQUIRE(limit_decimals(f_vec(i), 6) == vec_should_be(i)); }
	}

	SECTION("Solving the PDE should succeed") {
		VectorXd solution = solver.solve();
		VectorXd sol_should_be(4);
		//sol_should_be << 0.185185, 0.185185, 0.092592, 0.092592;
		sol_should_be << 0.453125, 0.3125, 0.125, 0.28125;
		for (int i = 0; i < 4; i++) { REQUIRE(limit_decimals(solution(i), 6) == sol_should_be(i)); }
	}

	SECTION("Solving the PDE after refinement should succeed") {
		solver.refine();
		solver.refine();
		//solver.show();
		VectorXd solution = solver.solve();
		cout << "Showing solution" << endl;
		cout << solution << endl;
		//Element<2, 3, Point <2, double> > temp_el = solver.get_mesh().get_top_mesh_node()->data;
		//temp_el.show();
		//cout << "Show the refined mesh" << endl;
		//solver.show();
		//VectorXd sol_should_be(4);
		//sol_should_be << 0.185185, 0.185185, 0.092592, 0.092592;
		//sol_should_be << 0.453125, 0.3125, 0.125, 0.28125;
		//for (int i = 0; i < 4; i++) { REQUIRE(limit_decimals(solution(i), 6) == sol_should_be(i)); }
	}

	/*SECTION("Getting the total solution points should succeed") {
		VectorXd solution = solver.solve();
		MatrixXd sol_values = solver.get_solution_values(solution);
		cout << "values" << endl;
		cout << sol_values << endl;
		//solver.show();
	}

	SECTION("Refining Mesh should improve the values") {
		solver.refine();
		solver.refine();
		//solver.refine();
		VectorXd sol = solver.solve();
		MatrixXd values = solver.get_solution_values(sol);
		cout << "values" << endl;
		cout << values << endl;
	}*/

}