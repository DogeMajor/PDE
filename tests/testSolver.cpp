#include "../include/point.h"
#include "../include/Function.h"
#include "../include/Vertex.h"
#include "../include/element.h"
#include "../include/ElementFactory.h"
#include "../include/PDE.h"
#include "../include/Solver.h"
#include "../include/HelpfulTools.h"
#include "../include/TestingTools.h"
#include <math.h>


using namespace std;

double PI = 3.1415926535;

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

double f_kern_const(VectorXd coords) {//max val should be roughly 0.07
	return 1.0;
}

double analytic_sol(VectorXd coords) {//Particle in a N-Dim box...
	double result = (1/(2*PI*PI))*sin(coords[0] * PI);
	for (int i = 1; i < coords.rows(); i++) {
		result *= sin(coords[i] * PI);
	}
	return result;
}

//N-dim box's boundary
/*bool bound_cond(VectorXd coords) {
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
	//if (coords[1] == 1.0) { return 1; }
	return 0;
}

VectorXd bound_normal(VectorXd coords) {
	int sz = coords.size();
	VectorXd result = VectorXd::Zero(sz);
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) {
			result(i) = (coords[i] == 0.0)? -1: 1;
			return result;
		}
	}
	return result;
}

bool point_bound_cond(Point<2, double> coords) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) { return true; }
	}
	return false;
}

bool point_bound_is_inside(Point<2, double> coords) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] <= 0.0) || (coords[i] >= 1.0)) { return false; }
	}
	return true;
}

double point_bound_val(Point<2, double> coords) {
	//if (coords[1] == 1.0) { return 1; }
	return 0;
}

Point<2, double> point_bound_normal(Point<2, double> coords) {//Not needed!
	vector<double> normal = { 0.0,0.0 };
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) {
			normal[i] = (coords[i] == 0.0) ? -1 : 1;
			return Point<2, double>(normal);
		}
	}
	return Point<2, double>(normal);
}


*/

double error_norm(MatrixXd sol_values) {
	double norm = 0;
	int sz = sol_values.cols() - 1;
	VectorXd loc;
	for (int i = 0; i < sol_values.rows(); i++) {
		loc = sol_values.row(i).head(sz);
		norm = norm + pow(sol_values(i, sz) - analytic_sol(loc), 2);
	}
	return norm * (1 / double(sol_values.rows()));
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
	vector<Vertex<2, VectorXd> *> nodes(3, nullptr);
	location << 0.0, 0.0;
	nodes[0] = new Vertex<2, VectorXd>(location);
	location << 1.0, 0.0;
	nodes[1] = new Vertex<2, VectorXd>(location);
	location << 1.0, 1.0;
	nodes[2] = new Vertex<2, VectorXd>(location);
	vector<SimplexFunction <VectorXd> > funcs(3);
	VectorXd coeffs(3);
	coeffs << -1, 0, 1;
	funcs[0].coeff = coeffs;
	coeffs << 1, -1, 0;
	funcs[1].coeff = coeffs;
	coeffs << 0, 1, 0;
	funcs[2].coeff = coeffs;
	Element <2, 3, VectorXd> element(nodes, funcs);

	vector<Vertex<2, VectorXd> *> nodes2;
	nodes2.push_back(nodes[0]);
	nodes2.push_back(nodes[2]);
	location << 0, 1;
	nodes2.push_back(new Vertex<2, VectorXd>(location));
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
	srand(time(NULL));
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
	Vertex<2, Point <2, double> > n_1(point1);
	Vertex<2, Point <2, double> > n_2(point2);
	Vertex<2, Point <2, double> > n_3(point3);
	vector<Vertex <2, Point <2, double>  > * > node_vec(3, nullptr);
	node_vec[0] = new Vertex<2, Point <2, double> >(point1);
	node_vec[1] = new Vertex<2, Point <2, double> >(point2);
	node_vec[2] = new Vertex<2, Point <2, double> >(point3);
	ElementFactory<2, 3, Point <2, double> > factory;
	Element<2, 3, Point <2, double> > el1 = factory.build(node_vec);
	vector<Vertex <2, Point <2, double> > *> node_vec2;//(3, nullptr);
	node_vec2.push_back(node_vec[0]);
	node_vec2.push_back(node_vec[2]);
	node_vec2.push_back(new Vertex<2, Point <2, double> >(point4));
	Element<2, 3, Point <2, double> > el2 = factory.build(node_vec2);
	Mesh<2, 3, Point <2, double> > mesh(el1);
	mesh.push(el2);
	//mesh.reset_indices(boundaries);
	//mesh.get_top().show();
	//Element<2, 3, Point <2, double> > element2 = factory.build(node_vec);//Too lazy to find out what the funcs would be...

	Mesh<2, 3, Point <2, double> >* mesh_ptr;
	mesh_ptr = &mesh;

	BilinearFunction bl_fn;
	bl_fn.mat = MatrixXd::Identity(2, 2);
	PDE<2, Point <2, double> > pde(bl_fn, f_kern_sin);
	PDE<2, Point <2, double> > pde2(bl_fn, f_kern_const);
	BoundaryConditions<Point <2, double> > boundaries = { point_bound_cond, point_bound_is_inside, point_bound_val, point_bound_normal, 0.0000001 };

	Solver<2, Point <2, double> > solver(pde, mesh_ptr, boundaries);
	MatrixXd STIFFNESS_MAT(4, 4);
	STIFFNESS_MAT << 1, 0, -.5, -.5, 0, 1, -.5, -0.5, -.5, -.5, 1, 0, -.5, -.5, 0, 1;

	SECTION("Test solving the pde2") {
		Solver<2, Point <2, double>>  solver2(pde2, mesh_ptr, boundaries);
		REQUIRE(solver.get_stiffness_matrix(3) == STIFFNESS_MAT);
		solver2.refine();
		solver2.refine_and_adjust_to_boundary();
		solver2.refine();

		VectorXd refined_sol2 = solver2.solve();
		cout << refined_sol2 << endl;
		cout << refined_sol2.maxCoeff() << endl;
		cout << endl;
		REQUIRE(0.070 < refined_sol2.maxCoeff());
		REQUIRE(refined_sol2.maxCoeff() < 0.075);
		REQUIRE(refined_sol2.minCoeff() == 0);
		MatrixXd ref_values = solver2.get_solution_values(refined_sol2);
		cout << ref_values << endl;
	}
	
	/*SECTION("Getting sparse stiffness matrix should succeed") {
		map<array<int, 2>, double> sparse_map = solver.get_sparse_stiffness_map();
		SparseMatrix<double> test = solver.get_sparse_stiffness_matrix(3);
		MatrixXd to_dense = test.toDense();
		cout << to_dense << endl;
		cout << solver.get_stiffness_matrix(3) << endl;
		REQUIRE(to_dense == STIFFNESS_MAT);
	}

	SECTION("Getting sparse stiffness map should succeed") {
		int n = mesh.get_max_outer_index();
		map<array<int, 2>, double> s_map = solver.get_sparse_stiffness_map();
		show_map(s_map);
		REQUIRE(s_map[{0, 0}] == STIFFNESS_MAT(0, 0));
		REQUIRE(s_map[{1, 2}] == STIFFNESS_MAT(1, 2));
	}

	SECTION("Getting sparse inner stiffness matrix should succeed") {
		solver.refine();
		int max_outer = mesh.get_max_outer_index();
		SparseMatrix<double> stiff_base = solver.get_sparse_stiffness_matrix(max_outer);
		MatrixXd hehehe = stiff_base.toDense();
		cout << hehehe << endl;
		int sz = mesh.get_max_inner_index() + 1;
		cout << sz << endl;
		cout << stiff_base.coeff(0,1) << endl;
		SparseMatrix<double> inner_stiff = solver.get_sparse_inner_stiffness_matrix(solver.get_sparse_stiffness_map());
		REQUIRE(inner_stiff.coeff(0, 0) == 4);
		REQUIRE(inner_stiff.cols() == 1);
		REQUIRE(inner_stiff.rows() == 1);
	}

	SECTION("Getting sparse boundary matrix should succeed") {
		solver.refine();
		SparseMatrix<double> b_mat = solver.get_sparse_boundary_matrix(solver.get_sparse_stiffness_map());
		cout << b_mat.toDense() << endl;
		//cout << "ehehe!" << endl;
		//cout << solver.get_boundary_matrix(solver.get_stiffness_matrix(8));
	}

	SECTION("Test get_stiffness matrix(MatrixXd)") {
		MatrixXd stiffness_mat = solver.get_stiffness_matrix(3);
		REQUIRE(stiffness_mat == STIFFNESS_MAT);
		solver.refine();
		cout << solver.get_stiffness_matrix(8) << endl;
	}

	SECTION("Solving the PDE should succeed") {
		//solver.refine();
		//solver.refine();
		Seeder timer = Seeder();
		int nanosecs = 0;
		solver.refine();
		solver.refine();
		solver.refine();
		nanosecs = timer.get_nanoseconds();
		cout << "Duration in millisecs for refine(): " << timer.get_milliseconds() << endl;

		VectorXd solution = solver.solve();
		cout << "Duration in milli secs for solve(): " << timer.get_milliseconds()  << endl;

		cout << solution << endl;
		REQUIRE(solution.size() == mesh.get_max_outer_index() + 1);
		MatrixXd calc_values = solver.get_solution_values(solution);
		cout << calc_values << endl;
		cout << "Max:" << solution.maxCoeff() << " Min: " << solution.minCoeff() << endl;
		double error_squared = error_norm(calc_values);
		cout << error_squared / calc_values.rows() << endl;
		double avg = solution.mean();
		cout << "Avg. relative error norm in" << sqrt(error_squared);
		REQUIRE(sqrt(error_squared) < 0.005);
		
	}

	SECTION("Test get_boundary_coeffs") {
		solver.refine();
		VectorXd boundary_should_be(8);
		boundary_should_be << 0, 0, 0, 0, 0, 0, 0, 0;
		REQUIRE(boundary_should_be == solver.get_boundary_coeffs());
	}

	SECTION("Refining Mesh should improve the values") {
		//solver.refine();
		//solver.refine();
		solver.refine();
		solver.refine();
		VectorXd sol = solver.solve();
		cout << "Max:" << sol.maxCoeff() << " Min: " << sol.minCoeff() << endl;
		VectorXd max_r(2);
		max_r << 0.5, 0.5;
		cout << analytic_sol(max_r) << endl;
		MatrixXd values = solver.get_solution_values(sol);
		cout << "values" << endl;
		cout << values << endl;
		//mesh.get_top().show();
	}*/
	

}