#include "../include/Point.h"
#include "../include/Function.h"
#include "../include/Vertex.h"
#include "../include/Element.h"
#include "../include/ElementFactory.h"
#include "../include/PDE.h"
#include "../include/DAO.H"
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

double f_kern_const(VectorXd coords) {//max val of the solution should be roughly 0.07
	return 1.0;
}

double f_kern_zero(VectorXd coords) {
	return 0.0;
}

double analytic_sol(VectorXd coords) {//Particle in a N-Dim box...
	int dim = coords.size();
	double result = (1/(dim*PI*PI))*sin(coords[0] * PI);
	for (int i = 1; i < coords.rows(); i++) {
		result *= sin(coords[i] * PI);
	}
	return result;
}


double circle_sol(VectorXd coords) {
	if (coords[0] == 0 && coords[1] == 0) { return 0; }
	double squared_r = coords[0] * coords[0] + coords[1] * coords[1];
	double angle = acos(coords[0] / sqrt(squared_r));
	return pow(squared_r, 5)*cos(10 * angle);
}

double c_error_norm(MatrixXd sol_values) {
	double norm = 0;
	int sz = sol_values.cols() - 1;
	VectorXd loc;
	for (int i = 0; i < sol_values.rows(); i++) {
		loc = sol_values.row(i).head(sz);
		norm = norm + pow(sol_values(i, sz) - circle_sol(loc), 2);
	}
	cout << "error sum" << norm << endl;
	return norm * (1 / double(sol_values.rows()));
}


double error_norm(MatrixXd sol_values) {
	double norm = 0;
	int sz = sol_values.cols() - 1;
	VectorXd loc;
	for (int i = 0; i < sol_values.rows(); i++) {
		loc = sol_values.row(i).head(sz);
		norm = norm + pow(sol_values(i, sz) - analytic_sol(loc), 2);
	}
	cout << "error sum" << norm << endl;
	return norm * (1 / double(sol_values.rows()));
}

#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"


/*TEST_CASE("Test Solver with Point -based Mesh") {
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
		solver2.refine();
		solver2.refine();
		solver2.refine();
		VectorXd refined_sol2 = solver2.solve();

		REQUIRE(0.070 < refined_sol2.maxCoeff());
		REQUIRE(refined_sol2.maxCoeff() < 0.075);
		REQUIRE(refined_sol2.minCoeff() == 0);
		MatrixXd ref_values = solver2.get_solution_values(refined_sol2);
		cout << ref_values << endl;
	}
	
	SECTION("Getting sparse stiffness matrix should succeed") {
		int tot_sz = mesh_ptr->get_max_outer_index() + 1;
		SparseMatrix<double> test = solver.get_sparse_stiffness_matrix(tot_sz);
		MatrixXd to_dense = test.toDense();
		REQUIRE(to_dense == STIFFNESS_MAT);
	}

	SECTION("Getting sparse stiffness matrix should succeed") {
		solver.refine();
		int max_outer = mesh.get_max_outer_index()+1;
		SparseMatrix<double> stiff_base = solver.get_sparse_stiffness_matrix(max_outer);
		MatrixXd tot_stiffness = stiff_base.toDense();
		REQUIRE(tot_stiffness.coeff(0, 0) == 4);
		REQUIRE(tot_stiffness.cols() == 9);
		REQUIRE(tot_stiffness.rows() == 9);
	}

	SECTION("Solving the PDE should succeed") {
		//solver.refine();
		//solver.refine();
		Timer timer = Timer();
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
	}

	SECTION("Solving PDE in a unit circle domain should succeed") {
		VectorXd temp_loc(2);
		temp_loc << 1, 0;
		//In case we need to compare vertices later on with the ones in Mesh!
		Vertex<2, VectorXd > c_v1(temp_loc);
		temp_loc << 0, 1;
		Vertex<2, VectorXd > c_v2(temp_loc);
		temp_loc << -1, 0;
		Vertex<2, VectorXd > c_v3(temp_loc);
		temp_loc << 0, -1;
		Vertex<2, VectorXd > c_v4(temp_loc);


		vector<Vertex <2, VectorXd> * > c_vertices1(3, nullptr);
		temp_loc << 1, 0;
		c_vertices1[0] = new Vertex <2, VectorXd>(temp_loc);
		temp_loc << 0, 1;
		c_vertices1[1] = new Vertex <2, VectorXd>(temp_loc);
		temp_loc << -1, 0;
		c_vertices1[2] = new Vertex <2, VectorXd>(temp_loc);
		ElementFactory<2, 3, VectorXd> factory;
		Element<2, 3, VectorXd> c_el1 = factory.build(c_vertices1);
		vector<Vertex <2, VectorXd> *> c_vertices2;//(3, nullptr);
		c_vertices2.push_back(c_vertices1[0]);
		c_vertices2.push_back(c_vertices1[2]);
		temp_loc << 0, -1;
		c_vertices2.push_back(new Vertex<2,VectorXd >(temp_loc));


		Element<2, 3, VectorXd> c_el2 = factory.build(c_vertices2);
		Mesh<2, 3, VectorXd> c_mesh(c_el1);
		c_mesh.push(c_el2);
		//mesh.reset_indices(boundaries);
		//mesh.get_top().show();
		//Element<2, 3, Point <2, double> > element2 = factory.build(node_vec);//Too lazy to find out what the funcs would be...

		Mesh<2, 3, VectorXd>* c_mesh_ptr;
		c_mesh_ptr = &c_mesh;
		
		BilinearFunction c_bl_fn;
		c_bl_fn.mat = MatrixXd::Identity(2, 2);
		PDE<2, VectorXd> c_pde(c_bl_fn, f_kern_zero);
		
		BoundaryConditions<VectorXd> c_boundaries = { c_cond, c_is_inside, c_val, c_normal, 0.000001 };

		Solver<2, VectorXd> c_solver(c_pde, c_mesh_ptr, c_boundaries);

		SECTION("Solving pde in unit circle should succeed") {
			c_solver.refine();
			c_solver.refine();
			c_solver.refine();
			c_solver.refine();
			//c_solver.refine();
			//c_mesh.show();
			cout << c_mesh.get_max_inner_index() << endl;
			cout << c_mesh.get_max_outer_index() << endl;

			//cout << c_solver.get_sparse_inner_stiffness_matrix(c_stiffness_map).toDense() << endl;
			c_solver.save_grid("refined_grid.txt");
			VectorXd temp(2);
			temp << 1, 0;
			cout << "Cos(c_v1)" << c_boundaries.val(temp);

			VectorXd c_sol = c_solver.solve();
			cout << c_sol << endl;
			cout << c_sol.maxCoeff() << endl;
			cout << endl;
			
			MatrixXd c_values = c_solver.get_solution_values(c_sol);
			cout << c_values << endl;
			c_solver.save_grid("c_mesh_grid.txt");
			double c_error_avg = c_error_norm(c_values);
			REQUIRE(c_error_avg < 0.05);
			cout << "Average error" << c_error_avg;
			cout << "circle_sol at 0,1" << circle_sol(c_v2.get_location());
			cout << "circle_sol at 0,0.5" << circle_sol(0.5*c_v2.get_location());
		}
	}
}*/

TEST_CASE("solving 3-D poisson should succeed") {
	VectorXd cube_center(3);
	cube_center << 0.5, 0.5, 0.5;
	VectorXd cube_lengths(3);
	cube_lengths << 1, 1, 1;
	VectorXd p_loc(3);

	vector<Vertex<3, VectorXd> *> p_vertices(4, nullptr);
	p_loc << 0, 0, 0;
	p_vertices[0] = new Vertex<3, VectorXd>(p_loc);
	p_loc << 1, 0.5, 0.5;
	p_vertices[1] = new Vertex<3, VectorXd>(p_loc);
	p_loc << 0.5, 1, 0.5;
	p_vertices[2] = new Vertex<3, VectorXd>(p_loc);
	p_loc << 0, 0, 1;
	p_vertices[3] = new Vertex<3, VectorXd>(p_loc);

	ElementFactory<3, 4, VectorXd> p_factory;
	Element<3, 4, VectorXd> p_element = p_factory.build(p_vertices);
	p_element.set_indices(-1);
	p_element.set_index_maps();

	BoundaryConditions<VectorXd> p_boundaries = { bound_cond, bound_is_inside, bound_val, bound_normal, 0.000001 };
	ElementDivider <3, 4, VectorXd> p_divider(p_boundaries);

	MeshDAO<3, 4, VectorXd> p_mesh_dao;
	SolverDAO<3, 4, VectorXd> p_solver_dao;
	Mesh<3, 4, VectorXd> p_mesh(p_element);

	map< array<int, 2>, int> P_EDGES_MAP;
	int A, B;
	int I = 0;
	for (int i = 0; i < 4; i++) {
		A = p_element[i].get_index();
		for (int j = i + 1; j < 4; j++) {
			B = p_element[j].get_index();
			P_EDGES_MAP.insert(pair< array<int, 2>, int>({ min(A,B), max(A,B) }, 1));
			I++;
		}
	}
	map<array<int, 2>, int> p_empty_edges;
	Mesh<3, 4, VectorXd> p_empty_mesh;
	Mesh <3, 4, VectorXd>* p_mesh_ptr = &p_empty_mesh;

	BilinearFunction p_bl_fn;
	p_bl_fn.mat = MatrixXd::Identity(3, 3);
	PDE<3, VectorXd> p_pde(p_bl_fn, f_kern_const);
	Solver<3, VectorXd> p_solver(p_pde, p_mesh_ptr, p_boundaries);
	Timer timer = Timer();

	SECTION("Filling mesh with simplices covering unit box with origo (0.5,0.5,0.5) should succeed") {
		p_solver.fill_mesh_covering_box(cube_center, cube_lengths);
		p_solver.refine();
		p_solver.refine();
		int MAX_INNER = p_mesh_ptr->get_max_inner_index();
		int MAX_OUTER = p_mesh_ptr->get_max_outer_index();
		int MAX_NODE = p_mesh_ptr->how_many_nodes();

		/*SECTION("Get f_vec should succeed") {
			VectorXd p_f_vec = p_solver.get_f_vec(MAX_INNER);
			//cout << p_f_vec << endl;
			REQUIRE(p_f_vec.maxCoeff() < 0.022);
			REQUIRE(p_f_vec.minCoeff() > 0.008);

		}

		SECTION("Getting partial f_vec should succeed") {
			VectorXd f_vec1 = p_solver.get_f_vec_part(0, MAX_NODE/2);
			//cout << f_vec1 << endl;
			REQUIRE(f_vec1.maxCoeff() < 0.022);
			REQUIRE(f_vec1.minCoeff() == 0.0);
			VectorXd f_vec2 = p_solver.get_f_vec_part(MAX_NODE / 2 + 1, MAX_NODE);
			//cout << f_vec2 << endl;
			REQUIRE(f_vec2.maxCoeff() < 0.022);
			REQUIRE(f_vec2.minCoeff() == 0.0);
			VectorXd f_vec_tot = f_vec1 + f_vec2;
			VectorXd difference = f_vec_tot - p_solver.get_f_vec(MAX_INNER);
			//cout << difference << endl;
			REQUIRE(difference.maxCoeff() < 0.002);

		}

		SECTION("Get f_vec2 should succeed") {
			int start_time = timer.get_milliseconds();
			VectorXd p_f_vec_async = p_solver.get_f_vec_async(4);
			int t_async = timer.get_milliseconds() - start_time;
			cout << "Calculating f_vec asynchronously takes this much time in ms: " << t_async << endl;
			cout << p_f_vec_async << endl;
			REQUIRE(p_f_vec_async.maxCoeff() < 0.022);
			REQUIRE(p_f_vec_async.minCoeff() > 0.008);
			cout << "Max index" << p_mesh_ptr->how_many_nodes() -1 << endl;
			start_time = timer.get_milliseconds();
			VectorXd p_f_vec_normal = p_solver.get_f_vec(MAX_INNER);
			int t_normal = timer.get_milliseconds() - start_time;
			cout << "Calculating f_vec normally takes this much time in ms: " << t_normal << endl;
			cout << p_f_vec_normal << endl;
			VectorXd async_diff = p_f_vec_normal - p_f_vec_async;
			REQUIRE(async_diff.maxCoeff() < 0.002);
			REQUIRE(t_async < 0.6*t_normal);
		}*/

		SECTION("Refining and solving should succeed") {
			
			REQUIRE(p_mesh_ptr->get_max_inner_index() == 26);
			REQUIRE(p_mesh_ptr->get_max_outer_index() == 124);
			VectorXd p_sol = p_solver.solve();
			MatrixXd p_grid = p_mesh_dao.get_grid_values(p_mesh_ptr);
			MatrixXd p_sol_values = p_solver_dao.get_solution_values(p_mesh_ptr, p_sol);
			REQUIRE(abs(p_sol.maxCoeff() - 0.05) < 0.01);
			REQUIRE(error_norm(p_sol_values) < 4 * pow(10, - 5));
		}

		SECTION("Testing the speed of Solver") {

			//p_solver.refine();
			cout << timer.get_milliseconds() << endl;
			VectorXd p_sol2 = p_solver.solve();
			cout << timer.get_milliseconds() << endl;
		}

	}

}