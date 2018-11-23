#include "../include/Point.h"
#include "../include/Mesh.h"
#include "../include/TestingTools.h"
#include "../include/DAO.h"
#include <math.h>

#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"


TEST_CASE("Test the real Mesh with Elements based on Points") {

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
	vector<Vertex<2, Point <2, double>  > * > node_vec(3, nullptr);
	node_vec[0] = new Vertex<2, Point <2, double> >(point1);
	node_vec[1] = new Vertex<2, Point <2, double> >(point2);
	node_vec[2] = new Vertex<2, Point <2, double> >(point3);
	ElementFactory<2, 3, Point <2, double> > factory;
	Element<2, 3, Point <2, double> > el1 = factory.build(node_vec);
	vector<Vertex<2, Point <2, double> > *> node_vec2;//(3, nullptr);
	node_vec2.push_back(node_vec[0]);
	node_vec2.push_back(node_vec[2]);
	node_vec2.push_back(new Vertex<2, Point <2, double> >(point4));
	Element<2, 3, Point <2, double> > el2 = factory.build(node_vec2);

	Mesh<2, 3, Point <2, double> > el_mesh(el1);
	el_mesh.push(el2);
	
	BoundaryConditions<Point<2, double> > boundaries;
	boundaries.cond_fn = point_bound_cond;
	boundaries.is_inside_fn = point_bound_is_inside;
	boundaries.val = point_bound_val;
	boundaries.accuracy = 0.0001;
	el_mesh.set_element_divider(boundaries);
	el_mesh.reset_indices(boundaries);

	Mesh<2, 3, Point <2, double> > circle_mesh(el1);
	circle_mesh.push(el2);
	
	BoundaryConditions<Point<2, double> > circle_boundaries;
	circle_boundaries.cond_fn = p_circle_cond;
	circle_boundaries.is_inside_fn = p_circle_is_inside;
	circle_boundaries.val = p_circle_val;
	circle_boundaries.accuracy = 0.001;
	circle_mesh.set_element_divider(circle_boundaries);
	circle_mesh.reset_indices(circle_boundaries);

	MeshDAO<2, 3, Point <2, double> > el_mesh_dao;

	vector<double> corner_vec = {0.5, 1.15};
	Point <2, double> corner(corner_vec);
	REQUIRE(circle_boundaries.is_inside(corner));
	vector<double> ORIG = { .5, .5 };
	//cout << sqrt(dist_squared<2, vector<double> >(point2.get_value(), ORIG));

	SECTION("Mesh can be initialized with default constructor") {//OK
		Mesh<2, 3, Point <2, double> > empty_mesh;
		REQUIRE(empty_mesh.how_many() == 3);
	}

	SECTION("One can access top element") {//OK
		REQUIRE(el_mesh.get_top() == el2);
	}

	SECTION("One can access last element") {//OK
		REQUIRE(el_mesh.get_last() == el1);
	}

	SECTION("One can pop and push at the top") {//OK
		REQUIRE(el_mesh.pop());
		REQUIRE(el_mesh.get_top() == el1);
		REQUIRE(el_mesh.push(el2));
		REQUIRE(el_mesh.get_top() == el2);
		REQUIRE(el_mesh.get_last() == el1);
	}

	SECTION("One can pop a MeshNode from anywhere in the Mesh") {//OK!!
		MeshNode<Element<2, 3, Point <2, double> > >* top_mesh_node = el_mesh.get_top_mesh_node();
		REQUIRE(el_mesh.pop(top_mesh_node));
		REQUIRE(el_mesh.get_last() == el2);
		top_mesh_node = el_mesh.get_top_mesh_node();
		REQUIRE(el_mesh.pop(top_mesh_node) == false);
		//Rebuild the mesh!!
		REQUIRE(el_mesh.pop());
		REQUIRE(el_mesh.push(el1));
		REQUIRE(el_mesh.push(el2));
		REQUIRE(el_mesh.how_many_nodes() == 2);
		REQUIRE(el_mesh.get_last() == el1);
	}

	SECTION("One can push a MeshNode from anywhere in the Mesh") {//OK!!
		MeshNode<Element<2, 3, Point <2, double> > >* top_mesh_node = el_mesh.get_top_mesh_node();
		REQUIRE(el_mesh.push(top_mesh_node, el2));//Should be pushed at the middle place
		Element<2, 3, Point <2, double> > new_el(el1);
		REQUIRE(el_mesh.push(nullptr, new_el));
		REQUIRE(el_mesh.how_many_nodes() == 4);
		REQUIRE(el_mesh.get_last() == el1);
		//Rebuild the mesh!!
		REQUIRE(el_mesh.pop() == true);
		REQUIRE(el_mesh.pop(top_mesh_node) == true);
		REQUIRE(el_mesh.how_many_nodes() == 2);
	}

	SECTION("One can push and pop a MeshNode from anywhere in the Mesh") {//OK!!
		vector<Vertex<2, Point <2, double>  > * > node_vec_D;
		node_vec_D.push_back(new Vertex<2, Point <2, double> >(point6));
		node_vec_D.push_back(node_vec[1]);
		node_vec_D.push_back(node_vec[2]);
		Element<2, 3, Point <2, double> > el4 = factory.build(node_vec_D);
		
		vector<Vertex<2, Point <2, double>  > * > node_vec_C;
		node_vec_C.push_back(new Vertex<2, Point <2, double> >(point5));
		node_vec_C.push_back(node_vec[1]);
		node_vec_C.push_back(node_vec_D[0]);
		Element<2, 3, Point <2, double> > el3 = factory.build(node_vec_C);

		Mesh<2, 3, Point <2, double> > large_mesh(el1);

		SECTION("Pushing and popping mesh nodes succeeds everywhere in the mesh") {
			large_mesh.push(el4);
			large_mesh.push(large_mesh.get_top_mesh_node(), el3);
			large_mesh.push(large_mesh.get_top_mesh_node()->next, el2);
			REQUIRE(large_mesh.get_element(0) == el4);
			REQUIRE(large_mesh.get_element(1) == el3);
			REQUIRE(large_mesh.get_element(2) == el2);
			REQUIRE(large_mesh.get_element(3) == el1);

			REQUIRE(large_mesh.pop(large_mesh.get_top_mesh_node()->next));
			REQUIRE(large_mesh.get_element(2) == el1);
			REQUIRE(large_mesh.get_element(1) == el3);
			REQUIRE(large_mesh.get_element(0) == el4);

			REQUIRE(large_mesh.pop(large_mesh.get_top_mesh_node()->next->next) == false);
			REQUIRE(large_mesh.pop(large_mesh.get_top_mesh_node()->next));
			REQUIRE(large_mesh.get_element(1) == el3);
			REQUIRE(large_mesh.get_element(0) == el4);

			REQUIRE(large_mesh.pop());
			REQUIRE(large_mesh.get_element(0) == el3);

			REQUIRE(large_mesh.pop());
			REQUIRE(large_mesh.how_many_nodes() == 0);

			REQUIRE(large_mesh.pop() == false);
		}
	}

	SECTION("Refining the mesh should succeed") {
		REQUIRE(el_mesh.how_many_nodes() == 2);
		cout << "How many nodes totally exist" << el1.how_many() <<endl;
		el_mesh.refine();
		cout << "How many nodes totally exist after refinement" << el1.how_many() << endl;
		
		el_mesh.reset_indices(boundaries);
		REQUIRE(el_mesh.get_max_inner_index() == 0);
		REQUIRE(el_mesh.get_max_outer_index() == 8);
		
		REQUIRE(el_mesh.how_many_nodes() == 8);
		REQUIRE(el_mesh.get_element(2)[0].get_location()[0] == 0);
		REQUIRE(el_mesh.get_element(2)[0].get_location()[1] == 1.0);
		REQUIRE(el_mesh.get_element(2)[1].get_location()[0] == 0);
		REQUIRE(el_mesh.get_element(2)[1].get_location()[1] == 0.5);
	}

	SECTION("Refining the mesh with boundary adjustments should succeed") {
		REQUIRE(el_mesh.how_many_nodes() == 2);
		cout << "How many nodes totally exist" << el1.how_many() << endl;
		el_mesh.refine();
		cout << "How many nodes totally exist after refinement" << el1.how_many() << endl;

		el_mesh.reset_indices(boundaries);
		MatrixXd grid = el_mesh_dao.get_grid_values(&el_mesh);
		cout << grid << endl;
		REQUIRE(el_mesh.get_max_inner_index() == 0);
		REQUIRE(el_mesh.get_max_outer_index() == 8);
		el_mesh.refine();
		el_mesh.reset_indices(boundaries);
		REQUIRE(el_mesh.get_max_inner_index() == 8);
		REQUIRE(el_mesh.get_max_outer_index() == 24);
		
		REQUIRE(el_mesh.how_many_nodes() == 32);
		cout << el_mesh_dao.get_grid_values(&el_mesh) << endl;

	}

	SECTION("SEtting edge sharings should work") {
		map<array<int, 2>, int> edge_sharings;
		el_mesh.set_edge_sharings(el1, edge_sharings);
		el_mesh.set_edge_sharings(el2, edge_sharings);
		for (auto iter = edge_sharings.begin(); iter != edge_sharings.end(); iter++) {
			cout << iter->second << " ";
			if (iter != edge_sharings.begin()) {
				REQUIRE(iter->second == 1);
			}
		}
	}

	SECTION("Setting indices (for Nodes inside of Elements!) in the mesh should succeed") {
		
		el_mesh.reset_indices(boundaries);
		REQUIRE(el_mesh.get_last()[1].get_index() == 3);
		REQUIRE(el_mesh.get_top()[0].get_index() == 0);
	}

	SECTION("Filling a circle with triangles in a mesh should succeed") {
		REQUIRE(circle_mesh.how_many_nodes() == 2);
		cout << "How many nodes totally exist" << el1.how_many() << endl;
		circle_mesh.refine();
		circle_mesh.reset_indices(circle_boundaries);
		circle_mesh.refine();
		circle_mesh.reset_indices(circle_boundaries);
		cout << "How many nodes totally exist after refinement" << el1.how_many() << endl;

		MatrixXd circle_grid = el_mesh_dao.get_grid_values(&circle_mesh);
		double dist = 0;
		double radius = sqrt(0.5);

		for (int i =0; i < circle_mesh.get_max_inner_index() + 1; i++) {
			dist = sqrt(pow(circle_grid(i, 0) - 0.5, 2) + pow(circle_grid(i, 1) - 0.5, 2));
			REQUIRE(dist < radius);
		}

		for (int i = circle_mesh.get_max_inner_index()+1; i < circle_grid.rows(); i++) {
			dist = sqrt(pow(circle_grid(i, 0) - 0.5, 2) + pow(circle_grid(i, 1) - 0.5, 2));
			REQUIRE(abs(dist - radius) < circle_boundaries.accuracy);
		}
	}
}


TEST_CASE("Generating a mesh based on 3-D Element<3,4,Vertex<2, VectorXd> > should succeed") {
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

	Mesh<3, 4, VectorXd> p_mesh(p_element);
	p_mesh.set_element_divider(p_boundaries);
	p_mesh.reset_indices(p_boundaries);
	
	MeshDAO<3, 4, VectorXd> p_mesh_dao;

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


	SECTION("Refining mesh should succeed") {
		cout << p_mesh.how_many_nodes() << endl;
		cout << p_mesh.get_max_outer_index() << endl;
		cout << p_mesh_dao.get_grid_values(&p_mesh) << endl;
		p_mesh.refine();
		p_mesh.reset_indices(p_boundaries);
		cout << p_mesh_dao.get_grid_values(&p_mesh) << endl;
		MatrixXd cube_grid = p_mesh_dao.get_grid_values(&p_mesh);
		p_mesh_dao.save_matrix("cube_grid.txt", cube_grid);
		cout << p_mesh.get_max_inner_index() << p_mesh.get_max_outer_index() << endl;
		REQUIRE(p_mesh.get_max_inner_index() == 0);
		REQUIRE(p_mesh.get_max_outer_index() == 9);
	}

}