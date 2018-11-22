#include "../include/Point.h"
#include "../include/Mesh.h"
#include "../include/MeshFactory.h"
#include "../include/TestingTools.h"
#include <math.h>

using namespace std;


#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


/*TEST_CASE("Test the real Mesh with Elements based on Points") {

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
	//el1.show();
	//el2.show();
	Mesh<2, 3, Point <2, double> > el_mesh(el1);
	el_mesh.push(el2);

	BoundaryConditions<Point<2, double> > boundaries;
	boundaries.cond_fn = point_bound_cond;
	boundaries.is_inside_fn = point_bound_is_inside;
	boundaries.val = point_bound_val;
	boundaries.accuracy = 0.0001;
	el_mesh.set_element_divider(boundaries);
	el_mesh.reset_indices(boundaries);
	//el_mesh.show();



	Mesh<2, 3, Point <2, double> > circle_mesh(el1);
	circle_mesh.push(el2);


	BoundaryConditions<Point<2, double> > circle_boundaries;
	circle_boundaries.cond_fn = p_circle_cond;
	circle_boundaries.is_inside_fn = p_circle_is_inside;
	circle_boundaries.val = p_circle_val;
	circle_boundaries.accuracy = 0.001;
	circle_mesh.set_element_divider(circle_boundaries);
	circle_mesh.reset_indices(circle_boundaries);

	vector<double> corner_vec = {0.5, 1.15};
	Point <2, double> corner(corner_vec);
	REQUIRE(circle_boundaries.is_inside(corner));
	//circle_mesh.show();
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
		//el_mesh.show();
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
			//large_mesh.show();
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
		MatrixXd grid = el_mesh.get_grid_values();
		cout << grid << endl;
		REQUIRE(el_mesh.get_max_inner_index() == 0);
		REQUIRE(el_mesh.get_max_outer_index() == 8);
		el_mesh.refine();
		el_mesh.reset_indices(boundaries);

		REQUIRE(el_mesh.how_many_nodes() == 32);
		cout << el_mesh.get_grid_values() << endl;

		el_mesh.save_matrix("grid.txt", grid);
		auto e_sharings = el_mesh.get_edge_sharings();
		//for (auto iter = e_sharings.begin(); iter != e_sharings.end(); iter++) {
			//cout << iter->first[0] << " " << iter->first[1] << endl;
			//cout << iter->second << endl;
		//}
	}

	SECTION("SEtting edge sharings should work") {
		map<array<int, 2>, int> edge_sharings;
		el_mesh.set_edge_sharings(el1, edge_sharings);
		el_mesh.set_edge_sharings(el2, edge_sharings);
		//for (auto iter = edge_sharings.begin(); iter != edge_sharings.end(); iter++) {
			//cout << iter->first[0] << " " << iter->first[1] << endl;
			//cout << iter->second << endl;
		//}
	}

	SECTION("Setting indices (for Nodes inside of Elements!) in the mesh should succeed") {

		el_mesh.reset_indices(boundaries);
		el_mesh.show();
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
		//circle_mesh.refine();
		cout << "How many nodes totally exist after refinement" << el1.how_many() << endl;

		//circle_mesh.reset_indices(circle_boundaries);
		MatrixXd circle_grid = circle_mesh.get_grid_values();
		cout << circle_grid << endl;
	}


}*/


TEST_CASE("MeshFactory should be able to generate meshes") {

	MeshFactory<3, 4, VectorXd> mesh_factory;
	SimplexDivisions<3> simplex_divisions;
	VectorXd cube_point(3);
	cube_point << 0.5, 0.5, 0.5;
	VectorXd cube_lengths(3);
	cube_lengths << 1, 1, 1;
	ElementFactory<3, 4, VectorXd> el_factory;

	SECTION("Getting box coordinates should succeed") {
		VectorXd mp(3);
		mp << 0.5, 0.5, 0.5;
		VectorXd lengths(3);
		lengths << 1, 1, 1;
		MatrixXd cube_coords = get_box_coordinates<3>(mp, lengths);
		REQUIRE(cube_coords(7, 0) == 1);
		REQUIRE(cube_coords(7, 1) == 1);
		REQUIRE(cube_coords(7, 2) == 1);
		REQUIRE(cube_coords.rows() == 8);
		REQUIRE(cube_coords.cols() == 3);
		VectorXd mp2(2);
		mp2 << 0.5, 1;
		VectorXd lengths2(2);
		lengths2 << 1, 2;
		MatrixXd box_coords = get_box_coordinates<2>(mp2, lengths2);
		REQUIRE(box_coords(3, 0) == 1);
		REQUIRE(box_coords(3, 1) == 2);
		REQUIRE(box_coords.rows() == 4);
		REQUIRE(box_coords.cols() == 2);
		//cout << cube_coords << endl;
	}

	SECTION("Generating box vertices should succeed") {
		VectorXd mid_point(3);
		mid_point << 0.5, 0.5, 0.5;
		VectorXd side_lengths(3);
		side_lengths << 1, 1, 1;
		VectorXd last_vert = side_lengths;
		vector <Vertex<3, VectorXd>* > box_verts = mesh_factory.build_box_vertices(mid_point, side_lengths);
		REQUIRE(box_verts[7]->get_location() == last_vert);
		REQUIRE(box_verts.size());
		for (int i = 0; i < box_verts.size(); i++) {
			delete box_verts[i];
		}
	}

	SECTION("Picking Simplex vertices from box vertices should succeed ") {
		vector <Vertex<3, VectorXd>* > box_vertices = mesh_factory.build_box_vertices(cube_point, cube_lengths);
		vector <Vertex<3, VectorXd>* > simplexverts0 = mesh_factory.get_simplex_vertices(0, box_vertices);

		Element<3, 4, VectorXd> el0 = el_factory.build(simplexverts0);
		//el0.show();
		
		vector <Vertex<3, VectorXd>* > simplexverts2 = mesh_factory.get_simplex_vertices(2, box_vertices);

		Element<3, 4, VectorXd> el2 = el_factory.build(simplexverts2);
		int vertice_amount = box_vertices[0]->how_many();
		Element<3, 4, VectorXd> el5 = mesh_factory.build_simplex(5, box_vertices);
		REQUIRE(vertice_amount == box_vertices[0]->how_many());

		vector<Element<3, 4, VectorXd> > els;
		for (int i = 0; i < 6; i++) {
			els.push_back(mesh_factory.build_simplex(i, box_vertices));
		}
		//els[5].show();
		VectorXd avg_loc = VectorXd::Zero(3);
		for (int i = 0; i < 6; i++) {
			cout <<endl<< i << ": " <<els[i].get_volume() << endl;
			REQUIRE(limit_decimals(els[i].get_volume(),4) == 0.1666);
			avg_loc += els[i].get_avg_location();
			els[i].show();
			//cout << els[i].get_avg_location() << endl;
		}
		avg_loc = (1/6.0)*avg_loc;
		REQUIRE(avg_loc == cube_point);
		/*Mesh<3, 4, VectorXd> mesh = mesh_factory.build_mesh(box_vertices);
		BoundaryConditions<VectorXd> p_boundaries = { bound_cond, bound_is_inside, bound_val, bound_normal, 0.000001 };
		//ElementDivider <3, 4, VectorXd> p_divider(p_boundaries);

		mesh.set_element_divider(p_boundaries);
		mesh.reset_indices(p_boundaries);
		mesh.show();
		//el5.show();
		//el2.show();

		for (int i = 0; i < box_vertices.size(); i++) {
			delete box_vertices[i];
		}*/
	}

	/*SECTION("Generating Simplex Elements covering a (unit)box [1]^3 should succeed") {
		
		vector<Element<3, 4, VectorXd> > simplexes = mesh_factory.cover_box_with_simplices(cube_point, cube_lengths);
		cout << simplexes.size();
		simplexes[0].show();
	}*/

}
