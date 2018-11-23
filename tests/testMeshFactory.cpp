#include "../include/Point.h"
#include "../include/Mesh.h"
#include "../include/MeshFiller.h"
#include "../include/DAO.h"
#include "../include/TestingTools.h"
#include <math.h>

using namespace std;

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


TEST_CASE("MeshFiller should be able to generate meshes") {

	MeshFiller<3, 4, VectorXd> mesh_filler;
	BoxDivisions<3> box_divisions;
	VectorXd cube_point(3);
	cube_point << 0.5, 0.5, 0.5;
	VectorXd cube_lengths(3);
	cube_lengths << 1, 1, 1;
	ElementFactory<3, 4, VectorXd> el_factory;
	BoundaryConditions<VectorXd> boundaries = { bound_cond, bound_is_inside, bound_val, bound_normal, 0.000001 };


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
		vector <Vertex<3, VectorXd>* > box_verts = mesh_filler.build_box_vertices(mid_point, side_lengths);
		REQUIRE(box_verts[7]->get_location() == last_vert);
		REQUIRE(box_verts.size());
		for (int i = 0; i < box_verts.size(); i++) {
			delete box_verts[i];
		}
	}

	SECTION("Picking Simplex vertices from box vertices should succeed ") {
		vector <Vertex<3, VectorXd>* > box_vertices = mesh_filler.build_box_vertices(cube_point, cube_lengths);
		vector <Vertex<3, VectorXd>* > simplexverts0 = mesh_filler.get_simplex_vertices(0, box_vertices);

		Element<3, 4, VectorXd> el0 = el_factory.build(simplexverts0);
		//el0.show();
		
		vector <Vertex<3, VectorXd>* > simplexverts2 = mesh_filler.get_simplex_vertices(2, box_vertices);

		Element<3, 4, VectorXd> el2 = el_factory.build(simplexverts2);
		int vertice_amount = box_vertices[0]->how_many();
		Element<3, 4, VectorXd> el5 = mesh_filler.build_simplex(5, box_vertices);
		REQUIRE(vertice_amount == box_vertices[0]->how_many());

		vector<Element<3, 4, VectorXd> > els;
		for (int i = 0; i < 6; i++) {
			els.push_back(mesh_filler.build_simplex(i, box_vertices));
		}
		//els[5].show();
		Mesh<3, 4, VectorXd> mesh = Mesh<3, 4, VectorXd>();
		Mesh<3, 4, VectorXd>* mesh_ptr = &mesh;
		VectorXd avg_loc = VectorXd::Zero(3);
		for (int i = 0; i < 6; i++) {
			REQUIRE(limit_decimals(els[i].get_volume(),4) == 0.1666);
			avg_loc += els[i].get_avg_location();
			//els[i].show();
			//mesh.push(els[i]);
			//cout << els[i].get_avg_location() << endl;
		}
		avg_loc = (1/6.0)*avg_loc;
		REQUIRE(avg_loc == cube_point);
		
		mesh_filler.build_mesh(box_vertices, mesh_ptr, boundaries);
		mesh_ptr->refine();
		mesh_ptr->reset_indices(boundaries);

		MatrixXd grid = mesh_ptr->get_grid_values();
		MeshDAO<3, 4, VectorXd> mesh_dao;
		Delimiters delimiters;
		map<string, string> mathematica_delimiter = delimiters.get_delimiters(Mathematica);
		cout << mathematica_delimiter["start_row"] << endl;
		MatrixXd test_saving(3, 3);
		test_saving << 1, 2, 3, 4,5,6,7,8,9;
		mesh_dao.save_matrix("testingDAO.txt", test_saving);
		mesh_dao.save_grid("test_grid.txt", mesh_ptr);
		//mesh_ptr->save_matrix("cube_grid.txt", grid);
		//mesh.show();
		//Mesh<3, 4, VectorXd>* same_mesh_ptr = mesh_filler.build_mesh(box_vertices);
		//same_mesh_ptr->show();
		/*Mesh<3, 4, VectorXd> mesh = mesh_filler.build_mesh(box_vertices);
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
		
		vector<Element<3, 4, VectorXd> > simplexes = mesh_filler.cover_box_with_simplices(cube_point, cube_lengths);
		cout << simplexes.size();
		simplexes[0].show();
	}*/

}
