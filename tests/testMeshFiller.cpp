#include "../include/Point.h"
#include "../include/Mesh.h"
#include "../include/MeshFiller.h"
#include "../include/DAO.h"
#include "../include/TestingTools.h"
#include <math.h>

using namespace std;

#define CATCH_CONFIG_MAIN
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
	VectorXd temp_loc(3);

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
	}

	SECTION("Generating box vertices should succeed") {
		VectorXd mid_point(3);
		mid_point << 0.5, 0.5, 0.5;
		VectorXd side_lengths(3);
		side_lengths << 1, 1, 1;
		VectorXd last_vert = side_lengths;
		vector <Vertex<3, VectorXd>* > box_verts = mesh_filler.build_box_vertices(mid_point, side_lengths);
		REQUIRE(box_verts[7]->get_location() == last_vert);
		REQUIRE(box_verts.size() == 8);
		for (int i = 0; i < box_verts.size(); i++) {
			delete box_verts[i];
		}
	}

	SECTION("Picking Simplex vertices from box vertices should succeed ") {
		vector <Vertex<3, VectorXd>* > box_vertices = mesh_filler.build_box_vertices(cube_point, cube_lengths);
		vector <Vertex<3, VectorXd>* > simplexverts0 = mesh_filler.get_simplex_vertices(0, box_vertices);
		Element<3, 4, VectorXd> el0 = el_factory.build(simplexverts0);
		
		temp_loc << 1, 0, 0;
		REQUIRE(el0[0].get_location() == temp_loc);
		temp_loc << 0, 0, 1;
		REQUIRE(el0[3].get_location() == temp_loc);
		temp_loc << .25, .25, .25;
		REQUIRE(el0.get_avg_location() == temp_loc);

		int vertice_amount = box_vertices[0]->how_many();
		Element<3, 4, VectorXd> el5 = mesh_filler.build_simplex(5, box_vertices);
		REQUIRE(vertice_amount == box_vertices[0]->how_many());
		temp_loc << 1, 0, 1;
		REQUIRE(el5[0].get_location() == temp_loc);
		temp_loc << .75, .75, .75;
		REQUIRE(el5.get_avg_location() == temp_loc);
	}

	SECTION("Filling mesh with simplices covering a box should succeed"){
		MeshDAO<3, 4, VectorXd> dao = MeshDAO<3, 4, VectorXd>();
		vector <Vertex<3, VectorXd>* > b_vertices = mesh_filler.build_box_vertices(cube_point, cube_lengths);
		
		SECTION("Generating simplices to cover the cube should succeed") {
			vector<Element<3, 4, VectorXd> > els;
			for (int i = 0; i < 6; i++) {
				els.push_back(mesh_filler.build_simplex(i, b_vertices));
			}
			VectorXd avg_loc = VectorXd::Zero(3);
			for (int i = 0; i < 6; i++) {
				REQUIRE(limit_decimals(els[i].get_volume(), 4) == 0.1666);
				avg_loc += els[i].get_avg_location();
			}
			avg_loc = (1 / 6.0)*avg_loc;
			REQUIRE(avg_loc == cube_point);
		}
		SECTION("Filling the Mesh should succeed") {
			Mesh<3, 4, VectorXd> mesh = Mesh<3, 4, VectorXd>();
			Mesh<3, 4, VectorXd>* mesh_ptr = &mesh;
			mesh_filler.build_mesh(b_vertices, mesh_ptr, boundaries);
			mesh_ptr->refine();
			mesh_ptr->reset_indices(boundaries);
			MatrixXd grid = dao.get_grid_values(mesh_ptr);
			REQUIRE(grid.rows() == pow(3,3));
			REQUIRE(grid.row(0)(0,0) == .5);
			REQUIRE(grid.row(0)(0,1) == .5);
			REQUIRE(grid.row(0)(0,2) == .5);
		}
	}
}
