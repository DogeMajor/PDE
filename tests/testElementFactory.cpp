#include "../include/ElementFactory.h"
//#include "../include/point.h"
//#include "../include/HelpfulTools.h"

using namespace std;

#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test ElementFactory" ) {

	VectorXd location(2);
	location << 0.0, 0.0;
	Vertex<2, VectorXd> vertex1(location);
	location << 1.0, 0.0;
	Vertex<2, VectorXd> vertex2(location);
	location << 1.0, 1.0;
	Vertex<2, VectorXd> vertex3(location);
	vector<Vertex<2, VectorXd> *> vertices(3, nullptr);
	location << 0.0, 0.0;
	vertices[0] = new Vertex<2, VectorXd>(location);
	location << 1.0, 0.0;
	vertices[1] = new Vertex<2, VectorXd>(location);
	location << 1.0, 1.0;
	vertices[2] = new Vertex<2, VectorXd>(location);

    ElementFactory <2, 3, VectorXd> factory;


    SECTION( "Test constructing ElementFactory" ){
        ElementFactory <2, 3, VectorXd> new_factory;
    }

    SECTION( "Test building an element" ){
        Element <2, 3, VectorXd> el = factory.build(vertices);
        location << 0.0, 0.0;
        REQUIRE( el[0].get_location() == location );
        location << 1.0, 1.0;
        REQUIRE( el[2].get_location() == location );
        SimplexFunction<VectorXd> first_fn = el.get_function(0);
        VectorXd first_coeffs(3);
        first_coeffs<< -1.0, 0.0, 1.0;
        REQUIRE( first_coeffs == first_fn.coeff );
    }

	SECTION("Test building several elements with same vertices") {
		Element <2, 3, VectorXd> el_2 = factory.build(vertices);
		location << 1.0, 1.0;
		REQUIRE(el_2[2].get_location() == location);
	}

    SECTION( "Test get_inv_matrix" ){
        MatrixXd M_inv = factory.get_inv_matrix(vertices);
        MatrixXd expected_M(3,3);
        expected_M << -1,1,0,0,-1,1,1,0,0;
        REQUIRE( expected_M == M_inv );
    }

    SECTION( "Test build_function" ){
        MatrixXd M(3,3);
        M << -1,1,0,0,-1,1,1,0,0;
        SimplexFunction<VectorXd> first_fn = factory.build_function(M,0);
        REQUIRE( first_fn(vertex1.get_location()) == 1 );
        REQUIRE( first_fn(vertex2.get_location()) == 0 );
        REQUIRE( first_fn(vertex3.get_location()) == 0 );
    }

    SECTION( "Test build_functions()" ){
        vector <SimplexFunction<VectorXd> > fns = factory.build_functions(vertices);
        REQUIRE( fns.size() == 3 );
        REQUIRE( fns[0](vertex1.get_location()) == 1 );
        REQUIRE( fns[1](vertex2.get_location()) == 1 );
        REQUIRE( fns[2](vertex2.get_location()) == 0 );
    }

	SECTION("Building vertices out of locations vector should succeed" ) {
		vector<VectorXd> locations;
		locations.push_back(vertex1.get_location());
		locations.push_back(vertex2.get_location());
		vector<Vertex<2, VectorXd> > xtr_vertices = factory.build_vertices(locations);
		REQUIRE(xtr_vertices[0].get_location() == locations[0]);
		REQUIRE(xtr_vertices[1].get_location() == locations[1]);
	}

	SECTION("Building Element from locations vector should succeed") {
		vector<VectorXd> locs;
		locs.push_back(vertex1.get_location());
		locs.push_back(vertex2.get_location());
		locs.push_back(vertex3.get_location());
		Element <2, 3, VectorXd> el3 = factory.build(locs);
		el3[0].show();
		vertex1.show();
		REQUIRE(el3[0].get_location() == locs[0]);
		REQUIRE(el3[1].get_location() == locs[1]);
		REQUIRE(el3[2].get_location() == locs[2]);
	}

}
