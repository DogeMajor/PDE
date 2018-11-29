#include "../include/Point.h"
#include "../include/Vertex.h"
#include "../C++ libs/eigen/Eigen/Dense"

using namespace std;
using namespace Eigen;


#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test vertex template with 3-D double vector from Eigen lib" ) {
    VectorXd location(3);
    location << 1.0, 2.0, 3.0;
    Vertex<3,VectorXd> vertex(location);
	VectorXd sim_location(3);
    sim_location << 1.0, 2.0, 3.0;
    Vertex<3,VectorXd> similar_vertex(sim_location);

	SECTION("Test default constructor") {
		Vertex<3, VectorXd> empty_vertex;
		Vertex<3, VectorXd> empty_vertices[2];
		REQUIRE(empty_vertex.get_shared_elements() == 0);
		REQUIRE(empty_vertex.get_index() == -1);
		REQUIRE(empty_vertex.get_location().size() == 0);
		REQUIRE(empty_vertices[1].get_index() == -1);
	}
	
	SECTION("Test copy constructor") {
		int vertex_number = vertex.how_many();
		Vertex<3, VectorXd> copyed_vertex(vertex);
		REQUIRE(vertex.how_many() == vertex_number + 1);
		REQUIRE(vertex.get_index() == copyed_vertex.get_index());
		REQUIRE(vertex.get_shared_elements() == copyed_vertex.get_shared_elements());
		REQUIRE(vertex.get_location() == copyed_vertex.get_location());
	}

	SECTION("Test set_index") {
		REQUIRE(vertex.get_index() == -1);
		vertex.set_index(2);
		REQUIRE(vertex.get_index() == 2);
	}

	SECTION("Test set_shared_elements") {
		REQUIRE(vertex.get_shared_elements() == 0);
		vertex.set_shared_elements(1);
		REQUIRE(vertex.get_shared_elements() == 1);
	}

    SECTION( "Test get_location" ){
        VectorXd value = vertex.get_location();
		cout << value[0] << endl;
		cout << value[1] << endl;
		cout << value[2] << endl;
        REQUIRE( value == location );
    }

    SECTION( "Test get_index" ){
        REQUIRE( vertex.get_index() == -1 );
    }

    SECTION( "Test get_shared_elements" ){
        REQUIRE( vertex.get_shared_elements() == 0 );
    }

	SECTION("Test how_many()") {
		REQUIRE(2 == vertex.how_many());
	}

    SECTION( "Test == and != operators" ){
        REQUIRE( vertex == vertex );
        vertex.set_shared_elements(2);
        REQUIRE( vertex != similar_vertex );
        similar_vertex.set_shared_elements(2);
        vertex.set_index(-1);
        REQUIRE( similar_vertex == vertex );
    }

    SECTION( "Test assignment operator" ){
        int vertex_no = vertex.how_many();
        Vertex<3,VectorXd> assigned_vertex = vertex;
        REQUIRE( vertex.how_many() == vertex_no+1 );
        REQUIRE( vertex.get_index() == assigned_vertex.get_index() );
        REQUIRE( vertex.get_shared_elements() == assigned_vertex.get_shared_elements() );
        REQUIRE( vertex.get_location() == assigned_vertex.get_location() );
        assigned_vertex.set_index(666);
        assigned_vertex = assigned_vertex;
        REQUIRE( assigned_vertex.get_index() == 666 );
    }

	SECTION("Test show()") {
		cout << "Showing Vertex<3,VectorXd>" << endl;
		vertex.show();
	}

    SECTION( "Testing constructing withPoint <3,double>" ){
        vector <double> coords = {1.0, 2.0, 3.0};
        Point <3,double> point1(coords);

        Vertex<3, Point<3,double> > vertex1(point1);
		cout << "Showing Vertex<3, Point<3,double> >" << endl;
		Point<3, double> point_loc = vertex1.get_location();
		for (int i = 0; i < 3; i++) {
			cout << point_loc[i] << endl;
		}
		
		vertex1.show();
        REQUIRE( vertex1.get_location().get_value() == point1.get_value() );
        REQUIRE( vertex1.get_location() == point1 );
    }

}
