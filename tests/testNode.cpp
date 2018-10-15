#include "../include/point.h"
#include "../include/node.h"
#include "../C++ libs/eigen/Eigen/Dense"

using namespace std;
using namespace Eigen;


#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test Node template with 3-D double vector from Eigen lib" ) {
    VectorXd location(3);
    location << 1.0, 2.0, 3.0;
    Node <3,VectorXd> node(location);
    location << 1.0, 2.0, 3.0;
    Node <3,VectorXd> similar_node(location);

	SECTION("Test default constructor") {
		Node <3, VectorXd> empty_node;
		Node <3, VectorXd> empty_nodes[2];
		REQUIRE(empty_node.get_shared_elements() == 0);
		REQUIRE(empty_node.get_index() == -1);
		REQUIRE(empty_node.get_location().size() == 0);
		REQUIRE(empty_nodes[1].get_index() == -1);
	}
	
    SECTION( "Test get_location" ){
        VectorXd value = node.get_location();
        REQUIRE( value == location );
    }

    SECTION( "Test get_index" ){
        REQUIRE( node.get_index() == 0 );
    }

    SECTION( "Test get_shared_elements" ){
        REQUIRE( node.get_shared_elements() == 0 );
    }


    SECTION( "Test show()" ){
        node.show();
    }

    SECTION( "Test how_many()" ){
        REQUIRE( 2 == node.how_many() );
    }

    SECTION( "Test set_index" ){
        REQUIRE( node.get_index() == 0 );
        node.set_index(2);
        REQUIRE( node.get_index() == 2 );
    }

    SECTION( "Test set_shared_elements" ){
        REQUIRE( node.get_shared_elements() == 0 );
        node.set_shared_elements(1);
        REQUIRE( node.get_shared_elements() == 1 );
    }

    SECTION( "Test == and != operators" ){
        REQUIRE( node == node );
        node.set_shared_elements(2);
        REQUIRE( node != similar_node );
        similar_node.set_shared_elements(2);
        node.set_index(0);
        REQUIRE( similar_node == node );
    }


    SECTION( "Test assignment operator" ){
        int node_no = node.how_many();
        Node <3,VectorXd> assigned_node = node;
        REQUIRE( node.how_many() == node_no+1 );
        REQUIRE( node.get_index() == assigned_node.get_index() );
        REQUIRE( node.get_shared_elements() == assigned_node.get_shared_elements() );
        REQUIRE( node.get_location() == assigned_node.get_location() );
        assigned_node.set_index(666);
        assigned_node = assigned_node;
        REQUIRE( assigned_node.get_index() == 666 );
    }

    SECTION( "Test copy constructor" ){
        int node_number = node.how_many();
        Node <3,VectorXd> copyed_node(node);
        REQUIRE( node.how_many() == node_number+1 );
        REQUIRE( node.get_index() == copyed_node.get_index() );
        REQUIRE( node.get_shared_elements() == copyed_node.get_shared_elements() );
        REQUIRE( node.get_location() == copyed_node.get_location() );
    }

    SECTION( "Testing constructing with Point <double>" ){
        vector <double> coords = {1.0, 2.0, 3.0};
        Point <double> point1(coords);
        Node <3, Point<double> > node1(point1);
        REQUIRE( node1.get_location().get_value() == point1.get_value() );
        REQUIRE( node1.get_location() == point1 );
    }

}
