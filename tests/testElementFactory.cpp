#include "../include/element.h"
#include "../include/point.h"
#include "../include/HelpfulTools.h"

using namespace std;


double fn(VectorXd coords){
    return coords.transpose()*coords;
}

#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test ElementFactory" ) {

    VectorXd location(2);
    location << 0.0, 0.0;
    Node <2,VectorXd> node1(location);
    location << 1.0, 0.0;
    Node <2,VectorXd> node2(location);
    location << 1.0, 1.0;
    Node <2,VectorXd> node3(location);
    vector<Node <2,VectorXd> *> nodes;
    nodes.push_back(&node1);
	nodes.push_back(&node2);
	nodes.push_back(&node3);

    ElementFactory <2, 3, VectorXd> factory;


    SECTION( "Test constructing ElementFactory" ){
        ElementFactory <2, 3, VectorXd> new_factory;
    }

    SECTION( "Test building an element" ){
        Element <2, 3, VectorXd> el = factory.build(nodes);
        location << 0.0, 0.0;
        REQUIRE( el[0].get_location() == location );
        location << 1.0, 1.0;
        REQUIRE( el[2].get_location() == location );
        SimplexFunction<VectorXd> first_fn = el.get_function(0);
        VectorXd first_coeffs(3);
        first_coeffs<< -1.0, 0.0, 1.0;
        REQUIRE( first_coeffs == first_fn.coeff );
    }

    SECTION( "Test get_inv_matrix" ){
        MatrixXd M_inv = factory.get_inv_matrix(nodes);
        MatrixXd expected_M(3,3);
        expected_M << -1,1,0,0,-1,1,1,0,0;
        REQUIRE( expected_M == M_inv );
    }

    SECTION( "Test build_function" ){
        MatrixXd M(3,3);
        M << -1,1,0,0,-1,1,1,0,0;
        SimplexFunction<VectorXd> first_fn = factory.build_function(M,0);
        REQUIRE( first_fn(node1.get_location()) == 1 );
        REQUIRE( first_fn(node2.get_location()) == 0 );
        REQUIRE( first_fn(node3.get_location()) == 0 );
    }

    SECTION( "Test build_functions()" ){
        vector <SimplexFunction<VectorXd> > fns = factory.build_functions(nodes);
        REQUIRE( fns.size() == 3 );
        REQUIRE( fns[0](node1.get_location()) == 1 );
        REQUIRE( fns[1](node2.get_location()) == 1 );
        REQUIRE( fns[2](node2.get_location()) == 0 );
    }

}

