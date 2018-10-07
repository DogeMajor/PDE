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
    Node <2,VectorXd> *nodes[3];
    nodes[0] = &node1;
    nodes[1] = &node2;
    nodes[2] = &node3;

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
        SimplexFunction first_fn = el.get_function(0);
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
        SimplexFunction first_fn = factory.build_function(M,0);
        REQUIRE( first_fn(node1.get_location()) == 1 );
        REQUIRE( first_fn(node2.get_location()) == 0 );
        REQUIRE( first_fn(node3.get_location()) == 0 );
    }

    SECTION( "Test build_functions()" ){
        vector <SimplexFunction> fns = factory.build_functions(nodes);
        REQUIRE( fns.size() == 3 );
        REQUIRE( fns[0](node1.get_location()) == 1 );
        REQUIRE( fns[1](node2.get_location()) == 1 );
        REQUIRE( fns[2](node2.get_location()) == 0 );
    }

}

/*
TEST_CASE( "Test FunctionAnalyzer" ) {

    VectorXd location(2);
    location << 0.0, 0.0;
    Node <2,VectorXd> node1(location);
    location << 1.0, 0.0;
    Node <2,VectorXd> node2(location);
    location << 1.0, 1.0;
    Node <2,VectorXd> node3(location);
    Node <2,VectorXd> *nodes[3];
    nodes[0] = &node1;
    nodes[1] = &node2;
    nodes[2] = &node3;
    ElementFactory <2, 3, VectorXd> factory;
    Element <2, 3, VectorXd> el = factory.build(nodes);

    SimplexFunction fn_a = el.get_function(0);
    SimplexFunction fn_b = el.get_function(1);
    SimplexFunction fn_c = el.get_function(2);

    BilinearFunction bl_f;
    bl_f.mat = MatrixXd::Identity(2,2);
    FunctionAnalyzer<2, 3, VectorXd> analyzer(bl_f, fn);

    SECTION( "Test constructing FunctionAnalyzer" ){
        FunctionAnalyzer<2, 3, VectorXd> new_analyzer;
    }

    SECTION( "Test sobolev_dot_product" ){
        REQUIRE( analyzer.sobolev_dot_product(element, fn_a, fn_a) == 0.5 );
        REQUIRE( analyzer.sobolev_dot_product(element, fn_a, fn_b) == -0.5 );
        REQUIRE( analyzer.sobolev_dot_product(element, fn_a, fn_c) == 0.0 );
        REQUIRE( analyzer.sobolev_dot_product(element, fn_b, fn_b) == 1.0 );
        REQUIRE( analyzer.sobolev_dot_product(element, fn_b, fn_c) == -0.5 );
        REQUIRE( analyzer.sobolev_dot_product(element, fn_c, fn_c) == 0.5 );
    }

    SECTION( "Test sobolev_f" ){
        REQUIRE( limit_decimals(analyzer.sobolev_f(element, fn_a),7) == 0.0925925 );
    }

}*/

