#include "../include/element.h"
#include "../include/point.h"
#include "../include/mesh.h"
#include "../include/baseMesh.h"
#include "../include/FunctionHandler.h"

#include <math.h>

using namespace std;

double limit_decimals(double number, int decimals){
    double N = pow(10, decimals);
    return double(int(number * N)) / N;
}

double fn(VectorXd coords){
    return coords.transpose()*coords;
}

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"

/*
TEST_CASE( "Test FunctionGenerator" ) {


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

    Element <2, 3, VectorXd> element(nodes);
    element.increase_shared_elements();
    FunctionGenerator<2, 3, VectorXd> gen;


    SECTION( "Test constructing FunctionGenerator" ){
        FunctionGenerator<2, 3, VectorXd> new_gen;
    }

    SECTION( "Test get_inv_matrix" ){
        MatrixXd M_inv = gen.get_inv_matrix(element);
        cout << M_inv << endl;
        REQUIRE( M_inv(0,0) == -1 );
        REQUIRE( M_inv(1,0) == 0 );
        REQUIRE( M_inv(2,0) == 1 );
    }

    SECTION( "Test build_function" ){
        MatrixXd M = gen.get_inv_matrix(element);
        SimplexFunction first_fn = gen.build_function(M,0);
        REQUIRE( first_fn(node1.get_location()) == 1 );
        REQUIRE( first_fn(node2.get_location()) == 0 );
        REQUIRE( first_fn(node3.get_location()) == 0 );
    }


}
*/

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

    Element <2, 3, VectorXd> element(nodes);
    element.increase_shared_elements();
    MatrixXd simplex_mat = element.get_simplex_matrix(element);
    FunctionGenerator<2, 3, VectorXd> gen;
    MatrixXd M = gen.get_inv_matrix(element);
    SimplexFunction fn_a = gen.build_function(M,0);
    SimplexFunction fn_b = gen.build_function(M,1);
    SimplexFunction fn_c = gen.build_function(M,2);
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

}

