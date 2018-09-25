#include "../include/node.h"
#include "../include/element.h"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Dense"
#include <math.h>

using namespace std;
using namespace Eigen;

double limit_decimals(double number, int decimals){
    double N = pow(10, decimals);
    return double(int(number * N)) / N;
}


#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test Node template with 3-D double vector from Eigen lib" ) {

    VectorXd location(2);
    location << 0.0, 0.0;
    Node <VectorXd> node1(location);
    location << 1.0, 0.0;
    Node <VectorXd> node2(location);
    Node <VectorXd> *nodes[2];
    nodes[0] = &node1;
    nodes[1] = &node2;
    Element <2, VectorXd> element(nodes);

    SECTION( "Test get_location" ){
        VectorXd value = nodes[0]->get_location();
        REQUIRE( value(0) == 0.0 );
        REQUIRE( value(1) == 0.0 );
        value = nodes[1]->get_location();
        REQUIRE( value == location );
    }

    SECTION( "Test operator[]" ){
        Node <VectorXd> node = element[0];
        node.show();
        REQUIRE( element[0] == node1 );
        REQUIRE( element[0] != node2 );
        //Node <VectorXd> node_val = *(nodes[0]);
        //REQUIRE( element[0] == node_val );
    }

    SECTION( "Test get_neighbour_amount" ){
        REQUIRE( nodes[0]->get_neighbour_amount() == 0 );
    }


    /*SECTION( "Test show()" ){
        nodes[0]->show();
    }*/


}

