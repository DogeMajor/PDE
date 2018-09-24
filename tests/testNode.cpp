#include "../include/point.h"
#include "../include/node.h"
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

    VectorXd location(3);
    location << 1.0, 2.0, 3.0;
    Node <3, VectorXd> node(location);

    SECTION( "Test get_location" ){
        VectorXd value = node.get_location();
        REQUIRE( value(0) == location(0) );
    }

    SECTION( "Test get_index" ){
        REQUIRE( node.get_index() == 0 );
    }

    SECTION( "Test get_neighbour_amount" ){
        REQUIRE( node.get_neighbour_amount() == 0 );
    }


    SECTION( "Test show()" ){
        node.show();
    }

    SECTION( "Test set_index" ){
        REQUIRE( node.get_index() == 0 );
        node.set_index(2);
        REQUIRE( node.get_index() == 2 );
    }

    SECTION( "Test set_neighbours_no" ){
        REQUIRE( node.get_neighbour_amount() == 0 );
        node.set_neighbour_amount(1);
        REQUIRE( node.get_neighbour_amount() == 1 );
    }

}

