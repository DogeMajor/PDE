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


TEST_CASE( "Test Element template containing Node template initiated with 2-D double vector from Eigen lib" ) {

    VectorXd location(2);
    location << 0.0, 0.0;
    Node <2,VectorXd> node1(location);
    location << 1.0, 0.0;
    Node <2,VectorXd> node2(location);
    Node <2,VectorXd> *nodes[2];
    nodes[0] = &node1;
    nodes[1] = &node2;
    nodes[0]->show();
    nodes[1]->show();
    Element <2, 2, VectorXd> element(nodes);


    SECTION( "Test operator []" ){
        Node <2,VectorXd> node = element[0];
        REQUIRE( element[0].get_location() == nodes[0]->get_location() );
        REQUIRE( element[0].get_location() != nodes[1]->get_location() );
    }

    SECTION( "Test show()" ){
        cout << "showing element[0]" << endl;
        element[0].show();
        cout << "showing node1" << endl;
        nodes[0]->show();
    }

    SECTION( "Test copy constructor" ){
        Element <2, 2, VectorXd> copyed_element(element);
        REQUIRE( copyed_element[0] == element[0] );
        REQUIRE( copyed_element[1] == element[1] );
    }

    SECTION( "Test increase_shared_elements" ){
        int shared_els = element[0].get_shared_elements();
        element.increase_shared_elements();
        REQUIRE( shared_els+1 == element[0].get_shared_elements() );
    }

    SECTION( "Test set_indices" ){
        element.set_indices();
        REQUIRE( 0 == element[0].get_index() );
        REQUIRE( 1 == element[1].get_index() );
    }

    SECTION( "Test assignment operator" ){
        Element <2, 2, VectorXd> assigned_element = element;
        REQUIRE( assigned_element[0] == element[0] );
        REQUIRE( assigned_element[1] == element[1] );
    }


    SECTION( "Test operators == and !=" ){
        Element <2, 2, VectorXd> new_element(nodes);
        REQUIRE( new_element == element );
        Node <2,VectorXd> *reflected_nodes[2];
        reflected_nodes[0] = &node2;
        reflected_nodes[1] = &node1;
        Element <2, 2, VectorXd> reflected_element(reflected_nodes);
        REQUIRE( reflected_element != element );
    }

}

