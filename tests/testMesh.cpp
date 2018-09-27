#include "../include/node.h"
#include "../include/element.h"
#include "../include/mesh.h"

#include <math.h>

using namespace std;



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
    location << 1.0, 1.0;
    Node <2,VectorXd> node3(location);
    Node <2,VectorXd> *nodes[3];
    nodes[0] = &node1;
    nodes[1] = &node2;
    nodes[2] = &node3;
    Element <2,3, VectorXd> el(nodes);
    Element <2,3, VectorXd> *top;
    top = &el;
    Mesh <2,3, VectorXd> mesh(top);
    //cout << mesh.weak_form_element(1,1) << endl;
    //mesh.show();



    SECTION( "Test get_simplex_matrix" ){
        Matrix<double,2,2> mat = mesh.get_simplex_matrix(el);
        cout << mat << endl;
        REQUIRE( mat == MatrixXd::Identity(2,2) );
    }

    SECTION( "Test get_element_volume" ){
        double volume = mesh.get_element_volume(el);
        REQUIRE( volume == 1.0 );
    }


}


