#include "../include/point.h"
#include "../include/node.h"
#include "../include/element.h"
#include "../include/PDE.h"
#include "../include/HelpfulTools.h"
#include <math.h>

using namespace std;
using namespace Eigen;

double f_kern(VectorXd coords){
    return coords.transpose()*coords;
}

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test PDE" ) {

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

    vector <SimplexFunction <VectorXd> > funcs(3);
    VectorXd coeffs(3);
    coeffs << -1,0,1;
    funcs[0].coeff = coeffs;
    coeffs << 1,-1,0;
    funcs[1].coeff = coeffs;
    coeffs << 0,1,0;
    funcs[2].coeff = coeffs;
    Element <2, 3, VectorXd> element(nodes, funcs);
    BilinearFunction bl_fn;
    bl_fn.mat = MatrixXd::Identity(2,2);
    PDE<2, VectorXd>  pde(bl_fn, f_kern);

    SECTION( "Test constructing PDE" ){
        PDE<2, VectorXd>  new_pde(BilinearFunction bl_fn, Function f_kernel);
    }

    SECTION( "Test inner product A(.,.)" ){
        REQUIRE( pde.A(element, funcs[0], funcs[0]) == 0.5 );
        REQUIRE( pde.A(element, funcs[0], funcs[1]) == -0.5 );
        REQUIRE( pde.A(element, funcs[0], funcs[2]) == 0.0 );
        REQUIRE( pde.A(element, funcs[1], funcs[1]) == 1.0 );
        REQUIRE( pde.A(element, funcs[1], funcs[2]) == -0.5 );
        REQUIRE( pde.A(element, funcs[2], funcs[2]) == 0.5 );
    }


    SECTION( "Test inner product with f" ){
        REQUIRE( limit_decimals(pde.f(element, funcs[0]),7) == 0.0925925 );
    }
}