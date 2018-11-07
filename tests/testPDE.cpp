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

//N-dim box's boundary
bool bound_cond(VectorXd coords) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) {return true;}
	}
	return false;
}

double bound_val(VectorXd coords) {
	if (coords[1] == 1.0) { return 1; }
	return 0;
}

VectorXd bound_normal(VectorXd coords) {
	int sz = coords.size();
	VectorXd result = VectorXd::Zero(sz);
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) {
			result(i) = (coords[i] == 0.0) ? -1 : 1;
			return result;
		}
	}
	return result;
}

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test PDE" ) {

	VectorXd location(2);
	location << 0.0, 0.0;
	Node <2, VectorXd> node1(location);
	location << 1.0, 0.0;
	Node <2, VectorXd> node2(location);
	location << 1.0, 1.0;
	Node <2, VectorXd> node3(location);
	vector<Node<2, VectorXd> *> nodes(3, nullptr);
	location << 0.0, 0.0;
	nodes[0] = new Node<2, VectorXd>(location);
	location << 1.0, 0.0;
	nodes[1] = new Node<2, VectorXd>(location);
	location << 1.0, 1.0;
	nodes[2] = new Node<2, VectorXd>(location);
	vector<SimplexFunction <VectorXd> > funcs(3);
	VectorXd coeffs(3);
	coeffs << -1, 0, 1;
	funcs[0].coeff = coeffs;
	coeffs << 1, -1, 0;
	funcs[1].coeff = coeffs;
	coeffs << 0, 1, 0;
	funcs[2].coeff = coeffs;

	Element <2, 3, VectorXd> element(nodes, funcs);

    BilinearFunction bl_fn;
    bl_fn.mat = MatrixXd::Identity(2,2);
    PDE<2, VectorXd>  pde(bl_fn, f_kern);

	BoundaryConditions<VectorXd> boundaries = {bound_cond, bound_val, bound_normal };
	

	SECTION("BoundaryConditions for 2-D box should work") {
		PDE<2, VectorXd>  new_pde(BilinearFunction bl_fn, Function f_kernel);
	}

    SECTION( "Test constructing PDE" ){
		VectorXd temp(2);
		temp << 0.5, 0.5;
		REQUIRE(boundaries.cond(temp)==false);
		REQUIRE(boundaries.val(temp) == 0.0);
		temp << 0.5, 1;
		REQUIRE(boundaries.cond(temp) == true);
		REQUIRE(boundaries.val(temp) == 1.0);
    }

    SECTION( "Test inner product A(.,.)" ){
        REQUIRE( pde.A(element, funcs[0], funcs[0]) == 0.5 );
        REQUIRE( pde.A(element, funcs[0], funcs[1]) == -0.5 );
        REQUIRE( pde.A(element, funcs[0], funcs[2]) == 0.0 );
        REQUIRE( pde.A(element, funcs[1], funcs[1]) == 1.0 );
        REQUIRE( pde.A(element, funcs[1], funcs[2]) == -0.5 );
        REQUIRE( pde.A(element, funcs[2], funcs[2]) == 0.5 );
    }

	SECTION("Test f(.,.)") {
		//cout << pde.f(element, funcs[0]) << endl;
		//cout << pde.f(element, funcs[1]) << endl;
		cout << pde.f(element, funcs[2]) << endl;

	}


    SECTION( "Test inner product with f" ){
        REQUIRE( limit_decimals(pde.f(element, funcs[0]),7) == 0.0925925 );
    }
}
