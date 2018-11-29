#include "../include/Point.h"
#include "../include/Function.h"
#include "../include/TestingTools.h"
#include "../include/HelpfulTools.h"

using namespace std;

#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test SimplexFunction<T>" ) {
    vector <SimplexFunction <VectorXd> > funcs(3);
    VectorXd coeffs(3);
    coeffs << -1,0,1;
    funcs[0].coeff = coeffs;
    coeffs << 1,-1,0;
    funcs[1].coeff = coeffs;
    coeffs << 0,1,0;
    funcs[2].coeff = coeffs;


    SECTION( "Test default constructor()" ){
        SimplexFunction <VectorXd> func;
        VectorXd c(3);
        c << 1,2,3;
        func.coeff = c;
        REQUIRE( func.coeff == c );
    }

    SECTION( "Test gradient(coords)" ){
        VectorXd loc(2);
        loc << 1,0;
        REQUIRE( funcs[0].gradient(loc) == funcs[0].coeff.head(2) );
        loc << 0,7;
        REQUIRE( funcs[0].gradient(loc) == funcs[0].coeff.head(2) );
    }

	SECTION("Test gradient(coords)") {
		cout << "Gradient of funcs[0] " << funcs[0].gradient() << endl;
		REQUIRE(funcs[0].gradient() == funcs[0].coeff.head(2));
		REQUIRE(funcs[0].gradient() == funcs[0].coeff.head(2));
	}

    SECTION( "Test operator (T coords)" ){
        VectorXd x(2);
        x << 1.0,2.0;
        REQUIRE( funcs[2](x) == 2 );
        REQUIRE( funcs[1](x) == -1 );
    }

    SECTION( "Test operators == and !=" ){
        SimplexFunction <VectorXd> func0;
        func0.coeff = funcs[0].coeff;
        REQUIRE( funcs[0] == funcs[0] );
        REQUIRE( funcs[0] == func0 );
        REQUIRE( funcs[0] != funcs[1] );
    }

    SECTION( "Test default constructor() with T == Point<double>" ){
        SimplexFunction <Point <2, double> > point_func;
        VectorXd fn_coeff(3);
        fn_coeff << 1.5, 0.5, 1.0;
        point_func.coeff = fn_coeff;
        vector<double> vec = {1,2};
        Point <2, double> point(vec);
        REQUIRE( point_func.coeff == fn_coeff );
        REQUIRE( point_func(point) == 3.5 );
        REQUIRE( point_func.gradient(point) == point_func.coeff.head(2) );
    }
}


TEST_CASE( "Test BilinearFunction" ) {
    BilinearFunction form;
    MatrixXd coeffs(2,2);
    coeffs << -1,0,1,1;
    form.mat = coeffs;


    SECTION( "Test constructor" ){
        BilinearFunction bl_f;
        MatrixXd M(2,2);
        M << -1,0,1,1;
        bl_f.mat = M;
        cout << M << endl;
        REQUIRE( bl_f.mat == M );
    }

    SECTION( "Test operator (x,y)" ){
        VectorXd a(2);
        a << 1.0,2.0;
        VectorXd b(2);
        b << 0.0,1.0;
        REQUIRE( form(a,b) == 2 );
    }

	SECTION("Test operator (x,y) with I as mat") {
		BilinearFunction form2;
		form2.mat = MatrixXd::Identity(2,2);
		VectorXd c(2);
		c << 1.0, -2.0;
		VectorXd d(2);
		d << -2.0, 0.0;
		REQUIRE(form2(c, c) == 5);
		REQUIRE(form2(d, d) == 4);
	}

}

TEST_CASE("Test BoundaryConditions") {

	SECTION("BoundaryConditions can be initialized with VectorXd") {
		BoundaryConditions<VectorXd> boundaries;
		boundaries.cond_fn = bound_cond;
		boundaries.is_inside_fn = bound_is_inside;
		boundaries.val = bound_val;
		boundaries.normal = bound_normal;
		boundaries.accuracy = 0.0001;
		VectorXd z(2);
		z << 1, 0;
		VectorXd normal(2);
		normal << 1, 0;
		REQUIRE(boundaries.cond(z) == 1);
		REQUIRE(boundaries.val(z) == 0);
		REQUIRE(boundaries.normal(z) == normal);
		SECTION("cond function gives right results near the boundary") {
			z << 0.00011, 0.5;
			REQUIRE(boundaries.cond(z) == 0);
			z << 0.99989, 0.5;
			REQUIRE(boundaries.cond(z) == 0);
			z << 1.00011, 0.5;
			REQUIRE(boundaries.cond(z) == 0);
			z << 0.5, 1.00011;
			REQUIRE(boundaries.cond(z) == 0);
			z << 1.00009, 0.5;
			REQUIRE(boundaries.cond(z) == 1);
			z << 0.99999, 0.5;
			REQUIRE(boundaries.cond(z) == 1);
		}
	}
	
	SECTION("BoundaryConditions can be initialized with Point <2, double> point(vec); ") {
		vector<double> r = { 1,0.5 };
		Point <2, double> p(r);
		BoundaryConditions<Point <2, double> > point_boundaries;
		point_boundaries.cond_fn = point_bound_cond;
		point_boundaries.is_inside_fn = point_bound_is_inside;
		point_boundaries.val = point_bound_val;
		point_boundaries.normal = point_bound_normal;
		point_boundaries.accuracy = 0.000001;
		REQUIRE(point_boundaries.cond(r)==1);
		REQUIRE(point_boundaries.val(r)==0);
		Point <2, double> normal_p = point_boundaries.normal(p);
		
		REQUIRE(normal_p[0] == 1);
		REQUIRE(normal_p[1] == 0);
		normal_p = point_boundaries.normal(0.5*p);
		REQUIRE(normal_p[0] == 0);
		REQUIRE(normal_p[1] == 0);
		REQUIRE(point_boundaries.is_inside(p) == 0);
		REQUIRE(point_boundaries.is_inside(0*p) == 0);
		REQUIRE(point_boundaries.is_inside(0.1*p) == 1);

	}
}

