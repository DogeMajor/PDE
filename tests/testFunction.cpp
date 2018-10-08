#include "../include/point.h"
#include "../include/Function.h"

using namespace std;
using namespace Eigen;


#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
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
        SimplexFunction <Point <double> > point_func;
        VectorXd fn_coeff(3);
        fn_coeff << 1.5, 0.5, 1.0;
        point_func.coeff = fn_coeff;
        vector<double> vec = {1,2};
        Point <double> point(vec);
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

}



