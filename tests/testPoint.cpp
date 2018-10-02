#include "../include/point.h"
#include <math.h>
#include <complex>

using namespace std;


#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test Point class" ) {

    vector<double> val;
    val.push_back(1.0);
    val.push_back(2.0);
    Point <double> point(val);

    SECTION( "Test default constructor" ){
        Point <double> default_point;
        vector <double> default_vec= default_point.get_value();
        cout << default_vec.size();
    }

    SECTION( "Test copy constructor" ){
        Point <double> copied_point(point);
        REQUIRE( copied_point.get_value() == val );
    }

    SECTION( "Test assignment operator" ){
        Point <double> assigned_point = point;
        REQUIRE( assigned_point.get_value() == val );
    }

    SECTION( "Test equal to unequal to operators" ){
        Point <double> same_point(point);
        Point <double> not_same_point;
        REQUIRE( same_point == point );
        REQUIRE( not_same_point != point );
    }

    SECTION( "Test get_functions" ){
        vector <double> value = point.get_value();
        REQUIRE( value[0] == 1.0 );
        REQUIRE( value[1] == 2.0 );
        REQUIRE( point.get_dimension() == 2 );
    }

    SECTION( "Test [] operator" ){

        REQUIRE( point[0] == 1.0 );
        REQUIRE( point[1] == 2.0 );
    }


    SECTION( "Test show()" ){
        point.show();
    }

    SECTION( "Test constructing with complex vector" ){

        vector <complex <double> > val;
        val.push_back(std::complex<double> (1,0));
        val.push_back(std::complex<double> (1,1));//1+i
        Point <complex <double> > c_point(val);
        c_point.show();
    }

}

