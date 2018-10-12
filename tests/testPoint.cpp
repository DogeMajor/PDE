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
    Point <double> point_0;
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
        Point <double> fresh_point;
        fresh_point = point;
        REQUIRE( fresh_point.get_value() == val );
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

	SECTION("Test + and - operators") {
		vector<double> zero(2);
		vector<double> twice(2);
		twice[0] = 2 * val[0];
		twice[1] = 2 * val[1];
		Point <double> similar_point(point);
		Point <double> zero_point = point-similar_point;
		Point <double> twice_point = point + similar_point;
		//Point <double> not_same_point;
		REQUIRE(zero_point.get_value() == zero);
		REQUIRE(twice_point.get_value() == twice);
	}

	SECTION("Test multiplication (with a constant) operators") {
		vector<double> two_times(2);
		vector<double> origo(2);
		two_times[0] = 2 * val[0];
		two_times[1] = 2 * val[1];
		
		Point <double> origo_p = 0*point;
		Point <double> doubled_point = point*2;
		//Point <double> not_same_point;
		REQUIRE(origo_p.get_value() == origo);
		REQUIRE(doubled_point.get_value() == two_times);
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
        REQUIRE( c_point[0].real() == 1 );
        REQUIRE( c_point[0].imag() == 0 );
        REQUIRE( c_point[1].real() == 1 );
        REQUIRE( c_point[1].imag() == 1 );
    }

}

