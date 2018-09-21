#include "../include/point.h"
#include <math.h>

using namespace std;


double limit_decimals(double number, int decimals){
    double N = pow(10, decimals);
    return double(int(number * N)) / N;
}


#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test Point class" ) {

    vector<double> val;
    val.push_back(1.0);
    val.push_back(2.0);
    Point <double> point(val);

    SECTION( "Test get_functions" ){
        vector <double> value = point.get_value();
        REQUIRE( value[0] == 1.0 );
        REQUIRE( value[1] == 2.0 );
        int index = point.get_index();
        REQUIRE( index == -1 );
    }

    SECTION( "Test set_functions" ){
        point.set_index(1);
        REQUIRE( 1 == point.get_index() );
    }

}

