#include "poisson.h"
#include <math.h>

Vector func(Vector x){
    Vector y = MatrixXd::Zero(x.rows(), 1);
    for(int i = 0; i < x.rows(); i++){
        y(i) = (3*x(i) + pow(x(i),2))*exp(x(i));
        }
    return y;
    }


double limit_decimals(double number, uint decimals){
    double N = pow(10, decimals);
    return double(int(number * N)) / N;
}


#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test cutting off decimals off doubles" ) {
        double number = 1.123456;

        REQUIRE( limit_decimals(number, 3) == 1.123d );
        REQUIRE( limit_decimals(number, 2) == 1.12d );
}
/*
TEST_CASE( "Poisson: 1-dimensional case" ) {
    double max_error = pow(-12, 10);
    VectorXi dims = MatrixXi::Zero(1, 1);
    dims(0) = 5;
    MatrixXd domain = MatrixXd::Zero(2, 1);
    domain(0,0) = 0.0;
    domain(1,0) = 1.0;
    Poisson poisson = Poisson(max_error, dims, domain, func);
    Vector zero_vector = MatrixXd::Zero(dims(0), 1);
    Vector test_vector = MatrixXd::Identity(dims(0), 1);

    SECTION( "solve -method works" ) {
        poisson.show();
        Vector solution = poisson.solve();

        REQUIRE( solution.rows() == 5 );
        REQUIRE ( limit_decimals(solution(1), 4) == 0.3053d );

        poisson.show();
        }


    SECTION( "Poisson::derivative" ) {
        Vector derivative = poisson.derivative(zero_vector);

        REQUIRE( derivative.rows() == 5 );
        REQUIRE( limit_decimals(derivative(2), 10) == 0.0d );

        Vector u = MatrixXd::Zero(dims(0), 1);
        double change = 0.1;
        for(uint i = 0; i < dims(0); i++){
            u(i) = domain(0) + pow((1+i)*change, 2);
            }
        derivative = poisson.derivative(u);

        REQUIRE( limit_decimals(derivative(0), 2) == -0.72d );
        REQUIRE( derivative.rows() == 5 );
        REQUIRE( limit_decimals(derivative(4), 2) == 12.24d);

    }
}
*/

TEST_CASE( "Poisson: 2-dimensional case" ) {
    double max_error = pow(-12, 10);
    VectorXi dims = MatrixXi::Zero(2, 1);
    dims(0) = 5;
    dims(1) = 3;
    MatrixXd domain = MatrixXd::Zero(2, 2);
    domain(0,0) = 0.0;
    domain(1,0) = 1.0;
    domain(0,1) = 0.0;
    domain(1,1) = 1.0;
    Poisson poisson = Poisson(max_error, dims, domain, func);
    Vector zero_vector = MatrixXd::Zero(dims(0), 1);
    Vector test_vector = MatrixXd::Identity(dims(0), 1);

    SECTION( "solve -method works" ) {
        poisson.show();
        Vector solution = poisson.solve();

        REQUIRE( solution.rows() == 5 );
        REQUIRE ( limit_decimals(solution(1), 4) == 0.3053d );

        poisson.show();
        }


    SECTION( "Poisson::derivative" ) {
        Vector derivative = poisson.derivative(zero_vector);

        REQUIRE( derivative.rows() == 5 );
        REQUIRE( limit_decimals(derivative(2), 10) == 0.0d );

        Vector u = MatrixXd::Zero(dims(0), 1);
        double change = 0.1;
        for(uint i = 0; i < dims(0); i++){
            u(i) = domain(0) + pow((1+i)*change, 2);
            }
        derivative = poisson.derivative(u);

        REQUIRE( limit_decimals(derivative(0), 2) == -0.72d );
        REQUIRE( derivative.rows() == 5 );
        REQUIRE( limit_decimals(derivative(4), 2) == 12.24d);

    }
}

