#include "../include/poisson.h"


Vector func(Vector x){
    Vector y = MatrixXd::Zero(x.rows(), 1);
    for(int i = 0; i < x.rows(); i++){
        y(i) = (3*x(i) + pow(x(i),2))*exp(x(i));
        }
    return y;
}

double two_D_func(Vector x){
    return -x(0)*x(1) + x(0)*x(0);
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
    poisson.set_matrix();
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
    dims(0) = 3;
    dims(1) = 2;
    MatrixXd domain = MatrixXd::Zero(2, 2);
    domain(0,0) = 0.0;
    domain(1,0) = 1.0;
    domain(0,1) = 0.0;
    domain(1,1) = 1.0;
    Poisson poisson = Poisson(max_error, dims, domain, func, two_D_func);
    //poisson.show();
    poisson.set_matrix();
    Vector zero_vector = MatrixXd::Zero(dims(0), 1);
    Vector test_vector = MatrixXd::Identity(dims(0), 1);

    SECTION( "Test set_matrix" ){
        SparseMatrix <double> mat(3*2,3*2);
        mat = poisson.get_diff_matrix();
        REQUIRE( mat.coeff(1,0) == -16.0 );
        REQUIRE( mat.coeff(1,1) == 50.0 );
        REQUIRE( mat.coeff(1,2) == -16.0 );
        REQUIRE( mat.coeff(1,3) == 0.0 );
        REQUIRE( limit_decimals(mat.coeff(1,4),1) == -9.0d );
    }


    SECTION( "test get_shift_number" ) {
        uint shift_no = poisson.get_shift_number(1);
        REQUIRE( shift_no == 3 );
        shift_no = poisson.get_shift_number(0);
        REQUIRE( shift_no == 1 );

        }


    SECTION( "test get_under_dim" ) {
        uint under_dim = poisson.get_under_dim(1);
        REQUIRE( under_dim == 2 );
        under_dim = poisson.get_under_dim(2);
        REQUIRE( under_dim == 3 );
        }

    SECTION( "test to_coords" ) {
        VectorXi coords = poisson.to_coords(0);
        REQUIRE( coords(0) == 0 );
        REQUIRE( coords(1) == 0 );

        coords = poisson.to_coords(4);
        REQUIRE( coords(0) == 1 );
        REQUIRE( coords(1) == 1 );
        coords = poisson.to_coords(7);

        REQUIRE( coords(0) == 1 );
        REQUIRE( coords(1) == 2 );
        }

    SECTION( "test to_index" ) {
        VectorXi coords(2);
        coords << 0,0;
        uint index = poisson.to_index(coords);
        REQUIRE( index == 0 );
        coords << 1,2;
        index = poisson.to_index(coords);
        REQUIRE( index == 7 );
        }


    SECTION( "test eval_func" ) {
        VectorXi coords(2);
        coords << 0,0;
        double result = poisson.eval_func(coords);
        VectorXd x(2);
        x << 0.25, 1.0/3.0;
        double result_1_1 = two_D_func(x);
        REQUIRE( result == result_1_1 );
        }

    SECTION( "test vectorize_scalar_func" ) {//not OK
        VectorXd f;
        f = poisson.vectorize_scalar_func();
        REQUIRE( f.rows() == 6 );
        VectorXd x(2);
        x << 0.5, 2.0/3.0;
        double result_2_2 = two_D_func(x);
        REQUIRE( f.coeff(4) == result_2_2);
        x << 0.25, 1.0/3.0;
        double result_1_1 = two_D_func(x);
        REQUIRE( f.coeff(0) == result_1_1 );
        }


/*
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

    }*/
}

