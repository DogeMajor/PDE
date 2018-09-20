#include "../include/poisson.h"
using namespace std;


double one_D_func(Vector x){
    return (3*x(0) + pow(x(0),2))*exp(x(0));
}

double two_D_func(Vector x){
    return -x(0)*x(1) + x(0)*x(0);
}


double PI = 3.14159;

double three_D_func(Vector x){
    double result = 1;
    for(int i=0; i<x.rows(); i++){
        result *= sin(x(i)*PI);
    }
    return result;
}


double limit_decimals(double number, int decimals){
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


TEST_CASE( "Poisson: 1-dimensional case" ) {
    double max_error = pow(-12, 10);
    VectorXi dims = MatrixXi::Zero(1, 1);
    dims(0) = 5;
    MatrixXd domain = MatrixXd::Zero(2, 1);
    domain(0,0) = 0.0;
    domain(1,0) = 1.0;
    Poisson poisson = Poisson(max_error, dims, domain, one_D_func);
    poisson.set_matrix();
    Vector zero_vector = MatrixXd::Zero(dims(0), 1);
    Vector test_vector = MatrixXd::Identity(dims(0), 1);

    SECTION( "solve -method works" ) {
        //poisson.show();
        Vector solution = poisson.solve();
        REQUIRE( solution.rows() == 5 );
        REQUIRE ( limit_decimals(solution(1), 4) == 0.3053d );
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



TEST_CASE( "Poisson: 3-dimensional case" ) {
    double max_error = pow(-12, 10);
    VectorXi dims = MatrixXi::Zero(3, 1);
    dims(0) = 4;
    dims(1) = 5;
    dims(2) = 6;
    MatrixXd domain = MatrixXd::Zero(2, 3);
    domain(0,0) = 0.0;
    domain(1,0) = 1.0;
    domain(0,1) = 0.0;
    domain(1,1) = 1.0;
    domain(0,2) = 0.0;
    domain(1,2) = 1.0; //A cube
    Poisson poisson = Poisson(max_error, dims, domain, three_D_func);
    //poisson.show();
    poisson.set_matrix();
    int all_dims = dims(0)*dims(1)*dims(2);
    Vector zero_vector = MatrixXd::Zero(all_dims, 1);
    Vector test_vector = MatrixXd::Identity(all_dims, 1);

    SECTION( "Test set_matrix" ){
        SparseMatrix <double> mat(4*5*6,4*5*6);
        mat = poisson.get_diff_matrix();
        REQUIRE( limit_decimals(mat.coeff(1,0),2) == -24.99 ); //h(0)^-2
        REQUIRE( limit_decimals(mat.coeff(1,1),1) == 220.0 );  //2*SUM_i h(i)^-2
        REQUIRE( limit_decimals(mat.coeff(1,2),2) == -24.99 );
        REQUIRE( limit_decimals(mat.coeff(1,3),1) == 0.0 );
        REQUIRE( limit_decimals(mat.coeff(1,1+1*4),1) == -36.0 );
        REQUIRE( limit_decimals(mat.coeff(1,1+1*4*5),1) == -49.0 );
    }


    SECTION( "test to_coords" ) {
        VectorXi coords = poisson.to_coords(0);
        REQUIRE( coords(0) == 0 );
        REQUIRE( coords(1) == 0 );
        REQUIRE( coords(2) == 0 );

        coords = poisson.to_coords(8);
        REQUIRE( coords(0) == 0 );
        REQUIRE( coords(1) == 2 );
        REQUIRE( coords(2) == 0 );

        coords = poisson.to_coords(54);
        REQUIRE( coords(0) == 2 );
        REQUIRE( coords(1) == 3 );
        REQUIRE( coords(2) == 2 );
        }

    SECTION( "test to_index(gives real index, i.e. not only inner points!" ) {
        VectorXi coords(3);
        coords << 0,0,0;
        int index = poisson.to_index(coords);
        REQUIRE( index == 0 );
        coords << 1,2,0;
        index = poisson.to_index(coords);
        REQUIRE( index == 9 );
        }


    SECTION( "test eval_func" ) {
        VectorXi coords(3);
        coords << 0,0,0;
        double result = poisson.eval_func(coords);
        VectorXd x(3);
        x << 0.20, 1.0/6.0, 1.0/7.0;
        double result_1_1_1 = three_D_func(x);
        REQUIRE( result == result_1_1_1 );
        }

    SECTION( "test vectorize_scalar_func" ) {
        VectorXd f;
        f = poisson.vectorize_scalar_func();
        REQUIRE( f.rows() == 4*5*6 );
        VectorXd x(3);
        x << 0.20, 1.0/6.0, 1.0/7.0;
        double result_1_1_1 = three_D_func(x);
        REQUIRE( f.coeff(0) == result_1_1_1);
        }


    SECTION( "solve -method works" ) {
        //max val of solution should be roughly 1/(PI^6*2*3)
        VectorXd solution = poisson.solve();
        REQUIRE( solution.rows() == 4*5*6 );
        REQUIRE ( limit_decimals(solution(1), 8) == 0.00638159 );
        VectorXd deriv = poisson.derivative(solution);
        VectorXd f = poisson.vectorize_scalar_func();
        REQUIRE( limit_decimals(deriv(0), 5) == limit_decimals(f(0), 5) );
        REQUIRE( limit_decimals(deriv(3),5) == limit_decimals(f(3),5) );
        }


    SECTION( "Poisson::derivative" ) {
        //max deriv. for 3-D sin-function should be roughly 59.218
        VectorXd derivative = poisson.derivative(zero_vector);
        REQUIRE( derivative.rows() == 4*5*6 );
        REQUIRE( limit_decimals(derivative(2), 10) == 0.0d );
        VectorXd x(4*5*6);
        x.setZero();
        VectorXd change(3);
        change << 0.20, 1.0/6.0, 1.0/7.0;
        VectorXi coeff(3);
        double PI = 3.14159;
        for(int i = 0; i < 4*5*6; i++){//f(x,y)=sin(PI*x)sin(PI*y)sin(PI*z), max 1 and min 0
            coeff = poisson.to_coords(i);
            x(i) = 1;
            for(int j=0; j<coeff.rows(); j++){
                x(i) *= sin((coeff(j)+1)*change(j)*PI);
                }
            }
        derivative = poisson.derivative(x);
        REQUIRE( limit_decimals(derivative(0), 2) == 3.68d );
        REQUIRE( derivative.rows() == 4*5*6 );
        REQUIRE( limit_decimals(derivative(4), 2) == 6.38d);
    }
}


TEST_CASE( "Poisson: 2-dimensional case" ) {
    double max_error = pow(-12, 10);
    VectorXi dims = MatrixXi::Zero(2, 1);
    dims(0) = 3;
    dims(1) = 5;
    MatrixXd domain = MatrixXd::Zero(2, 2);
    domain(0,0) = 0.0;
    domain(1,0) = 1.0;
    domain(0,1) = 0.0;
    domain(1,1) = 1.0;
    Poisson poisson = Poisson(max_error, dims, domain, two_D_func);
    poisson.set_matrix();
    int all_dims = dims(0)*dims(1);
    Vector zero_vector = MatrixXd::Zero(all_dims, 1);
    Vector test_vector = MatrixXd::Identity(all_dims, 1);


    SECTION( "Test set_matrix" ){
        SparseMatrix <double> mat(3*5,3*5);
        mat = poisson.get_diff_matrix();
        REQUIRE( mat.coeff(1,0) == -16.0 );
        REQUIRE( limit_decimals(mat.coeff(1,1),1) == 104.0 );
        REQUIRE( limit_decimals(mat.coeff(1,2),1) == -16.0 );
        REQUIRE( limit_decimals(mat.coeff(1,3),1) == 0 );
        REQUIRE( limit_decimals(mat.coeff(1,4),1) == -36.0 );
    }


    SECTION( "test get_shift_number" ) {
        uint shift_no = poisson.get_shift_number(1);
        REQUIRE( shift_no == 3 );
        shift_no = poisson.get_shift_number(0);
        REQUIRE( shift_no == 1 );
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


    SECTION( "test eval_func" ) {
        VectorXi coords(2);
        coords << 0,0;
        double result = poisson.eval_func(coords);
        VectorXd x(2);
        x << 0.25, 1.0/6.0;
        double result_1_1 = two_D_func(x);
        REQUIRE( result == result_1_1 );
        }

    SECTION( "test vectorize_scalar_func" ) {
        VectorXd f;
        f = poisson.vectorize_scalar_func();
        REQUIRE( f.rows() == 15 );
        VectorXd x(2);
        x << 0.5, 2.0/3.0;
        double result_2_4 = two_D_func(x);
        REQUIRE( f.coeff(10) == result_2_4);
        x << 0.25, 1.0/3.0;
        double result_1_2 = two_D_func(x);
        REQUIRE( f.coeff(3) == result_1_2 );
        }


    SECTION( "test to_index" ) {
        VectorXi coords(2);
        coords << 0,0;
        uint index = poisson.to_index(coords);
        REQUIRE( index == 0 );
        coords << 2,4;
        index = poisson.to_index(coords);
        REQUIRE( index ==  14);
        }


    SECTION( "solve -method works" ) {
        VectorXd solution = poisson.solve();
        REQUIRE( solution.rows() == 15 );
        REQUIRE ( limit_decimals(solution(1), 8) == 0.00457382d );
        VectorXd deriv = poisson.derivative(solution);
        VectorXd f = poisson.vectorize_scalar_func();
        REQUIRE( limit_decimals(deriv(0), 5) == limit_decimals(f(0), 5) );
        REQUIRE( limit_decimals(deriv(3),5) == limit_decimals(f(3),5) );
        }


    SECTION( "Poisson::derivative" ) {
        VectorXd derivative = poisson.derivative(zero_vector);
        REQUIRE( derivative.rows() == 15 );
        REQUIRE( limit_decimals(derivative(2), 10) == 0.0d );

        VectorXd x(15);
        x.setZero();
        VectorXd change(2);
        change << 0.25, 0.166667;
        VectorXi coeff(2);
        double PI = 3.14159;
        for(int i = 0; i < 15; i++){//f(x,y) =sin(PI*x)sin(PI*y)
            coeff = poisson.to_coords(i);
            x(i) = 1;
            for(int j=0; j<coeff.rows(); j++){
                x(i) *= sin((coeff(j)+1)*change(j)*PI);
                }
            }
        derivative = poisson.derivative(x);
        REQUIRE( limit_decimals(derivative(0), 2) == 6.72d );
        REQUIRE( derivative.rows() == 15 );
        REQUIRE( limit_decimals(derivative(4), 2) == 16.47d);

    }
}
