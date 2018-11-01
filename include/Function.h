#ifndef FUNCTION_H
#define FUNCTION_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include "point.h"
#include <math.h>
#include <functional>
#include <vector>

using namespace std;
using namespace Eigen;
typedef double (* Function)(VectorXd x);

template <typename T>
struct SimplexFunction{
    VectorXd coeff;
    bool operator==(const SimplexFunction &s) const{
        return (coeff == s.coeff);
    }
    bool operator!=(const SimplexFunction &s) const{
        return (coeff != s.coeff);
    }
    double operator()(T coords){
        double result = 0.0;
        for(int i=0; i<coeff.size()-1; i++){result += coeff[i]*coords[i];}
        return result + double(coeff.tail(1)[0]);
        }
    VectorXd gradient(T coords){//Does not depend on coords for simplex functions
        int coeff_size = coeff.size();
        return coeff.head(coeff_size-1);
    }

	VectorXd gradient() {//Does not depend on coords for simplex functions
		int coeff_size = coeff.size();
		return coeff.head(coeff_size-1);
	}
};


struct BilinearFunction{
    MatrixXd mat;
    double operator()(VectorXd x, VectorXd y) const{
        double result = (x.transpose())*mat*y;
        return result;
        }
};


template <typename T>
struct BoundaryConditions {
	typedef bool(*ConditionFn)(T x);
	typedef double(*ValueFn)(T x);
	ConditionFn cond;
	ValueFn val;
};

#endif
