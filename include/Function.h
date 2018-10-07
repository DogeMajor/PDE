#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include <math.h>
#include <functional>
#include <vector>

using namespace std;
using namespace Eigen;
typedef double (* Function)(VectorXd x);

struct SimplexFunction{
    VectorXd coeff;
    bool operator==(const SimplexFunction &s) const{
        return (coeff == s.coeff);
    }
    bool operator!=(const SimplexFunction &s) const{
        return (coeff != s.coeff);
    }
    double operator()(VectorXd coords){
        double result = dot_product(coeff,coords);
        return result + double(coeff.tail(1)[0]);
        }
    double dot_product(VectorXd a, VectorXd b){
        int min_size = min(a.size(), b.size());
        double result = 0.0;
        for(int i=0; i<min_size; i++){result += a[i]*b[i];}
        return result;
        }
    VectorXd gradient(VectorXd coords){
        int coeff_size = coeff.size();
        return coeff.head(coeff_size-1);
    }
};


struct BilinearFunction{
    MatrixXd mat;
    double operator()(VectorXd x, VectorXd y){
        double result = (x.transpose())*mat*y;
        return result;
        }
};

