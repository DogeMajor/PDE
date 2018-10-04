#ifndef PDE_H
#define PDE_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Core"
#include "node.h"
#include "element.h"
#include "FunctionHandler.h"
#include "Mesh.h"
#include <math.h>
#include <functional>

#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"
#include <iostream>
#include <math.h>

using namespace Eigen;

typedef double (* Function)(VectorXd x);
typedef Eigen::Triplet<double> Tri;


class PDE
{
    public:
        PDE(BilinearFunction bl_fn, Function fn);
        void show() const;
    protected:
        BilinearFunction A;
        Function f;
        vector <double> boundary_conds;

};

#endif




#endif
