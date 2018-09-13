#include "poisson.h"
using namespace std;


Poisson::Poisson(double m_error, VectorXi dims, MatrixXd dom, Function fn){
    max_error = m_error;
    dimensions = dims;
    domain = dom;
    h = MatrixXd::Zero(dims.rows(), 1);
    for(int i = 0; i < domain.cols(); i++){
        h(i) = (domain(1,i) - domain(0,i))/(dims(i) + 1);
        }
    func = fn;
    //Fix the matrix generation including several dimensions!!!
    std::vector<T> tripletList;
    tripletList.reserve(3*dimensions(0));
    tripletList.push_back(T(0,0,2));
    tripletList.push_back(T(0,1,-1));
    for(int i = 1; i < dimensions(0) - 1; i++){
        tripletList.push_back(T(i,i,2));
        tripletList.push_back(T(i,i-1,-1));
        tripletList.push_back(T(i,i+1,-1));
        }
    tripletList.push_back(T(dimensions(0)-1,dimensions(0)-2,-1));
    tripletList.push_back(T(dimensions(0)-1,dimensions(0)-1,2));
    SparseMatrix<double> mat(dimensions(0), dimensions(0));
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    A = mat;
    }

MatrixXd Poisson::get_diff_matrix() const{
    return A;
}

Vector Poisson::derivative(Vector u) const{
    Vector v = A*u;
    return v*pow(h(0), -2);
}

Vector Poisson::solve(){
    Vector x = MatrixXd::Zero(dimensions(0), 1);
    for(uint i = 0; i < dimensions(0); i++){
        x(i) = domain(0) + (1+i)*h(0);
        }
    Vector f = func(x);

    Vector u = MatrixXd::Zero(dimensions(0), 1);
    // fill A and b
    ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
    cg.compute(A);
    u = cg.solve(f);
    cout << "#iterations:     " << cg.iterations() << endl;
    cout << "estimated error: " << cg.error()      << endl;
    f = A*u;
    u = cg.solve(f);
    return pow(h(0),2)*u;
}

void Poisson::show() const{
    cout << "A: " << endl;
    cout << A << endl;
    cout << "dimension: " << dimensions(0) << endl;
    cout << "domain: " << domain << endl;
    cout << "h: " << h << endl;
}
