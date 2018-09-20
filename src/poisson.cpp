#include "../include/poisson.h"

using namespace std;


Poisson::Poisson(double m_error, VectorXi dims, MatrixXd dom, Function scalar_fn){
    max_error = m_error;
    dimensions = dims;
    domain = dom;
    h = MatrixXd::Zero(dims.rows(), 1);
    for(int i = 0; i < domain.cols(); i++){
        h(i) = (domain(1,i) - domain(0,i))/(dims(i) + 1);
        }
    scalar_func = scalar_fn;
    }

void Poisson::set_matrix(){
    int all_dims = get_all_dims();
    std::vector<T> tripletList;
    tripletList.reserve((2*dimensions.rows()+1)*all_dims);
    double h_factor = 0;
    for(int i = 0; i < h.rows(); i++){
        h_factor += pow(h(i),-2);
        }

    for(int n = 0; n < all_dims; n++){
        tripletList.push_back(T(n,n,2*h_factor));
        for(int dim = 0; dim < dimensions.rows(); dim++){
            int shift_no = get_shift_number(dim);
            int under_vector_len = dimensions(0)*shift_no;
            if(n-shift_no >= 0 and (int(n-shift_no)/under_vector_len == int(n)/under_vector_len)){
                tripletList.push_back(T(n,n-shift_no,-pow(h(dim),-2)));
                }
            if(n+shift_no <= all_dims-1 and (int(n+shift_no)/under_vector_len == int(n)/under_vector_len)){
                tripletList.push_back(T(n,n+shift_no,-pow(h(dim),-2)));
                }
            }
        }

    SparseMatrix<double> mat(all_dims, all_dims);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    A = mat;
}

int Poisson::get_all_dims(){
    return dimensions(dimensions.rows()-1)*get_shift_number(dimensions.rows()-1);
}

int Poisson::get_shift_number(int max_dim){
    int shift_no = 1;
    for(int i = 0; i < max_dim; i++){
        shift_no *= dimensions(i);
        }
    return shift_no;
}

int Poisson::to_index(VectorXi coords){
    int index = 0;
    int shift = 0;
    for(int i = 0; i < coords.rows(); i++){
        shift = get_shift_number(i);
        index += coords(i)*shift;
        }
    return index;
}

VectorXi Poisson::to_coords(int index){
    VectorXi coords(dimensions.rows());
    int denominator = get_shift_number(dimensions.rows()-1);
    for(int dim = dimensions.rows()-1; dim > 0; dim--){
        coords(dim) = (index/denominator);
        index = index%denominator;
        denominator /= dimensions(dim-1);
        }
    coords(0) = index;
    return coords;
}

double Poisson::eval_func(VectorXi coords){
    Vector x(coords.rows());
    for(int i=0; i<coords.rows(); i++){
        x(i) = domain(0,i) + (1+coords(i))*h(i);
        }
    return scalar_func(x);
}

VectorXd Poisson::vectorize_scalar_func(){
    VectorXi coords(dimensions.rows());
    int all_dims = get_all_dims();
    VectorXd f(all_dims);
    for(int n=0; n<all_dims; n++){
        coords = to_coords(n);
        f(n) = eval_func(coords);
        }
    return f;
}

SparseMatrix<double> Poisson::get_diff_matrix() const{
    return A;
}

VectorXd Poisson::derivative(Vector u) const{
    VectorXd v = A*u;
    return v;
}

VectorXd Poisson::solve(){

    VectorXd f = vectorize_scalar_func();
    int all_dims = get_all_dims();
    VectorXd u(all_dims);
    ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
    cg.compute(A);
    u = cg.solve(f);
    //cout << "#iterations:     " << cg.iterations() << endl;
    //cout << "estimated error: " << cg.error()      << endl;
    f = A*u;
    u = cg.solve(f);
    return u;
}

void Poisson::show() const{
    cout << "A: " << endl;
    cout << A << endl;
    cout << "dimensions: " << dimensions << endl;
    cout << "domain: " << domain << endl;
    cout << "h: " << h << endl;
}
