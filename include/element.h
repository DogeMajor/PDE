#ifndef ELEMENT_H
#define ELEMENT_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include "node.h"
#include <vector>
#include <utility>
#include <map>
#include <array>
#include <memory>
#include "Function.h"

using namespace std;
using namespace Eigen;


int factorial(int n){
    return (n != 0)? n*factorial(n-1) : 1;
}


//We simply want to use the already existing nodes and don't need to worry about garbage collection.
template <int Dim, int N, typename T>
class Element{

public:
    Element();
    Element(Node<Dim,T> *nod[N], vector <SimplexFunction <T> > funcs);
    Element(const Element &el);//copy constructor
    ~Element();
    void increase_shared_elements();
    void set_indices();
	Node<Dim, T>** get_nodes();
    SimplexFunction<T> get_function(int node_no);
    Node<Dim,T> operator[](int i);
    Element<Dim,N,T>& operator=(Element &el);
    bool operator==(const Element &el) const;
    bool operator!=(const Element &el) const;
    Matrix<double, Dim, Dim> get_simplex_matrix(Element &el) const;
    //vector <Element <Dim,N,T> > divide();
	vector <pair <int[2], T> > get_midpoints();
	Node<Dim, T>** get_midpoint_nodes(vector <pair <int[2], T> > &mid_points);
	vector < Element <Dim, N, T> > get_vertex_els(Node<Dim, T> *mid_nodes[N]);
	//vector < Element <Dim, N, T> > get_inner_els(Node<Dim, T> *mid_nodes[N]);
	map< array<int, 2>, Node<Dim, T>* > get_new_nodes();
    double get_volume() const;
    void show() const;

private:
    Node<Dim,T>* nodes[N];
    vector <SimplexFunction <T> > functions;

};

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(){
    //for(int i=0; i<N; i++){
        //nodes[i] = new Node<Dim,T>;
    //}
    //increase_shared_elements();
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(Node<Dim,T> *nod[N], vector <SimplexFunction <T> > funcs){
    for(int i=0; i<N; i++){//If node has no shared_elements it must be a new one!
        if(nod[i]->get_shared_elements() <= 0) {nodes[i] = new Node<Dim,T>(*nod[i]);}
        else {nodes[i] = nod[i];}
    }
    increase_shared_elements();
    functions = funcs;
}


template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(const Element &el){
    for(int i=0; i<N; i++){
        nodes[i] = el.nodes[i];
    }
    increase_shared_elements();
    functions = el.functions;
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::~Element(){
    int shared_elements = 0;
    for(int i=0; i<N; i++){
        shared_elements = nodes[i]->get_shared_elements();
        nodes[i]->set_shared_elements(shared_elements-1);
        if(shared_elements-1 <= 0){
            delete nodes[i];
        }
    }
    functions.clear();
    cout << "Element destroyed!" << endl;
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::increase_shared_elements(){
    int shared_elements = 0;
    for(int i=0; i<N; i++){
        shared_elements = nodes[i]->get_shared_elements();
        nodes[i]->set_shared_elements(shared_elements+1);
    }
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::set_indices(){
    for(int i=0; i<N; i++){
        nodes[i]->set_index(i);
    }
}

template <int Dim, int N, typename T>
Node<Dim, T> ** Element<Dim, N, T>::get_nodes() {
	return nodes;
}

template <int Dim, int N, typename T>
SimplexFunction<T> Element<Dim,N,T>::get_function(int node_no){
     return functions[node_no];
}

template <int Dim, int N, typename T>
Node<Dim,T> Element<Dim,N,T>::operator[](int i){
    return *(nodes[i]);
}

template <int Dim, int N, typename T>
Element<Dim,N,T>& Element<Dim,N,T>::operator=(Element &el){
    if(*this != el){
        for(int i=0; i<N; i++){
            if(nodes[i]->get_shared_elements() <= 0) {delete nodes[i];}
        }
        for(int i=0; i<N; i++){
            nodes[i] = el.nodes[i];
        }
    }
    functions = el.functions;
    return *this;
}

template <int Dim, int N, typename T>
bool Element<Dim,N,T>::operator==(const Element &el) const{
    bool same_nodes = true;
    for(int i=0; i<N; i++){
        same_nodes = same_nodes && bool(*(nodes[i]) == *(el.nodes[i]));
        }
    bool same_function_coeff = (functions == el.functions);
    bool both_funcs_lacking = ((functions.size() == 0) && (el.functions.size() == 0));
    bool same_funcs = same_function_coeff || both_funcs_lacking;
    return same_nodes && same_funcs;
}

template <int Dim, int N, typename T>
bool Element<Dim,N,T>::operator!=(const Element &el) const{
    return !(*this == el);
}

template <int Dim, int N, typename T>
Matrix<double, Dim, Dim> Element<Dim,N,T>::get_simplex_matrix(Element &el) const{
    //For a simplex, N == Dim+1
    MatrixXd simplex_mat = MatrixXd::Zero(Dim,Dim);
    for(int col=0; col<Dim; col++){
        for(int row=0; row<Dim; row++){
            simplex_mat(row, col) = el[row+1].get_location()[col] - el[row].get_location()[col];
        }
    }
    return simplex_mat;
}

template <int Dim, int N, typename T>
double Element<Dim,N,T>::get_volume() const{
    Element<Dim,N,T> temp = *this;
    MatrixXd simplex_mat = (Dim == N-1)? get_simplex_matrix(temp): MatrixXd::Zero(Dim,Dim);
    return simplex_mat.determinant()/factorial(Dim);
}



/*template <int Dim, int N, typename T>
vector <Element <Dim, N, T> > Element<Dim, N, T>::divide() {
	vector <Element <Dim, N, T> > els;
	Node<Dim, T> *new_nodes[N];

	//Element <Dim, N, T> el_iter(*this);
	els.push_back(Element)
}*/

//1. item of pair gives the original point indices and the second in the midpoint
template <int Dim, int N, typename T>//OK
vector <pair <int[2], T > > Element<Dim, N, T>::get_midpoints() {
	vector <pair <int[2], T> >mid_points;
	pair <int[2], T> mid_point;
	T loc;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			mid_point.first[0] = i;
			mid_point.first[1] = j;
			mid_point.second = 0.5*(nodes[i]->get_location() + nodes[j]->get_location());
			mid_points.push_back(mid_point);
		}
	}
	return mid_points;
}

template <int Dim, int N, typename T> //OK
Node<Dim, T> ** Element<Dim, N, T>::get_midpoint_nodes(vector <pair <int[2], T> > &mid_points) {
	Node<Dim, T> *midpoint_nodes[(Dim*(Dim+1))/2];
	for (int i = 0; i < mid_points.size(); i++) {
		midpoint_nodes[i] = new Node<Dim, T>(mid_points[i].second);
	}
	return midpoint_nodes;
}

template <int Dim, int N, typename T>
vector < Element <Dim, N, T> > Element<Dim, N, T>::get_vertex_els(Node<Dim, T> *mid_nodes[N]) {
	vector <pair <int[2], T> > mid_points = get_midpoints();
	vector < Element <Dim, N, T> > vertex_els;
	for (int i = 0; i < N; i++) {
		cout << "To be coded...";
	}
	return vertex_els;
}

template <int Dim, int N, typename T>
map< array<int, 2>, Node<Dim, T> * > Element<Dim, N, T>::get_new_nodes() {
	//vector <pair <int[2], T> > m_points = get_midpoints();
	//Node<2, T> ** m_nodes = get_midpoint_nodes(m_points);
	map< array<int, 2>, Node<Dim, T> * > node_map;
	T loc;
	
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			//node_map[{i, j}] = m_nodes[i + j - 1];
			loc = 0.5*(nodes[i]->get_location() + nodes[j]->get_location());
			//node_map[{i, j}] = new Node<Dim, VectorXd>(*m_nodes[i + j - 1]);
			Node<Dim, T>* temp_node = new Node<Dim, T>(loc);
			temp_node->show();
			node_map.insert(pair< array<int, 2>, Node<Dim, T>* >({ i, j }, temp_node));
			//node_map[{i, j}] = new Node<Dim, T>(loc);
		}
	}
	return node_map;

}


template <int Dim, int N, typename T>
void Element<Dim,N,T>::show() const{
    cout <<"#elements: " << N << endl;
    for(int i=0; i<N; i++){
        nodes[i]->show();
    }
    cout <<"#elements: " << functions.size() << endl;
    for(int i=0; i<functions.size(); i++){
        cout << "Function coefficients for node no " << i <<endl;
        for(int j=0; j<functions[i].coeff.rows(); j++){
            cout << functions[i].coeff[j] <<" " <<endl;
        }
        cout << endl;
    }
}


template <int Dim, int N, typename T>
class ElementFactory{

public:
    ElementFactory(){}
    ~ElementFactory(){}
    Element<Dim, N, T> build(Node<Dim,T> *nod[N]);
    MatrixXd get_inv_matrix(Node<Dim,T>* nodes[]);
    SimplexFunction<T> build_function(MatrixXd M, int node_no);
    vector <SimplexFunction <T> > build_functions(Node<Dim,T>* nodes[]);
};

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementFactory<Dim,N,T>::build(Node<Dim,T> *nod[N]){
    vector <SimplexFunction <T> > funcs = build_functions(nod);
    Element<Dim, N, T> el(nod, funcs);
    return el;
}

template <int Dim, int N, typename T>
MatrixXd ElementFactory<Dim,N,T>::get_inv_matrix(Node<Dim,T>* nodes[]){
    MatrixXd M(Dim+1, Dim+1);
    for(int node=0; node<Dim+1; node++){
        for(int col=0; col<Dim; col++){
            M(node,col) = nodes[node]->get_location()[col];
        }
        M(node,Dim) = 1;
    }
    return M.inverse();
}

template <int Dim, int N, typename T>
SimplexFunction<T> ElementFactory<Dim,N,T>::build_function(MatrixXd M, int node_no){
    VectorXd unit_vec = VectorXd::Zero(Dim+1);
    unit_vec(node_no) = 1;
    VectorXd coeffs = M*unit_vec;
    SimplexFunction<T> fn_obj;
    fn_obj.coeff = coeffs;
    return fn_obj;
}


template <int Dim, int N, typename T>
vector <SimplexFunction <T> > ElementFactory<Dim,N,T>::build_functions(Node<Dim,T>* nodes[]){
    MatrixXd M_inv = get_inv_matrix(nodes);
    vector <SimplexFunction <T> > functions;
    for(int i=0; i<Dim+1; i++){
        functions.push_back(build_function(M_inv, i));
    }
    return functions;
}



template <int Dim, int N, typename T>
class VectorElement {

public:
	VectorElement();
	VectorElement(vector < Node <Dim, T>* > nodes_vec, vector <SimplexFunction <T> > funcs);
	VectorElement(vector < Node <Dim, T> > &nodes_vec, vector <SimplexFunction <T> > funcs);
	~VectorElement();
	void increase_shared_elements();
	void show() const;

private:
	vector < Node <Dim, T>* > nodes;
	vector <SimplexFunction <T> > functions;

};

template <int Dim, int N, typename T>
VectorElement<Dim, N, T>::VectorElement() {
	//for(int i=0; i<N; i++){
		//nodes[i] = new Node<Dim,T>;
	//}
	//increase_shared_elements();
}

template <int Dim, int N, typename T>
VectorElement<Dim, N, T>::VectorElement(vector < Node<Dim, T>* > nodes_vec, vector <SimplexFunction <T> > funcs) {
	for (int i = 0; i < nodes_vec.size(); i++) { nodes.push_back(nodes_vec[i]); }
	//nodes = nodes_vec;
	increase_shared_elements();
	functions = funcs;
}

template <int Dim, int N, typename T>
VectorElement<Dim, N, T>::VectorElement(vector < Node<Dim, T> > &nodes_vec, vector <SimplexFunction <T> > funcs) {
	for (int i = 0; i < nodes_vec.size(); i++) { nodes.push_back(&nodes_vec[i]); }
	//nodes = nodes_vec;
	increase_shared_elements();
	functions = funcs;
}



template <int Dim, int N, typename T>
VectorElement<Dim, N, T>::~VectorElement() {
	int shared_elements = 0;
	/*for (int i = 0; i < N; i++) {
		shared_elements = nodes[i]->get_shared_elements();
		nodes[i]->set_shared_elements(shared_elements - 1);
		if (shared_elements - 1 <= 0) {
			delete nodes[i];
		}
	}*/
	functions.clear();
	cout << "Element destroyed!" << endl;
}

template <int Dim, int N, typename T>
void VectorElement<Dim, N, T>::increase_shared_elements() {
	int shared_elements = 0;
	for (int i = 0; i < N; i++) {
		shared_elements = nodes[i]->get_shared_elements();
		nodes[i]->set_shared_elements(shared_elements + 1);
	}
}

template <int Dim, int N, typename T>
void VectorElement<Dim, N, T>::show() const {
	cout << "#elements: " << N << endl;
	for (int i = 0; i < N; i++) {
		nodes[i]->show();
	}
	cout << "#elements: " << functions.size() << endl;
	for (int i = 0; i < functions.size(); i++) {
		cout << "Function coefficients for node no " << i << endl;
		for (int j = 0; j < functions[i].coeff.rows(); j++) {
			cout << functions[i].coeff[j] << " " << endl;
		}
		cout << endl;
	}
}



#endif
