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
#include <math.h>
#include "Function.h"
#include "HelpfulTools.h"

using namespace std;
using namespace Eigen;


int factorial(int n){
    return (n != 0)? n*factorial(n-1) : 1;
}


//We simply want to use the already existing nodes and don't need to worry about garbage collection.
template <int Dim, int N, typename T>
class Element : Counter<Element<Dim, N, T> > {

public:
    Element();
	Element(vector< Node <Dim, T>* > nodes_vec, vector<SimplexFunction <T> > funcs);
    Element(const Element &el);//copy constructor
    ~Element();
    void increase_shared_elements();
	void decrease_shared_elements();
    void set_indices();
	int how_many() const;
	vector <Node <Dim, T>* > get_nodes();
    SimplexFunction<T> get_function(int node_no);
    Node<Dim,T>& operator[](int i);
    Element<Dim,N,T>& operator=(const Element &el);
    bool operator==(const Element &el) const;
    bool operator!=(const Element &el) const;
    Matrix<double, Dim, Dim> get_simplex_matrix(Element &el) const;
	vector <pair <int[2], T> > get_midpoints();
	vector < Node <Dim, T>* >  get_midpoint_nodes();
    double get_volume() const;
    void show() const;

private:
	vector < Node <Dim, T>* > nodes; // with pointers #nodes does not increase when new element is added provided that nodes have already been built
    vector <SimplexFunction <T> > functions;
};

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element()
	: nodes(N, nullptr) {
    //for(int i=0; i<N; i++){
        //nodes[i] = new Node<Dim,T>;
    //}
    //increase_shared_elements();
}

template <int Dim, int N, typename T>
Element<Dim, N, T>::Element(vector < Node <Dim, T>* > nodes_vec, vector <SimplexFunction <T> > funcs) {
	nodes = nodes_vec;
	increase_shared_elements();
	functions = funcs;
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(const Element &el){
    //for(int i=0; i<N; i++){nodes[i] = el.nodes[i];}
	nodes = el.nodes;
    increase_shared_elements();
    functions = el.functions;
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::~Element(){
	if (nodes[0] != nullptr) {//If the first node is null then all of them are
		decrease_shared_elements();
		for (int i = 0; i < nodes.size(); i++) {
			if (nodes[i]->get_shared_elements() <= 0) { delete nodes[i]; cout << "Node no " << i << " destroyed" << endl; }
		}
	}
	nodes.clear();
    functions.clear();
    //cout << "Element destroyed!" << endl;
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
void Element<Dim, N, T>::decrease_shared_elements() {
	int shared_elements = 0;
	for (int i = 0; i < N; i++) {
		shared_elements = nodes[i]->get_shared_elements();
		nodes[i]->set_shared_elements(shared_elements - 1);
	}
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::set_indices(){
    for(int i=0; i<N; i++){
        nodes[i]->set_index(i);
    }
}

template <int Dim, int N, typename T>
int Element<Dim, N, T>::how_many() const{
	return objects_alive;
}

template <int Dim, int N, typename T>
vector < Node <Dim, T>* > Element<Dim, N, T>::get_nodes() {
	return nodes;
}

template <int Dim, int N, typename T>
SimplexFunction<T> Element<Dim,N,T>::get_function(int node_no){
     return functions[node_no];
}

template <int Dim, int N, typename T>
Node<Dim,T>& Element<Dim,N,T>::operator[](int i){
    return *(nodes[i]);
}

template <int Dim, int N, typename T>
Element<Dim,N,T>& Element<Dim,N,T>::operator=(const Element &el){
    if(*this != el){
		decrease_shared_elements();
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
vector < Node <Dim, T>* >  Element<Dim, N, T>::get_midpoint_nodes() {
	vector <pair <int[2], T> >mid_points = get_midpoints();
	vector < Node <Dim, T>* >  midpoint_nodes((Dim*(Dim+1))/2, nullptr);
	for (int i = 0; i < mid_points.size(); i++) {
		midpoint_nodes[i] = new Node<Dim, T>(mid_points[i].second);
	}
	return midpoint_nodes;
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::show() const{
    //cout <<"#elements: " << N << endl;
    for(int i=0; i<nodes.size(); i++){nodes[i]->show();}
    //cout <<"#elements: " << functions.size() << endl;
    /*for(int i=0; i<functions.size(); i++){
        cout << "Function coefficients for node no " << i <<endl;
        for(int j=0; j<functions[i].coeff.rows(); j++){
            cout << functions[i].coeff[j] <<" " <<endl;
        }
        cout << endl;
    }*/
}


template <int Dim, int N, typename T>
class ElementFactory{

public:
    ElementFactory(){}
    ~ElementFactory(){}
    Element<Dim, N, T> build(vector < Node <Dim, T>* > nodes_vec);
    MatrixXd get_inv_matrix(vector < Node <Dim, T>* > nodes_vec);
    SimplexFunction<T> build_function(MatrixXd M, int node_no);
    vector <SimplexFunction <T> > build_functions(vector < Node <Dim, T>* > nodes_vec);
};

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementFactory<Dim,N,T>::build(vector < Node <Dim, T>* > nodes_vec){
    vector <SimplexFunction <T> > funcs = build_functions(nodes_vec);
    Element<Dim, N, T> el(nodes_vec, funcs);
    return el;
}

template <int Dim, int N, typename T>
MatrixXd ElementFactory<Dim,N,T>::get_inv_matrix(vector < Node <Dim, T>* > nodes_vec){
    MatrixXd M(Dim+1, Dim+1);
    for(int node=0; node<Dim+1; node++){
        for(int col=0; col<Dim; col++){
            M(node,col) = nodes_vec[node]->get_location()[col];
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
vector <SimplexFunction <T> > ElementFactory<Dim,N,T>::build_functions(vector < Node <Dim, T>* > nodes_vec){
    MatrixXd M_inv = get_inv_matrix(nodes_vec);
    vector <SimplexFunction <T> > functions;
    for(int i=0; i<Dim+1; i++){
        functions.push_back(build_function(M_inv, i));
    }
    return functions;
}


template <int Dim, int N, typename T>
class ElementDivider {

public:
	ElementDivider() { factory = ElementFactory <Dim, N, T>(); }
	~ElementDivider() {}
	vector <Element <Dim, N, T>* > divide(Element <Dim, N, T>& el);
	Element<Dim, N, T> get_vertex_element(int I, vector <Node <Dim, T>* >  midpoint_nodes, map<array<int, 2>, int> midpoints_map, Element <Dim, N, T>& el);
	Element<Dim, N, T> get_inner_element(int I, vector <Node <Dim, T>* >  midpoint_nodes, map<array<int, 2>, int> midpoints_map);
	map<array<int, 2>, int> get_midpoints_map();
	T average_location(vector <Node <Dim, T>* >  chosen_nodes);
	Node <Dim, T>& nearest_node(T location, vector <Node <Dim, T>* >  nodes);

private:
	ElementFactory <Dim, N, T> factory;

};


template <int Dim, int N, typename T>
map<array<int, 2>, int> ElementDivider<Dim, N, T>::get_midpoints_map() {
	map<array<int, 2>, int> midpoints_map;
	int index = 0;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			midpoints_map.insert(pair< array<int, 2>, int>({ i, j }, index));
			index++;
		}
	}
	return midpoints_map;
}

template <int Dim, int N, typename T>
vector <Element <Dim, N, T>* > ElementDivider<Dim, N, T>::divide(Element <Dim, N, T>& el) {
	vector <Element <Dim, N, T>* > els;//( Dim*(Dim + 1)) / 2, nullptr);
	vector <Node <Dim, T>* >  midpoint_nodes = el.get_midpoint_nodes();
	map< array<int, 2>, int> midpoints_map = get_midpoints_map();
	//Diverse grejer
	for (int i = 0; i < N; i++) {
		els.push_back(new Element <Dim, N, T>(get_vertex_element(i, midpoint_nodes, midpoints_map, el)));
	}
	int inner_els_amount = pow(2, Dim) - N;
	for (int i = 0; i < inner_els_amount; i++) {
		els.push_back(new Element <Dim, N, T>(get_inner_element(i, midpoint_nodes, midpoints_map)));
	}
	return els;
}

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementDivider<Dim, N, T>::get_vertex_element(int I, vector <Node <Dim, T>* >  midpoint_nodes, map< array<int, 2>, int> midpoints_map, Element <Dim, N, T>& el) {
	vector < Node <Dim, T>* > added_nodes;
	added_nodes.push_back(new Node <Dim, T>(el[I]));//Not clear if new should be used...
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if (i == I || j == I) {
				added_nodes.push_back(midpoint_nodes[midpoints_map[{i, j}]]);
			}
		}
	}
	return factory.build(added_nodes);
}

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementDivider<Dim, N, T>::get_inner_element(int I, vector <Node <Dim, T>* >  midpoint_nodes, map< array<int, 2>, int> midpoints_map) {
	vector < Node <Dim, T>* > new_nodes;
	map< array<int, 2>, int> unused_nodes_map;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if (i == I || j == I) {new_nodes.push_back(midpoint_nodes[midpoints_map[{i, j}]]);}
			else {unused_nodes_map[{i, j}] = midpoints_map[{i, j}];}
		}
	}
	new_nodes.push_back(midpoint_nodes[unused_nodes_map.begin()->second]);//chooses first element
	return factory.build(new_nodes);
}

template <int Dim, int N, typename T>
T ElementDivider<Dim, N, T>::average_location(vector <Node <Dim, T>* >  chosen_nodes) {
	T loc = chosen_nodes[0]->get_location();
	for (int i = 1; i < chosen_nodes.size(); i++) {
		loc = loc + chosen_nodes[i]->get_location();
	}
	return loc * (1 / double(chosen_nodes.size()));
}

template <int Dim, int N, typename T>
Node <Dim, T>& ElementDivider<Dim, N, T>::nearest_node(T location, vector <Node <Dim, T>* >  nodes) {
	vector<double> distances;
	for (int i = 0; i < nodes.size(); i++) {
		distances.push_back(dist_squared<Dim, T>(nodes[i]->get_location(), location));
	}
	pair<int, double> smallest = find_smallest<vector<double> >(distances, nodes.size());
	return *nodes[smallest.first];
}

#endif
