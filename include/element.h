#ifndef ELEMENT_H
#define ELEMENT_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include "node.h"
#include <vector>
#include <map>
#include <array>
#include <math.h>
#include "Function.h"
#include "HelpfulTools.h"
#include "VolumeCalculator.h"

using namespace std;
using namespace Eigen;

//int factorial(int n){
    //return (n != 0)? n*factorial(n-1) : 1;
//}

//We simply want to use the already existing nodes and don't need to worry about garbage collection.
//This is the reason why destructor deletes only for nodes which share 0 elements!
template <int Dim, int N, typename T>
class Element : Counter<Element<Dim, N, T> > {

public:
    Element();
	Element(vector<Node <Dim, T>* > nodes_vec, vector<SimplexFunction <T> > funcs);
    Element(const Element &el);
    ~Element();
    void increase_shared_elements();
	void decrease_shared_elements();
    int set_indices(int index);//Every unique node gets an index bigger than this
	void set_all_indices_to(const int index);
	int set_inner_node_indices(int index, BoundaryConditions<T> conds);
	void set_outer_node_indices_to(const int index, BoundaryConditions<T> conds);
	//int set_outer_node_indices(int index, BoundaryConditions<T> conds);
	int how_many() const;
	vector <Node <Dim, T>* > get_nodes();
    SimplexFunction<T> get_function(int node_no);
    Node<Dim,T>& operator[](int i);
    Element<Dim,N,T>& operator=(const Element &el);
    bool operator==(const Element &el) const;
    bool operator!=(const Element &el) const;
    //Matrix<double, Dim, Dim> get_simplex_matrix(Element &el) const;

	map<array<int, 2>, T> get_midpoints_map();

	vector <pair <int[2], T> > get_midpoints();
	vector <Node <Dim, T>* >  get_midpoint_nodes();
	vector <Node <Dim, T>* >  get_midpoint_nodes(map<array<int, 2>, T> m_points_map);
    double get_volume() const;
	T get_avg_location();
	int nodes_size() const { return nodes.size(); }
    void show() const;

private:
	vector <Node <Dim, T>* > nodes; // with pointers #nodes does not increase when new element is added provided that nodes have already been built
    vector <SimplexFunction <T> > functions;
	VolumeCalculator<Dim, T> volume_calculator;
};


template <int Dim, int N, typename T>
Element<Dim,N,T>::Element()
	: nodes(N, nullptr) {
	volume_calculator = VolumeCalculator<Dim, T>();
}

template <int Dim, int N, typename T>
Element<Dim, N, T>::Element(vector <Node <Dim, T>* > nodes_vec, vector <SimplexFunction <T> > funcs) {
	nodes = nodes_vec;
	increase_shared_elements();
	functions = funcs;
	volume_calculator = VolumeCalculator<Dim, T>();
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(const Element &el){
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
int Element<Dim,N,T>::set_indices(int index){
    for(int i=0; i<N; i++){
		if (nodes[i]->get_index() == -1) {
			nodes[i]->set_index(index+1);
			index++;
		}
    }
	return index;
}

template <int Dim, int N, typename T>
void Element<Dim, N, T>::set_all_indices_to(const int index) {
	for (int i = 0; i < N; i++) { nodes[i]->set_index(index); }
}

template <int Dim, int N, typename T>//Ok!
int Element<Dim, N, T>::set_inner_node_indices(int index, BoundaryConditions<T> conds) {
	for (int i = 0; i < N; i++) {
		if ((nodes[i]->get_index() == -1) && (conds.cond(nodes[i]->get_location())==false)) {
			nodes[i]->set_index(index + 1);
			index++;
		}
	}
	return index;
}

template <int Dim, int N, typename T>//Ok!
void Element<Dim, N, T>::set_outer_node_indices_to(const int index, BoundaryConditions<T> conds) {
	for (int i = 0; i < N; i++) {
		if (conds.cond(nodes[i]->get_location()) == true) {
			nodes[i]->set_index(index);
		}
	}
}

/*template <int Dim, int N, typename T>//Not ok!
int Element<Dim, N, T>::set_outer_node_indices(int index, BoundaryConditions<T> conds) {
	for (int i = 0; i < N; i++) {
		if ((nodes[i]->get_index() == -1) && (conds.cond(nodes[i]->get_location()) == true)) {
			nodes[i]->set_index(index + 1);
			index++;
		}
	}
	return index;
}*/


template <int Dim, int N, typename T>
int Element<Dim, N, T>::how_many() const{
	return objects_alive;
}

template <int Dim, int N, typename T>
vector <Node <Dim, T>* > Element<Dim, N, T>::get_nodes() {
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
double Element<Dim,N,T>::get_volume() const{
	return volume_calculator.get_volume(nodes);
}

template <int Dim, int N, typename T>//OK
map<array<int, 2>, T> Element<Dim, N, T>::get_midpoints_map() {
	map<array<int, 2>, T> midpoints_map;
	int I,J = 0;
	T loc;
	for (int i = 0; i < N; i++) {
		I = nodes[i]->get_index();
		for (int j = i + 1; j < N; j++) {
			J = nodes[j]->get_index();
			loc = 0.5*(nodes[i]->get_location() + nodes[j]->get_location());
			midpoints_map.insert(pair<array<int, 2>, T>({ I, J }, loc));
		}
	}
	return midpoints_map;
}

//1. item of pair gives the original point indices and the second in the midpoint
template <int Dim, int N, typename T>//OK  Local map, indices are not node indices!!!
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
vector <Node <Dim, T>* >  Element<Dim, N, T>::get_midpoint_nodes() {
	vector <pair <int[2], T> >mid_points = get_midpoints();
	vector <Node <Dim, T>* >  midpoint_nodes((Dim*(Dim+1))/2, nullptr);
	for (int i = 0; i < mid_points.size(); i++) {
		midpoint_nodes[i] = new Node<Dim, T>(mid_points[i].second);
	}
	return midpoint_nodes;
}
template <int Dim, int N, typename T> //OK
vector <Node <Dim, T>* >  Element<Dim, N, T>::get_midpoint_nodes(map<array<int, 2>, T> m_points_map) {
	typedef map<array<int, 2>, T>::const_iterator PointsMapIter;
	vector <Node <Dim, T>* >  midpoint_nodes;
	T loc;
	for (PointsMapIter iter = m_points_map.begin();  iter != m_points_map.end(); iter++) {
		loc = iter->second;
		midpoint_nodes.push_back(new Node<Dim, T>(loc));
	}
	return midpoint_nodes;
}

template <int Dim, int N, typename T>
T Element<Dim, N, T>::get_avg_location() {
	T loc = nodes[0]->get_location();
	for (int i = 1; i < N; i++) {
		loc = loc + nodes[i]->get_location();
	}
	return (1 / double(N))*loc;
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::show() const{
    for(int i=0; i<nodes.size(); i++){nodes[i]->show();}
    for(int i=0; i<functions.size(); i++){
        cout << "Function coefficients for node no " << i <<endl;
        for(int j=0; j<functions[i].coeff.rows(); j++){
            cout << functions[i].coeff[j] <<" " <<endl;
        }
        cout << endl;
    }
}

#endif
