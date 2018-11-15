#ifndef ELEMENT_H
#define ELEMENT_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include "Vertex.h"
#include <vector>
#include <map>
#include <array>
#include <math.h>
#include "Function.h"
#include "HelpfulTools.h"
#include "VolumeCalculator.h"

using namespace std;
using namespace Eigen;


struct IndexMaps {
	typedef map<int, int> IndexMap;
	IndexMap map_to_local;
	IndexMap map_to_global;
	void erase_maps() {
		map_to_local.erase(map_to_local.begin(), map_to_local.end());
		map_to_global.erase(map_to_global.begin(), map_to_global.end());
	}

};

//We simply want to use the already existing vertices and don't need to worry about garbage collection.
//This is the reason why destructor deletes only for vertices which share 0 elements!
template <int Dim, int N, typename T>
class Element : Counter<Element<Dim, N, T> > {

public:
    Element();
	Element(vector<Vertex <Dim, T>* > vertices_vec, vector<SimplexFunction <T> > funcs);
    Element(const Element &el);
    ~Element();
    void increase_shared_elements();
	void decrease_shared_elements();
    int set_indices(int index);//Every unique Vertex gets an index bigger than this
	void set_all_indices_to(const int index);
	int set_inner_vertex_indices(int index, BoundaryConditions<T> conds);
	void set_outer_vertex_indices_to(const int index, BoundaryConditions<T> conds);

	int to_local(int global_index) { return index_maps.map_to_local[global_index]; }
	int to_global(int local_index) { return index_maps.map_to_global[local_index]; }

	void set_index_maps();
	void set_f_variation(int n, double var) { f_variations[n] = var; }
	vector<double> get_f_variations() { return f_variations; }
	double get_avg_f_variation() { return sum<double>(f_variations) / double(N); }
	//int set_outer_vertex_indices(int index, BoundaryConditions<T> conds);
	int how_many() const;
	vector <Vertex <Dim, T>* > get_vertices();
    SimplexFunction<T> get_function(int vertex_no);
    Vertex<Dim,T>& operator[](int i);
    Element<Dim,N,T>& operator=(const Element &el);
    bool operator==(const Element &el) const;
    bool operator!=(const Element &el) const;
    //Matrix<double, Dim, Dim> get_simplex_matrix(Element &el) const;

	map<array<int, 2>, T> get_midpoints_map();

	vector <pair <int[2], T> > get_midpoints();
	vector <Vertex <Dim, T>* >  get_midpoint_vertices();
	vector <Vertex <Dim, T>* >  get_midpoint_vertices(map<array<int, 2>, T> m_points_map);
    double get_volume() const;
	T get_avg_location();
	int vertices_size() const { return vertices.size(); }
    void show() const;

private:
	vector <Vertex <Dim, T>* > vertices; // with pointers #vertices does not increase when new element is added provided that vertices have already been built
    vector <SimplexFunction <T> > functions;
	IndexMaps index_maps;
	VolumeCalculator<Dim, T> volume_calculator;
	vector<double> f_variations;//Needed for refine algo in mesh!
};


template <int Dim, int N, typename T>
Element<Dim,N,T>::Element()
	: vertices(N, nullptr), f_variations(N) {
	volume_calculator = VolumeCalculator<Dim, T>();
}

template <int Dim, int N, typename T>
Element<Dim, N, T>::Element(vector <Vertex <Dim, T>* > vertices_vec, vector <SimplexFunction <T> > funcs) {
	vertices = vertices_vec;
	increase_shared_elements();
	functions = funcs;
	volume_calculator = VolumeCalculator<Dim, T>();
	f_variations = vector<double>(N);
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::Element(const Element &el){
	vertices = el.vertices;
    increase_shared_elements();
    functions = el.functions;
	f_variations = el.f_variations;
}

template <int Dim, int N, typename T>
Element<Dim,N,T>::~Element(){
	if (vertices[0] != nullptr) {//If the first Vertex is null then all of them are
		decrease_shared_elements();
		for (int i = 0; i < vertices.size(); i++) {
			if (vertices[i]->get_shared_elements() <= 0) { delete vertices[i]; cout << "Vertex no " << i << " destroyed" << endl; }
		}
	}
	vertices.clear();
    functions.clear();
	f_variations.clear();
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::increase_shared_elements(){
    int shared_elements = 0;
    for(int i=0; i<N; i++){
        shared_elements = vertices[i]->get_shared_elements();
        vertices[i]->set_shared_elements(shared_elements+1);
    }
}

template <int Dim, int N, typename T>
void Element<Dim, N, T>::decrease_shared_elements() {
	int shared_elements = 0;
	for (int i = 0; i < N; i++) {
		shared_elements = vertices[i]->get_shared_elements();
		vertices[i]->set_shared_elements(shared_elements - 1);
	}
}

template <int Dim, int N, typename T>
int Element<Dim,N,T>::set_indices(int index){
    for(int i=0; i<N; i++){
		if (vertices[i]->get_index() == -1) {
			vertices[i]->set_index(index+1);
			index++;
		}
    }
	return index;
}

template <int Dim, int N, typename T>
void Element<Dim, N, T>::set_all_indices_to(const int index) {
	for (int i = 0; i < N; i++) { vertices[i]->set_index(index); }
}

template <int Dim, int N, typename T>//Ok!
int Element<Dim, N, T>::set_inner_vertex_indices(int index, BoundaryConditions<T> conds) {
	for (int i = 0; i < N; i++) {
		if ((vertices[i]->get_index() == -1) && (conds.cond(vertices[i]->get_location())==false)) {
			vertices[i]->set_index(index + 1);
			index++;
		}
	}
	set_index_maps();
	return index;
}

template <int Dim, int N, typename T>//Ok!
void Element<Dim, N, T>::set_outer_vertex_indices_to(const int index, BoundaryConditions<T> conds) {
	for (int i = 0; i < N; i++) {
		if (conds.cond(vertices[i]->get_location()) == true) {
			vertices[i]->set_index(index);
		}
	}
}

template <int Dim, int N, typename T>
void Element<Dim, N, T>::set_index_maps() {
	index_maps.erase_maps();
	int I;
	for (int i = 0; i < N; i++) {
		I = vertices[i]->get_index();
		index_maps.map_to_global[i] = I;
		index_maps.map_to_local[I] = i;
	}
}
/*template <int Dim, int N, typename T>//Not ok!
int Element<Dim, N, T>::set_outer_vertex_indices(int index, BoundaryConditions<T> conds) {
	for (int i = 0; i < N; i++) {
		if ((vertices[i]->get_index() == -1) && (conds.cond(vertices[i]->get_location()) == true)) {
			vertices[i]->set_index(index + 1);
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
vector <Vertex <Dim, T>* > Element<Dim, N, T>::get_vertices() {
	return vertices;
}

template <int Dim, int N, typename T>
SimplexFunction<T> Element<Dim,N,T>::get_function(int vertex_no){
     return functions[vertex_no];
}

template <int Dim, int N, typename T>
Vertex<Dim,T>& Element<Dim,N,T>::operator[](int i){
    return *(vertices[i]);
}

template <int Dim, int N, typename T>
Element<Dim,N,T>& Element<Dim,N,T>::operator=(const Element &el){
    if(*this != el){
		decrease_shared_elements();
        for(int i=0; i<N; i++){
            if(vertices[i]->get_shared_elements() <= 0) {delete vertices[i];}
        }
        for(int i=0; i<N; i++){
            vertices[i] = el.vertices[i];
        }
    }
    functions = el.functions;
    return *this;
}

template <int Dim, int N, typename T>
bool Element<Dim,N,T>::operator==(const Element &el) const{
    bool same_vertices = true;
    for(int i=0; i<N; i++){
        same_vertices = same_vertices && bool(*(vertices[i]) == *(el.vertices[i]));
        }
    bool same_function_coeff = (functions == el.functions);
    bool both_funcs_lacking = ((functions.size() == 0) && (el.functions.size() == 0));
    bool same_funcs = same_function_coeff || both_funcs_lacking;
    return same_vertices && same_funcs;
}

template <int Dim, int N, typename T>
bool Element<Dim,N,T>::operator!=(const Element &el) const{
    return !(*this == el);
}

template <int Dim, int N, typename T>
double Element<Dim,N,T>::get_volume() const{
	return volume_calculator.get_volume(vertices);
}

template <int Dim, int N, typename T>//OK
map<array<int, 2>, T> Element<Dim, N, T>::get_midpoints_map() {
	map<array<int, 2>, T> midpoints_map;
	int I,J = 0;
	T loc;
	for (int i = 0; i < N; i++) {
		I = vertices[i]->get_index();
		for (int j = i + 1; j < N; j++) {
			J = vertices[j]->get_index();
			loc = 0.5*(vertices[i]->get_location() + vertices[j]->get_location());
			midpoints_map.insert(pair<array<int, 2>, T>({ I, J }, loc));
		}
	}
	return midpoints_map;
}

//1. item of pair gives the original point indices and the second in the midpoint
template <int Dim, int N, typename T>//OK  Local map, indices are not vertex indices!!!
vector <pair <int[2], T > > Element<Dim, N, T>::get_midpoints() {
	vector <pair <int[2], T> >mid_points;
	pair <int[2], T> mid_point;
	T loc;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			mid_point.first[0] = i;
			mid_point.first[1] = j;
			mid_point.second = 0.5*(vertices[i]->get_location() + vertices[j]->get_location());
			mid_points.push_back(mid_point);
		}
	}
	return mid_points;
}

template <int Dim, int N, typename T> //OK
vector <Vertex <Dim, T>* >  Element<Dim, N, T>::get_midpoint_vertices() {
	vector <pair <int[2], T> >mid_points = get_midpoints();
	vector <Vertex <Dim, T>* >  midpoint_vertices((Dim*(Dim+1))/2, nullptr);
	for (int i = 0; i < mid_points.size(); i++) {
		midpoint_vertices[i] = new Vertex<Dim, T>(mid_points[i].second);
	}
	return midpoint_vertices;
}
template <int Dim, int N, typename T> //OK
vector <Vertex <Dim, T>* >  Element<Dim, N, T>::get_midpoint_vertices(map<array<int, 2>, T> m_points_map) {
	typedef map<array<int, 2>, T>::const_iterator PointsMapIter;
	vector <Vertex <Dim, T>* >  midpoint_vertices;
	T loc;
	for (PointsMapIter iter = m_points_map.begin();  iter != m_points_map.end(); iter++) {
		loc = iter->second;
		midpoint_vertices.push_back(new Vertex<Dim, T>(loc));
	}
	return midpoint_vertices;
}

template <int Dim, int N, typename T>
T Element<Dim, N, T>::get_avg_location() {
	T loc = vertices[0]->get_location();
	for (int i = 1; i < N; i++) {
		loc = loc + vertices[i]->get_location();
	}
	return (1 / double(N))*loc;
}

template <int Dim, int N, typename T>
void Element<Dim,N,T>::show() const{
    for(int i=0; i<vertices.size(); i++){vertices[i]->show();}
    for(int i=0; i<functions.size(); i++){
        cout << "Function coefficients for Vertex no " << i <<endl;
        for(int j=0; j<functions[i].coeff.rows(); j++){
            cout << functions[i].coeff[j] <<" " <<endl;
        }
        cout << endl;
    }
}

#endif
