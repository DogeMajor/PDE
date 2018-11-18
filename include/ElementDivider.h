#ifndef ELEMENTDIVIDER_H
#define ELEMENTDIVIDER_H
#include "Element.h"
#include "ElementFactory.h"
#include "Function.h"

template <int Dim, int N, typename T>
class ElementDivider {

public:
	ElementDivider() {}
	ElementDivider(BoundaryConditions<T> bound);
	~ElementDivider() {}
	map<array<int, 2>, T> adjust_midpoints(Element<Dim, N, T> &el, map<array<int, 2>, T> m_map, map<array<int, 2>, int> &edges, int max_inner);
	T find_surface_point(T old_mid_loc, T avg, double accuracy);
	map<array<int, 2>, Vertex<Dim, T>* > get_mid_vertices_map(Element<Dim, N, T> &el, map<array<int, 2>, Vertex<Dim, T>* > &commons, map<array<int, 2>, int> &edges, int max_inner);
	vector <Element <Dim, N, T>* > divide(Element <Dim, N, T>& el, map<array<int, 2>, Vertex<Dim, T>* > &commons, map<array<int, 2>, int> &edges, int max_inner);
	vector <Element <Dim, N, T>* > divide_without_adjusting(Element <Dim, N, T>& el, map<array<int, 2>, Vertex<Dim, T>* > &commons, int max_inner);
	Element<Dim, N, T> get_corner_element(int I, vector <Vertex <Dim, T>* >  midpoint_vertices, Element <Dim, N, T>& el);
	Element<Dim, N, T> get_inner_element(int I, vector <Vertex <Dim, T>* >  midpoint_vertices);

private:
	ElementFactory <Dim, N, T> factory;
	map<array<int, 2>, int> midpoint_indices_map;
	BoundaryConditions<T> boundaries;

};

template <int Dim, int N, typename T>
ElementDivider<Dim, N, T>::ElementDivider(BoundaryConditions<T> bound) {
	factory = ElementFactory <Dim, N, T>();
	boundaries = bound;
	int index = 0;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			midpoint_indices_map.insert(pair< array<int, 2>, int>({ i, j }, index));
			index++;
		}
	}
}

template <int Dim, int N, typename T>
T ElementDivider<Dim, N, T>::find_surface_point(T old_mid_loc, T avg, double error_tolerance) {
	
	double max_error = 1;
	T direction = old_mid_loc - avg;
	T loc = old_mid_loc;
	int n = 1;
	
	while (max_error > error_tolerance) {
		if (boundaries.is_inside(loc)) {
			loc = loc + pow(2, -n)*direction;
		}
		else {
			loc = loc - pow(2, -n)*direction;
		}
		if (int(boundaries.is_inside(loc)) + int(boundaries.is_inside(old_mid_loc)) == 1) {//XOR!!!
			n++;
			//cout << n <<": "<< b.is_inside(old_mid_loc) << b.is_inside(loc) << endl;
			max_error = sqrt(dist_squared<Dim, T>(loc, old_mid_loc));
		}
		old_mid_loc = loc;
	}
	return loc;
}

template <int Dim, int N, typename T>//Changes new midpoints so that they are on the boundary if the element is a boundary element
map<array<int, 2>, T> ElementDivider<Dim, N, T>::adjust_midpoints(Element<Dim, N, T> &el, map<array<int, 2>, T> m_map, map<array<int, 2>, int> &edges, int max_inner) {
	int I, J;
	T avg = el.get_avg_location();
	T surface_location;
	double error_tolerance = boundaries.accuracy / double(20);
	for (map<array<int, 2>, T>::const_iterator iter = m_map.begin(); iter != m_map.end(); iter++) {
		I = iter->first[0];
		J = iter->first[1];
		//if ((I > max_inner) && (J > max_inner) && (!boundaries.cond(iter->second))) {//Meaning that the vertices are at the boundary surface!!
		if (edges[{min(I, J), max(I, J)}] == 1) {
			//surface_location = find_surface_point(iter->second, avg, error_tolerance);
			surface_location = (!boundaries.cond(iter->second))? find_surface_point(iter->second, avg, error_tolerance): iter->second;

			m_map[{I, J}] = surface_location;
		}
	}
	return m_map;
}

template <int Dim, int N, typename T>//Also adds new mid vertices to commons!!
map<array<int, 2>, Vertex<Dim, T>* > ElementDivider<Dim, N, T>::get_mid_vertices_map(Element<Dim, N, T> &el, map< array<int, 2>, Vertex<Dim, T>* > &commons, map<array<int,2>, int> &edges, int max_inner) {
	map<array<int, 2>, Vertex<Dim, T>* > vertices_map;
	map<array<int, 2>, T> m_map = adjust_midpoints(el, el.get_midpoints_map(), edges, max_inner);
	T loc;
	int I, J;
	for (int i = 0; i < N; i++) {
		I = el[i].get_index();
		for (int j = i + 1; j < N; j++) {
			J = el[j].get_index();
			loc = m_map[{I, J}];
			if (commons.count({ I,J }) > 0 ) {
				vertices_map[{I, J}] = commons[{I, J}];
			}
			else if (commons.count({ J,I }) > 0) {
				vertices_map[{I, J}] = commons[{J, I}];
			}
			else {
				vertices_map.insert(pair<array<int, 2>, Vertex<Dim, T>* >({ I,J }, new Vertex<Dim, T>(loc)));
				commons[{I, J}] = vertices_map[{I, J}];
			}
		}
	}
	return vertices_map;
}

template <int Dim, int N, typename T>
vector <Element <Dim, N, T>* > ElementDivider<Dim, N, T>::divide(Element <Dim, N, T>& el, map <array<int, 2>, Vertex<Dim, T>* > &commons, map<array<int, 2>, int> &edges, int max_inner) {
	vector <Element <Dim, N, T>* > els;//( Dim*(Dim + 1)) / 2, nullptr);
	vector <Vertex <Dim, T>* >  midpoint_vertices(N, nullptr);
	map <array<int, 2>, Vertex<Dim, T>* > m_vertices_map = get_mid_vertices_map(el, commons, edges, max_inner);
	int i, j, k;
	for (map< array<int, 2>, Vertex<Dim, T>* >::const_iterator iter = m_vertices_map.begin(); iter != m_vertices_map.end(); iter++) {
		i = el.to_local(iter->first[0]);
		j = el.to_local(iter->first[1]);
		k = midpoint_indices_map[{i, j}];
		midpoint_vertices[k] = iter->second;
	}
	
	for (int i = 0; i < N; i++) {
		els.push_back(new Element <Dim, N, T>(get_corner_element(i, midpoint_vertices, el)));
	}
	
	int inner_els_amount = pow(2, Dim) - N;
	for (int i = 0; i < inner_els_amount; i++) {
		els.push_back(new Element <Dim, N, T>(get_inner_element(i, midpoint_vertices)));
	}
	
	return els;
}

template <int Dim, int N, typename T>
vector <Element <Dim, N, T>* > ElementDivider<Dim, N, T>::divide_without_adjusting(Element <Dim, N, T>& el, map<array<int, 2>, Vertex<Dim, T>* > &commons, int max_inner) {

}


template <int Dim, int N, typename T>
Element<Dim, N, T> ElementDivider<Dim, N, T>::get_corner_element(int k, vector <Vertex <Dim, T>* >  midpoint_vertices, Element <Dim, N, T>& el) {
	vector <Vertex <Dim, T>* > added_vertices;
	int I, J;
	int K = el.to_global(k);
	added_vertices.push_back(&el[k]);
	for (int i = 0; i < N; i++) {
		int I = el.to_global(i);
		for (int j = i + 1; j < N; j++) {
			int I = el.to_global(i);
			if (i == k || j == k) {
				added_vertices.push_back(midpoint_vertices[midpoint_indices_map[{i, j}]]);
			}
		}
	}
	return factory.build(added_vertices);
}

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementDivider<Dim, N, T>::get_inner_element(int I, vector<Vertex <Dim, T>* >  midpoint_vertices) {
	vector <Vertex <Dim, T>* > new_vertices;
	map< array<int, 2>, int> unused_vertices_map;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			
			if (i == I || j == I) { new_vertices.push_back(midpoint_vertices[midpoint_indices_map[{i, j}]]); }
			else { unused_vertices_map[{i, j}] = midpoint_indices_map[{i, j}]; }
		}
	}
	new_vertices.push_back(midpoint_vertices[unused_vertices_map.begin()->second]);//chooses first element
	return factory.build(new_vertices);
}

#endif
