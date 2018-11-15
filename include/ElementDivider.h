#ifndef ELEMENTDIVIDER_H
#define ELEMENTDIVIDER_H
#include "element.h"
#include "ElementFactory.h"

template <int Dim, int N, typename T>
class ElementDivider {

public:
	ElementDivider();
	~ElementDivider() {}
	map<array<int, 2>, T> get_midlocation_map(Element<Dim, N, T> &el);
	map< array<int, 2>, Node<Dim, T>* > get_mid_nodes_map(Element<Dim, N, T> &el, map< array<int, 2>, Node<Dim, T>* > &commons);
	vector <Element <Dim, N, T>* > divide(Element <Dim, N, T>& el, map< array<int, 2>, Node<Dim, T>* > &commons);
	Element<Dim, N, T> get_vertex_element(int I, vector <Node <Dim, T>* >  midpoint_nodes, Element <Dim, N, T>& el);
	Element<Dim, N, T> get_inner_element(int I, vector <Node <Dim, T>* >  midpoint_nodes);

private:
	ElementFactory <Dim, N, T> factory;
	map<array<int, 2>, int> midpoint_indices_map;

};

template <int Dim, int N, typename T>
ElementDivider<Dim, N, T>::ElementDivider() {
	factory = ElementFactory <Dim, N, T>();
	int index = 0;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			midpoint_indices_map.insert(pair< array<int, 2>, int>({ i, j }, index));
			index++;
		}
	}
}

template <int Dim, int N, typename T>
map<array<int, 2>, T>  ElementDivider<Dim, N, T>::get_midlocation_map(Element<Dim, N, T> &el) {
	map<array<int, 2>, T> m_map = el.get_midpoints_map();
	return m_map;
}

template <int Dim, int N, typename T>//Also adds new mid nodes to commons!!
map< array<int, 2>, Node<Dim, T>* > ElementDivider<Dim, N, T>::get_mid_nodes_map(Element<Dim, N, T> &el, map< array<int, 2>, Node<Dim, T>* > &commons) {
	map< array<int, 2>, Node<Dim, T>* > nodes_map;
	map<array<int, 2>, T> m_map = el.get_midpoints_map();
	T loc;
	int I, J;
	for (int i = 0; i < N; i++) {
		I = el[i].get_index();
		for (int j = i + 1; j < N; j++) {
			J = el[j].get_index();
			loc = m_map[{I, J}];
			if (commons.count({ I,J }) > 0 ) {
				nodes_map[{I, J}] = commons[{I, J}];
			}
			else if (commons.count({ J,I }) > 0) {
				nodes_map[{I, J}] = commons[{J, I}];
			}
			else {
				nodes_map.insert(pair<array<int, 2>, Node<Dim, T>* >({ I,J }, new Node<Dim, T>(loc)));
				commons[{I, J}] = nodes_map[{I, J}];
			}
		}
	}
	return nodes_map;
}

template <int Dim, int N, typename T>
vector <Element <Dim, N, T>* > ElementDivider<Dim, N, T>::divide(Element <Dim, N, T>& el, map< array<int, 2>, Node<Dim, T>* >  &commons) {
	vector <Element <Dim, N, T>* > els;//( Dim*(Dim + 1)) / 2, nullptr);
	vector <Node <Dim, T>* >  midpoint_nodes(N, nullptr);
	map< array<int, 2>, Node<Dim, T>* > m_nodes_map = get_mid_nodes_map(el, commons);
	int i, j, k;
	for (map< array<int, 2>, Node<Dim, T>* >::const_iterator iter = m_nodes_map.begin(); iter != m_nodes_map.end(); iter++) {
		i = el.to_local(iter->first[0]);
		j = el.to_local(iter->first[1]);
		k = midpoint_indices_map[{i, j}];
		midpoint_nodes[k] = iter->second;
	}
	
	for (int i = 0; i < N; i++) {
		els.push_back(new Element <Dim, N, T>(get_vertex_element(i, midpoint_nodes, el)));
	}
	
	int inner_els_amount = pow(2, Dim) - N;
	for (int i = 0; i < inner_els_amount; i++) {
		els.push_back(new Element <Dim, N, T>(get_inner_element(i, midpoint_nodes)));
	}
	
	return els;
}

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementDivider<Dim, N, T>::get_vertex_element(int k, vector <Node <Dim, T>* >  midpoint_nodes, Element <Dim, N, T>& el) {
	vector <Node <Dim, T>* > added_nodes;
	//added_nodes.push_back(new Node <Dim, T>(el[I]));//Not clear if new should be used...
	int I, J;
	int K = el.to_global(k);
	added_nodes.push_back(&el[k]);
	for (int i = 0; i < N; i++) {
		int I = el.to_global(i);
		for (int j = i + 1; j < N; j++) {
			int I = el.to_global(i);
			if (i == k || j == k) {
				added_nodes.push_back(midpoint_nodes[midpoint_indices_map[{i, j}]]);
			}
		}
	}
	return factory.build(added_nodes);
}

template <int Dim, int N, typename T>
Element<Dim, N, T> ElementDivider<Dim, N, T>::get_inner_element(int I, vector<Node <Dim, T>* >  midpoint_nodes) {
	vector <Node <Dim, T>* > new_nodes;
	map< array<int, 2>, int> unused_nodes_map;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			
			if (i == I || j == I) { new_nodes.push_back(midpoint_nodes[midpoint_indices_map[{i, j}]]); }
			else { unused_nodes_map[{i, j}] = midpoint_indices_map[{i, j}]; }
		}
	}
	new_nodes.push_back(midpoint_nodes[unused_nodes_map.begin()->second]);//chooses first element
	return factory.build(new_nodes);
}

#endif
