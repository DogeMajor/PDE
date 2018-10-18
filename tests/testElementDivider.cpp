#include "../include/element.h"
#include "../include/point.h"
#include "../include/HelpfulTools.h"

using namespace std;


double fn(VectorXd coords) {
	return coords.transpose()*coords;
}

#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"


TEST_CASE("Test ElementDivider") {

	VectorXd location(2);
	location << 0.0, 0.0;
	Node <2, VectorXd> node1(location);
	location << 1.0, 0.0;
	Node <2, VectorXd> node2(location);
	location << 1.0, 1.0;
	Node <2, VectorXd> node3(location);
	vector<Node<2, VectorXd> *> nodes(3, nullptr);
	location << 0.0, 0.0;
	nodes[0] = new Node<2, VectorXd>(location);
	location << 1.0, 0.0;
	nodes[1] = new Node<2, VectorXd>(location);
	location << 1.0, 1.0;
	nodes[2] = new Node<2, VectorXd>(location);
	//cout << nodes[1]->how_many() << endl;
	vector<SimplexFunction <VectorXd> > funcs(3);
	VectorXd coeffs(3);
	coeffs << -1, 0, 1;
	funcs[0].coeff = coeffs;
	coeffs << 1, -1, 0;
	funcs[1].coeff = coeffs;
	coeffs << 0, 1, 0;
	funcs[2].coeff = coeffs;

	Element <2, 3, VectorXd> element(nodes, funcs);
	ElementDivider <2, 3, VectorXd> divider;

	map< array<int, 2>, int> MIDPOINTS_MAP;
	int I = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			MIDPOINTS_MAP.insert(pair< array<int, 2>, int>({ i, j }, I));
			I++;
		}
	}
	//show_map(MIDPOINTS_MAP);

	SECTION("dist_squared template function should work") {
		VectorXd a(3);
		a << 0, 1, 2;
		VectorXd b(3);
		b << 1, 2, 5;
		REQUIRE(dist_squared<3, VectorXd> (a, b) == 11.0);
		REQUIRE(dist_squared<3, VectorXd> (a, a) == 0.0);
	}


	SECTION("find_smallest (pair) template function should work") {
		VectorXd B(3);
		B << 0.3, 0.2, 0.5;
		pair<int, double> tiniest = find_smallest<VectorXd>(B, B.size());
		cout << tiniest.first << endl;
		REQUIRE(tiniest.first == 1);
		REQUIRE(tiniest.second == 0.2);
		vector <double> dists = { 1.62, 0.82, 0.02 };
		tiniest = find_smallest<vector<double> >(dists, dists.size());
		REQUIRE(tiniest.first == 2);
		REQUIRE(tiniest.second == 0.02);
	}

	SECTION( "Test constructing ElementDivider" ) {
		ElementDivider <2, 3, VectorXd> new_divider;
	}

	SECTION("Generating midpoint_map should succeed") {
		map<array<int, 2>, int> m_map = divider.get_midpoints_map();
		REQUIRE(m_map == MIDPOINTS_MAP);
	}

	//SECTION("Generating location_map should succeed") {
		//map<VectorXd, array<int, 2> > location_map = divider.get_locations_map(element);
		//cout << location_map.size() << endl;
		//REQUIRE(m_map == MIDPOINTS_MAP);
	//}
	

	SECTION("One can generate new vertex elements out of old element when refining the mesh") {
		vector <Node <2, VectorXd >* > mid_nodes = element.get_midpoint_nodes();
		Element <2, 3, VectorXd > el_AB = divider.get_vertex_element(0, mid_nodes, MIDPOINTS_MAP, element);
		REQUIRE(el_AB[0].get_location() == element[0].get_location());
		REQUIRE(el_AB.get_function(0)(el_AB[0].get_location()) == 1);
		REQUIRE(el_AB.get_function(1)(el_AB[0].get_location()) == 0);
		REQUIRE(el_AB[1].get_location() == 0.5*(element[0].get_location()+ element[1].get_location()));
		REQUIRE(el_AB[2].get_location() == 0.5*(element[0].get_location() + element[2].get_location()));
	}

	SECTION("One can generate new inner elements out of old element when refining the mesh") {
		vector <Node <2, VectorXd >* > m_nodes = element.get_midpoint_nodes();
		Element <2, 3, VectorXd > inner_el = divider.get_inner_element(0, m_nodes, MIDPOINTS_MAP);		
		REQUIRE(inner_el.get_function(0)(inner_el[0].get_location()) == 1);
		REQUIRE(inner_el.get_function(1)(inner_el[0].get_location()) == 0);
		REQUIRE(inner_el[0].get_location() == 0.5*(element[0].get_location() + element[1].get_location()));
		REQUIRE(inner_el[1].get_location() == 0.5*(element[0].get_location() + element[2].get_location()));
		REQUIRE(inner_el[2].get_location() == 0.5*(element[1].get_location() + element[2].get_location()));
	}

	SECTION("Generating new Elements should succeed") {
		vector <Element <2, 3, VectorXd >* > els = divider.divide(element);
		cout << els.size() << endl;
		els[0]->show();
		REQUIRE(els.size() == 4);
		REQUIRE(divider.get_inner_element(0, element.get_midpoint_nodes(), MIDPOINTS_MAP)[1].get_location() == (*els[3])[1].get_location());
		REQUIRE(divider.get_vertex_element(1, element.get_midpoint_nodes(), MIDPOINTS_MAP, element)[2].get_location() == (*els[1])[2].get_location());
	}

	SECTION("Calculating average location should succeed") {
		vector <Node <2, VectorXd >* > el_nodes = element.get_nodes();
		VectorXd avg_loc = divider.average_location(el_nodes);
		REQUIRE(limit_decimals(avg_loc[0],4) == 0.6666);
		REQUIRE(limit_decimals(avg_loc[1],4) == 0.3333);
	}

	SECTION("Finding nearest node should succeed") {
		vector <Node <2, VectorXd >* > el_nodes = element.get_nodes();
		VectorXd z(2);
		z << 0.9, 0.9;
		Node <2, VectorXd > nearest_node = divider.nearest_node(z, el_nodes);
		REQUIRE(nearest_node.get_location() == el_nodes[2]->get_location());
	}


}

