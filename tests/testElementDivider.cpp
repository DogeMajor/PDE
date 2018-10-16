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
	vector<Node <2, VectorXd>* > nodes;
	nodes.push_back(&node1);
	nodes.push_back(&node2);
	nodes.push_back(&node3);
	vector <SimplexFunction <VectorXd> > funcs(3);
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
		//I++;
	}
	show_map(MIDPOINTS_MAP);

	int AMOUNT = 4;
	map< array<int, 2>, int> M_MAP;
	int J = 0;
	for (int i = 0; i < AMOUNT; i++) {
		for (int j = i + 1; j < AMOUNT; j++) {
			M_MAP.insert(pair< array<int, 2>, int>({ i, j }, J));
			J++;
		}
		//J++;
	}
	//show_map(M_MAP);

	SECTION( "Test constructing ElementDivider" ) {
		ElementDivider <2, 3, VectorXd> new_divider;
	}

	SECTION( "One can generate new vertex elements out of old element while refining the mesh" ) {
		vector <Node <2, VectorXd >* > mid_nodes = element.get_midpoint_nodes();
		Element <2, 3, VectorXd > el_AB = divider.get_vertex_element(0, mid_nodes, MIDPOINTS_MAP, element);
		el_AB.show();
	}

}

