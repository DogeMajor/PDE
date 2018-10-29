#include "../include/element.h"
#include "../include/point.h"
#include "../include/HelpfulTools.h"

using namespace std;


double fn(VectorXd coords) {
	return coords.transpose()*coords;
}

#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"

TEST_CASE("Test ElementDivider with Point template -based Nodes") {
	vector <double> vec1 = { 0.0, 0.0 };
	vector <double> vec2 = { 1.0, 0.0 };
	vector <double> vec3 = { 1.0, 1.0 };
	vector <double> vec4 = { 0.0, 1.0 };
	vector <double> vec5 = { 2.0, 0.0 };
	vector <double> vec6 = { 2.0, 1.0 };
	Point<2,double> point1(vec1);
	Point<2,double> point2(vec2);
	Point<2,double> point3(vec3);
	Point<2,double> point4(vec4);
	Point<2,double> point5(vec5);
	Point<2,double> point6(vec6);
	Node<2, Point <2,double> > n_1(point1);
	Node<2, Point <2,double> > n_2(point2);
	Node<2, Point <2,double> > n_3(point3);
	vector<Node <2, Point <2,double>  > * > node_vec(3, nullptr);
	node_vec[0] = new Node<2, Point<2,double> >(point1);
	node_vec[1] = new Node<2, Point<2,double> >(point2);
	node_vec[2] = new Node<2, Point <2,double> >(point3);
	ElementFactory <2,3, Point<2,double> > factory;
	Element<2, 3, Point<2,double> > element = factory.build(node_vec);
	vector<Node <2, Point <2,double> > *> node_vec2;//(3, nullptr);
	node_vec2.push_back(node_vec[0]);
	node_vec2.push_back(node_vec[2]);
	node_vec2.push_back(new Node<2, Point <2,double> >(point4));
	Element<2, 3, Point<2,double> > el2 = factory.build(node_vec2);
	//Mesh<2, 3, Point <2,double> > el_mesh(el1);
	//el_mesh.push(el2);
	ElementDivider<2, 3, Point <2,double> > divider;
	element.set_indices(-1);
	el2.set_indices(2);
	el2.show();
	element.show();


	map< array<int, 2>, int> MIDPOINTS_MAP;
	int I = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			MIDPOINTS_MAP.insert(pair< array<int, 2>, int>({ i, j }, I));
			I++;
		}
	}

	typedef map<array<int, 2>, Point<2, double> >::const_iterator PointsMapIter;
	typedef map<array<int, 2>, Node<2, Point<2, double> >* >::const_iterator NodesMapIter;

	/*SECTION("dist_squared template function should work") {
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
			ElementDivider <2, 3, Point<2, double> > new_divider;
		}*/

	SECTION("Generating midpoint nodes when some are already generated should succeed") {
		//vector <Node <2, Point<2,double> > * > m_nods = element.get_midpoint_nodes();
		
		//vector<double> loc = m_nods[0]->get_location();
		//cout << loc[0] << endl;
		//Point<double> locations;
		
		//vector <Node<2, VectorXd> > xtr_nodes = factory.build_nodes(locations);
		//xtr_nodes[0].show();

		//Element<2, 3, VectorXd> inner_el = factory.build(m_nods);
		//inner_el.show();
		//typedef vector <Node<2, VectorXd>* > (Element <2, 3, VectorXd>::*ELPTR) ();
		//typedef vector <Node<2, Point<double> >* >(Element <2, 3, Point<double> >::*ELPTR) ();
		/*ELPTR el_ptr = &Element<2, 3, VectorXd>::get_nodes;
		vector <Node<2, VectorXd>* > new_nodes = (element.*el_ptr)();
		cout << "Hehe " << endl;
		new_nodes[0]->show();
		cout << new_nodes[0]->get_location();
		m_nods[0]->show();
		cout << "Hehe " << endl;
		cout << new_nodes.size() << endl;
		//for (int i = 0; i < m_nods.size(); i++) {
			//m_nods[i]->show();
		//}
		REQUIRE(m_nods.size() == 3);*/
	}

	SECTION("Generating map of midpoints should succeed") {
		cout << "El1" << endl;
		element.show();
		cout << "El2" << endl;
		el2.show();

		map<array<int, 2>, Point<2, double> > m_nods_map = element.get_midpoints_map();
		REQUIRE(m_nods_map.size() == 3);
		for (PointsMapIter iter = m_nods_map.begin(); iter != m_nods_map.end(); iter++) {
			REQUIRE(bool(iter->first[0] < iter->first[1]));//For this special case of indexing of nodes, not generally!!!
			REQUIRE(iter->second[0] >= 0);
			REQUIRE(iter->second[0] <= 1);
		}
	}

	SECTION("Generating map of midnodes should succeed") {//OK!
		map<array<int, 2>, Point<2, double> > points_map = element.get_midpoints_map();
		map<array<int, 2>, Node<2, Point<2, double> >* > commons;
		map<array<int, 2>, Node<2,Point<2, double> >* > m_p_node_map = divider.get_mid_nodes_map(element, commons);
		REQUIRE(m_p_node_map.size() == 3);
		for (NodesMapIter iter = m_p_node_map.begin(); iter != m_p_node_map.end(); iter++) {
			REQUIRE(bool(iter->first[0] < iter->first[1]));//For this special case of indexing of nodes, not generally!!!
			//iter->second->show();
			REQUIRE(iter->second->get_location() == points_map[iter->first]);
		}

		map<array<int, 2>, Point<2, double> > points_map2 = el2.get_midpoints_map();
		commons = m_p_node_map;
		map<array<int, 2>, Node<2, Point<2, double> >* > m_p_node_map2 = divider.get_mid_nodes_map(el2, commons);
		
		for (NodesMapIter iter = m_p_node_map2.begin(); iter != m_p_node_map2.end(); iter++) {
			iter->second->show();
			points_map2[iter->first].show();
			//REQUIRE(iter->second->get_location() == points_map2[iter->first]);
		}

	}

	/*SECTION("Generating new mid nodes from midpoints map should succeed") {
		map<array<int, 2>, Point<2, double> > m_p_map = element.get_midpoints_map();
		vector <Node<2, Point<2, double> >* > mid_nods = element.get_midpoint_nodes(m_p_map);
		//for (int i = 0; i < mid_nods.size(); i++) {
		//	mid_nods[i]->show();
		//}
		cout << "Showing midpoint of type Point <2, double> between nodes 0 and 1" << endl;
		m_p_map[{0, 1}].show();
		REQUIRE(mid_nods.size() == 3);
	}

	SECTION("Generating midpoint_map should succeed") {//OK!
		map<array<int, 2>, Point <2, double> > m_loc_map = divider.get_midlocation_map(element);
		cout << "Showing midlocations_map {0,1}" << endl;
		REQUIRE(m_loc_map.size() == 3);
		REQUIRE(m_loc_map[{0, 1}] == 0.5*(point1 + point2) );
		REQUIRE(m_loc_map[{0, 2}] == 0.5*(point1 + point3) );
		REQUIRE(m_loc_map[{1, 2}] == 0.5*(point2 + point3) );
		//REQUIRE(m_loc_map[{0, 2}] == point4)
		//REQUIRE(m_loc_map[{1, 2}] == point4)
	}


	SECTION("One can generate new vertex elements out of old element when refining the mesh") {
		vector <Node <2, Point <2, double> >* > mid_nodes = element.get_midpoint_nodes();
		Element <2, 3, Point <2, double> > el_AB = divider.get_vertex_element(0, mid_nodes, MIDPOINTS_MAP, element);
		REQUIRE(el_AB[0].get_location() == element[0].get_location());
		REQUIRE(el_AB.get_function(0)(el_AB[0].get_location()) == 1);
		REQUIRE(el_AB.get_function(1)(el_AB[0].get_location()) == 0);
		REQUIRE(el_AB[1].get_location() == 0.5*(element[0].get_location()+ element[1].get_location()));
		REQUIRE(el_AB[2].get_location() == 0.5*(element[0].get_location() + element[2].get_location()));
	}

	SECTION("One can generate new inner elements out of old element when refining the mesh") {
		vector <Node <2, Point <2, double> >* > m_nodes = element.get_midpoint_nodes();
		Element <2, 3, Point <2, double> > inner_el = divider.get_inner_element(0, m_nodes, MIDPOINTS_MAP);
		REQUIRE(inner_el.get_function(0)(inner_el[0].get_location()) == 1);
		REQUIRE(inner_el.get_function(1)(inner_el[0].get_location()) == 0);
		REQUIRE(inner_el[0].get_location() == 0.5*(element[0].get_location() + element[1].get_location()));
		REQUIRE(inner_el[1].get_location() == 0.5*(element[0].get_location() + element[2].get_location()));
		REQUIRE(inner_el[2].get_location() == 0.5*(element[1].get_location() + element[2].get_location()));
	}

	SECTION("Generating new Elements should succeed") {
		vector <Element <2, 3, Point <2,double> >* > els = divider.divide(element);
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
	}*/
}





/*TEST_CASE("Test ElementDivider with Node<2, VectorXd>") {

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

	ElementFactory<2, 3, VectorXd> factory;
	vector<Node <2, VectorXd> *> nodes2;//(3, nullptr);
	nodes2.push_back(nodes[0]);
	nodes2.push_back(nodes[2]);
	location << 0.0, 1.0;
	nodes2.push_back(new Node<2, VectorXd >(location));
	Element<2, 3, VectorXd> el2 = factory.build(nodes2);
	
	ElementDivider <2, 3, VectorXd> divider;

	map< array<int, 2>, int> MIDPOINTS_MAP;
	int I = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			MIDPOINTS_MAP.insert(pair< array<int, 2>, int>({ i, j }, I));
			I++;
		}
	}


	//cout << MIDPOINTS_MAP[{0, 1}] << endl;
	//cout << nodes[0]->get_location();

	//show_map(MIDPOINTS_MAP);

	/*SECTION("dist_squared template function should work") {
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

	SECTION("Generating midpoint nodes when some are already generated should succeed") {
		vector <Node<2, VectorXd>* > m_nods = divider.get_midpoint_nodes(element);
		vector<VectorXd> locations;
		VectorXd loc1(2);
		loc1 << 0.0, 0.0;
		locations.push_back(loc1);
		loc1 << 1.0, 0.0;
		locations.push_back(loc1);
		loc1 << 1.0, 1.0;
		locations.push_back(loc1);
		VectorXd temp(2);
		temp = locations[0];
		cout << temp;
		//vector <Node<2, VectorXd> > xtr_nodes = factory.build_nodes(locations);
		//xtr_nodes[0].show();
		
		//Element<2, 3, VectorXd> inner_el = factory.build(m_nods);
		//inner_el.show();
		//typedef vector <Node<2, VectorXd>* > (Element <2, 3, VectorXd>::*ELPTR) ();
		//typedef vector <Node<2, Point<double> >* >(Element <2, 3, Point<double> >::*ELPTR) ();
		/*ELPTR el_ptr = &Element<2, 3, VectorXd>::get_nodes;
		vector <Node<2, VectorXd>* > new_nodes = (element.*el_ptr)();
		cout << "Hehe " << endl;
		new_nodes[0]->show();
		cout << new_nodes[0]->get_location();
		m_nods[0]->show();
		cout << "Hehe " << endl;
		cout << new_nodes.size() << endl;
		//for (int i = 0; i < m_nods.size(); i++) {
			//m_nods[i]->show();
		//}
		REQUIRE(m_nods.size() == 3);*/
	//}

	/*SECTION("Generating map of midpoints should succeed") {
		typedef map<array<int, 2>, VectorXd>::const_iterator PointsMapIter;
		map<array<int, 2>, VectorXd> m_nods_map = element.get_midpoints_map();
		REQUIRE(m_nods_map.size() == 3);
		for (PointsMapIter iter = m_nods_map.begin(); iter != m_nods_map.end(); iter++) {
			REQUIRE(bool(iter->first[0] < iter->first[1]));//For this special case of indexing of nodes, not generally!!!
			REQUIRE(iter->second[0] >= 0);
			REQUIRE(iter->second[0] <= 1);
		}
	}

	SECTION("Generating new mid nodes from midpoints map should succeed") {
		map<array<int, 2>, VectorXd> m_p_map = element.get_midpoints_map();
		vector <Node<2, VectorXd>* > mid_nods = element.get_midpoint_nodes(m_p_map);
		for (int i = 0; i < mid_nods.size(); i++) {
			mid_nods[i]->show();
		}
		REQUIRE(mid_nods.size() == 3);
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


}*/

