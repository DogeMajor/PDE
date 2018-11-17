#include "../include/ElementDivider.h"
#include "../include/Point.h"
#include "../include/HelpfulTools.h"
#include "../include/Function.h"
#include "../include/TestingTools.h"

using namespace std;
using namespace Eigen;

/*
double fn(VectorXd coords) {
	return coords.transpose()*coords;
}

typedef map<array<int, 2>, Vertex<2, Point<2, double> >* >::const_iterator VerticesMapIter;
typedef map<array<int, 2>, Vertex<2, Point<2, double> >* > VerticesMap;

void show_map(VerticesMap &n_map) {
	for (VerticesMapIter iter = n_map.begin(); iter != n_map.end(); iter++) {
		cout << "Key and value" << endl;
		cout << iter->first[0] << ", " << iter->first[1] << endl;
		iter->second->show();
	}
}

typedef map<array<int, 2>, Point<2, double> >::const_iterator PointsMapIter;
typedef map<array<int, 2>, Point<2, double> > PointsMap;

void show_map(PointsMap &p_map) {
	for (PointsMapIter iter = p_map.begin(); iter != p_map.end(); iter++) {
		cout << "Key and value" << endl;
		cout << iter->first[0] << ", " << iter->first[1] << endl;
		cout << iter->second[0] << ", " << iter->second[1] << endl;
	}
}


bool point_bound_cond(Point<2, double> coords) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) { return true; }
	}
	return false;
}

bool point_bound_is_inside(Point<2, double> coords) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] <= 0.0) || (coords[i] >= 1.0)) { return false; }
	}
	return true;
}

double point_bound_val(Point<2, double> coords) {
	if (coords[1] == 1.0) { return 1; }
	return 0;
}
*/
#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"

TEST_CASE("Test ElementDivider with Point template -based vertices") {
	vector <double> vec1 = { 0.0, 0.0 };
	vector <double> vec2 = { 1.0, 0.0 };
	vector <double> vec3 = { 1.0, 1.0 };
	vector <double> vec4 = { 0.0, 1.0 };
	vector <double> vec5 = { 2.0, 0.0 };
	vector <double> vec6 = { 2.0, 1.0 };
	Point<2, double> point1(vec1);
	Point<2, double> point2(vec2);
	Point<2, double> point3(vec3);
	Point<2, double> point4(vec4);
	Point<2, double> point5(vec5);
	Point<2, double> point6(vec6);
	Vertex<2, Point <2, double> > n_1(point1);
	Vertex<2, Point <2, double> > n_2(point2);
	Vertex<2, Point <2, double> > n_3(point3);
	vector<Vertex<2, Point <2, double>  > * > vertex_vec(3, nullptr);
	vertex_vec[0] = new Vertex<2, Point<2, double> >(point1);
	vertex_vec[1] = new Vertex<2, Point<2, double> >(point2);
	vertex_vec[2] = new Vertex<2, Point <2, double> >(point3);
	ElementFactory <2, 3, Point<2, double> > factory;
	Element<2, 3, Point<2, double> > element = factory.build(vertex_vec);
	vector<Vertex<2, Point <2, double> > *> vertex_vec2;//(3, nullptr);
	vertex_vec2.push_back(vertex_vec[0]);
	vertex_vec2.push_back(vertex_vec[2]);
	vertex_vec2.push_back(new Vertex<2, Point <2, double> >(point4));
	Element<2, 3, Point<2, double> > el2 = factory.build(vertex_vec2);

	BoundaryConditions<Point<2, double> > boundaries;
	boundaries.cond_fn = point_bound_cond;
	boundaries.is_inside_fn = point_bound_is_inside;
	boundaries.val = point_bound_val;
	boundaries.accuracy = 0.001;

	ElementDivider<2, 3, Point <2, double> > divider(boundaries);
	element.set_indices(-1);
	element.set_index_maps();
	el2.set_indices(2);
	el2.set_index_maps();

	

	map< array<int, 2>, int> MIDPOINTS_MAP;
	int I = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			MIDPOINTS_MAP.insert(pair< array<int, 2>, int>({ i, j }, I));
			I++;
		}
	}


	SECTION("dist_squared template function should work") {
		VectorXd a(3);
		a << 0, 1, 2;
		VectorXd b(3);
		b << 1, 2, 5;
		REQUIRE(dist_squared<3, VectorXd>(a, b) == 11.0);
		REQUIRE(dist_squared<3, VectorXd>(a, a) == 0.0);
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

	SECTION("Test constructing ElementDivider") {
		ElementDivider <2, 3, Point<2, double> > new_divider(boundaries);
	}


	SECTION("Generating map of midpoints should succeed") {
		PointsMap m_nods_map = element.get_midpoints_map();
		REQUIRE(m_nods_map.size() == 3);
		for (PointsMapIter iter = m_nods_map.begin(); iter != m_nods_map.end(); iter++) {
			REQUIRE(bool(iter->first[0] < iter->first[1]));//For this special case of indexing of vertices, not generally!!!
			REQUIRE(iter->second[0] >= 0);
			REQUIRE(iter->second[0] <= 1);
		}
	}

	SECTION("Finding surfacepoint should succeed") {
		vector <double> old_vec = { 0.8, 0.8 };
		Point<2, double> old_loc(old_vec);
		Point<2, double> avg_loc = element.get_avg_location();
		Point<2, double> surface_loc = divider.find_surface_point(old_loc, avg_loc, 0.0001);
		REQUIRE(boundaries.cond(surface_loc) == true);
		vector <double> mp = { 0.5, 0.5 };
		REQUIRE(boundaries.cond(divider.find_surface_point(mp, avg_loc, 0.00001)) == true);
	}

	SECTION("Adjusting midpoints according to boundary conds should succeed "){
		PointsMap m_map = element.get_midpoints_map();
		PointsMap adjusted_midpoints = divider.adjust_midpoints(element, m_map, -1);
		show_map(adjusted_midpoints);
	}

	SECTION("Generating map of midvertices should succeed") {//OK!
		map<array<int, 2>, Point<2, double> > points_map = element.get_midpoints_map();
		int vertices_no = n_1.how_many();
		map<array<int, 2>, Vertex<2, Point<2, double> >* > commons;
		map<array<int, 2>, Vertex<2,Point<2, double> >* > m_p_vertex_map = divider.get_mid_vertices_map(element, commons, 8999);//Testing with "it's almost 9000 meme"
		REQUIRE(m_p_vertex_map.size() == 3 );
		REQUIRE(vertices_no + 3 == n_1.how_many() );
		for (VerticesMapIter iter = m_p_vertex_map.begin(); iter != m_p_vertex_map.end(); iter++) {
			REQUIRE(bool(iter->first[0] < iter->first[1]) );//For this special case of indexing of vertices, not generally!!!
			REQUIRE(iter->second->get_location() == points_map[iter->first] );
		}
		show_map(m_p_vertex_map);
		REQUIRE(commons == m_p_vertex_map );//First element, at t_0 commons is empty!!
		map<array<int, 2>, Point<2, double> > points_map2 = el2.get_midpoints_map();
		map<array<int, 2>, Vertex<2, Point<2, double> >* > m_p_vertex_map2 = divider.get_mid_vertices_map(el2, commons, 8999);
		REQUIRE(vertices_no + 5 == n_1.how_many());
		for (VerticesMapIter iter = m_p_vertex_map2.begin(); iter != m_p_vertex_map2.end(); iter++) {
			REQUIRE( iter->second->get_location() == points_map2[iter->first] );
		}
		for (VerticesMapIter iter = m_p_vertex_map2.begin(); iter != m_p_vertex_map2.end(); iter++) {
			REQUIRE( iter->second == commons[iter->first] );
		}
		REQUIRE(commons.size() == 5 );
	}

	SECTION("Generating new mid vertices from midpoints map should succeed") {
		map<array<int, 2>, Point<2, double> > m_p_map = element.get_midpoints_map();
		vector <Vertex<2, Point<2, double> >* > mid_nods = element.get_midpoint_vertices(m_p_map);

		cout << "Showing midpoint of type Point <2, double> between vertices 0 and 1" << endl;
		m_p_map[{0, 1}].show();
		REQUIRE(mid_nods.size() == 3 );
	}

	SECTION("One can generate new vertex elements out of old element when refining the mesh") {
		vector <Vertex<2, Point <2, double> >* > mid_vertices = element.get_midpoint_vertices();
		Element <2, 3, Point <2, double> > el_AB = divider.get_corner_element(0, mid_vertices, element);
		REQUIRE(el_AB[0].get_location() == element[0].get_location());
		REQUIRE(el_AB.get_function(0)(el_AB[0].get_location()) == 1);
		REQUIRE(el_AB.get_function(1)(el_AB[0].get_location()) == 0);
		REQUIRE(el_AB[1].get_location() == 0.5*(element[0].get_location()+ element[1].get_location()));
		REQUIRE(el_AB[2].get_location() == 0.5*(element[0].get_location() + element[2].get_location()));
	}

	SECTION("One can generate new inner elements out of old element when refining the mesh") {
		vector <Vertex<2, Point <2, double> >* > m_vertices = element.get_midpoint_vertices();
		Element <2, 3, Point <2, double> > inner_el = divider.get_inner_element(0, m_vertices);
		REQUIRE(inner_el.get_function(0)(inner_el[0].get_location()) == 1);
		REQUIRE(inner_el.get_function(1)(inner_el[0].get_location()) == 0);
		REQUIRE(inner_el[0].get_location() == 0.5*(element[0].get_location() + element[1].get_location()));
		REQUIRE(inner_el[1].get_location() == 0.5*(element[0].get_location() + element[2].get_location()));
		REQUIRE(inner_el[2].get_location() == 0.5*(element[1].get_location() + element[2].get_location()));
	}

	SECTION( "Generating new Elements should succeed" ) {
		map< array<int, 2>, Vertex<2, Point <2, double> >* >  coms;
		int vertices_at_t0 = n_1.how_many();
		
		vector <Element <2, 3, Point <2,double> >* > els = divider.divide(element, coms, 8999);
		
		REQUIRE(els.size() == 4);
		REQUIRE(vertices_at_t0 + 3 == n_1.how_many() );
		vector <Element <2, 3, Point <2, double> >* > els2 = divider.divide(el2, coms, 8999);
		REQUIRE(els2.size() == 4);
		REQUIRE(vertices_at_t0 + 5 == n_1.how_many() );
		coms.erase(coms.begin(), coms.begin()); 
	}

}





/*TEST_CASE("Test ElementDivider with Vertex<2, VectorXd>") {

	VectorXd location(2);
	location << 0.0, 0.0;
	Vertex<2, VectorXd> vertex1(location);
	location << 1.0, 0.0;
	Vertex<2, VectorXd> vertex2(location);
	location << 1.0, 1.0;
	Vertex<2, VectorXd> vertex3(location);
	vector<Vertex<2, VectorXd> *> vertices(3, nullptr);
	location << 0.0, 0.0;
	vertices[0] = new Vertex<2, VectorXd>(location);
	location << 1.0, 0.0;
	vertices[1] = new Vertex<2, VectorXd>(location);
	location << 1.0, 1.0;
	vertices[2] = new Vertex<2, VectorXd>(location);
	//cout << vertices[1]->how_many() << endl;
	vector<SimplexFunction <VectorXd> > funcs(3);
	VectorXd coeffs(3);
	coeffs << -1, 0, 1;
	funcs[0].coeff = coeffs;
	coeffs << 1, -1, 0;
	funcs[1].coeff = coeffs;
	coeffs << 0, 1, 0;
	funcs[2].coeff = coeffs;
	Element <2, 3, VectorXd> element(vertices, funcs);

	ElementFactory<2, 3, VectorXd> factory;
	vector<Vertex<2, VectorXd> *> vertices2;//(3, nullptr);
	vertices2.push_back(vertices[0]);
	vertices2.push_back(vertices[2]);
	location << 0.0, 1.0;
	vertices2.push_back(new Vertex<2, VectorXd >(location));
	Element<2, 3, VectorXd> el2 = factory.build(vertices2);
	
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
	//cout << vertices[0]->get_location();

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

	/*SECTION("Generating map of midpoints should succeed") {
		typedef map<array<int, 2>, VectorXd>::const_iterator PointsMapIter;
		map<array<int, 2>, VectorXd> m_nods_map = element.get_midpoints_map();
		REQUIRE(m_nods_map.size() == 3);
		for (PointsMapIter iter = m_nods_map.begin(); iter != m_nods_map.end(); iter++) {
			REQUIRE(bool(iter->first[0] < iter->first[1]));//For this special case of indexing of vertices, not generally!!!
			REQUIRE(iter->second[0] >= 0);
			REQUIRE(iter->second[0] <= 1);
		}
	}


	//SECTION("Generating location_map should succeed") {
		//map<VectorXd, array<int, 2> > location_map = divider.get_locations_map(element);
		//cout << location_map.size() << endl;
		//REQUIRE(m_map == MIDPOINTS_MAP);
	//}

	SECTION("Calculating average location should succeed") {
		vector <Vertex<2, VectorXd >* > el_vertices = element.get_vertices();
		VectorXd avg_loc = divider.average_location(el_vertices);
		REQUIRE(limit_decimals(avg_loc[0],4) == 0.6666);
		REQUIRE(limit_decimals(avg_loc[1],4) == 0.3333);
	}

	SECTION("Finding nearest vertex should succeed") {
		vector <Vertex<2, VectorXd >* > el_vertices = element.get_vertices();
		VectorXd z(2);
		z << 0.9, 0.9;
		Vertex<2, VectorXd > nearest_vertex = divider.nearest_vertex(z, el_vertices);
		REQUIRE(nearest_vertex.get_location() == el_vertices[2]->get_location());
	}


}*/

