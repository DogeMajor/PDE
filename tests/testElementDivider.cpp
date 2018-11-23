#include "../include/ElementDivider.h"
#include "../include/Point.h"
#include "../include/HelpfulTools.h"
#include "../include/Function.h"
#include "../include/TestingTools.h"

using namespace std;
using namespace Eigen;


#define CATCH_CONFIG_MAIN
#include "../C++ libs/catch/catch.hpp"

/*TEST_CASE("Test ElementDivider with Point template -based vertices") {
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

	vector <double> midvec12 = { 0.5, 0.0 };
	vector <double> midvec13 = { 0.5, 0.5 };
	vector <double> midvec23 = { 1.0, 0.5 };
	Point<2, double> midpoint12(midvec12);
	Point<2, double> midpoint13(midvec13);
	Point<2, double> midpoint23(midvec23);
	vector<Vertex<2, Point <2, double>  > * > midpoint_vertices(3, nullptr);
	midpoint_vertices[0] = new Vertex<2, Point<2, double> >(midpoint12);
	midpoint_vertices[1] = new Vertex<2, Point<2, double> >(midpoint13);
	midpoint_vertices[2] = new Vertex<2, Point<2, double> >(midpoint23);

	Vertex<2, Point <2, double> > n_1(point1);
	Vertex<2, Point <2, double> > n_2(point2);
	Vertex<2, Point <2, double> > n_3(point3);
	vector<Vertex<2, Point <2, double>  > * > vertex_vec(3, nullptr);
	vertex_vec[0] = new Vertex<2, Point<2, double> >(point1);
	vertex_vec[1] = new Vertex<2, Point<2, double> >(point2);
	vertex_vec[2] = new Vertex<2, Point <2, double> >(point3);
	ElementFactory <2, 3, Point<2, double> > factory;
	Element<2, 3, Point<2, double> > element = factory.build(vertex_vec);
	vector<Vertex<2, Point <2, double> > *> vertex_vec2;
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

	map<array<int, 2>, int> empty_edges;

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
		PointsMap adjusted_midpoints = divider.adjust_midpoints(element, m_map, empty_edges);
		cout << "Show adjusted mid points" << endl;
		show_map(adjusted_midpoints);
		REQUIRE(adjusted_midpoints[{0, 1}] == m_map[{0, 1}]);
		REQUIRE(adjusted_midpoints[{1, 2}] == m_map[{1, 2}]);
		
	}

	SECTION("Generating map of midvertices should succeed") {//OK!
		map<array<int, 2>, Point<2, double> > points_map = element.get_midpoints_map();
		int vertices_no = n_1.how_many();
		map<array<int, 2>, Vertex<2, Point<2, double> >* > commons;
		map<array<int, 2>, Vertex<2,Point<2, double> >* > m_p_vertex_map = divider.get_mid_vertices_map(element, commons, empty_edges);//Testing with "it's almost 9000 meme"
		REQUIRE(m_p_vertex_map.size() == 3 );
		REQUIRE(vertices_no + 3 == n_1.how_many() );
		for (VerticesMapIter iter = m_p_vertex_map.begin(); iter != m_p_vertex_map.end(); iter++) {
			REQUIRE(bool(iter->first[0] < iter->first[1]) );//For this special case of indexing of vertices, not generally!!!
			REQUIRE(iter->second->get_location() == points_map[iter->first] );
		}
		show_map(m_p_vertex_map);
		REQUIRE(commons == m_p_vertex_map );//First element, at t_0 commons is empty!!
		map<array<int, 2>, Point<2, double> > points_map2 = el2.get_midpoints_map();
		map<array<int, 2>, Vertex<2, Point<2, double> >* > m_p_vertex_map2 = divider.get_mid_vertices_map(el2, commons, empty_edges);
		REQUIRE(vertices_no + 5 == n_1.how_many());
		for (VerticesMapIter iter = m_p_vertex_map2.begin(); iter != m_p_vertex_map2.end(); iter++) {
			REQUIRE( iter->second->get_location() == points_map2[iter->first] );
		}
		for (VerticesMapIter iter = m_p_vertex_map2.begin(); iter != m_p_vertex_map2.end(); iter++) {
			REQUIRE( iter->second == commons[iter->first] );
		}
		REQUIRE(commons.size() == 5 );
	}

	SECTION("One can generate new vertex elements out of old element when refining the mesh") {
		Element <2, 3, Point <2, double> > el_AB = divider.get_corner_element(0, midpoint_vertices, element);
		REQUIRE(el_AB[0].get_location() == element[0].get_location());
		REQUIRE(el_AB.get_function(0)(el_AB[0].get_location()) == 1);
		REQUIRE(el_AB.get_function(1)(el_AB[0].get_location()) == 0);
		REQUIRE(el_AB[1].get_location() == 0.5*(element[0].get_location()+ element[1].get_location()));
		REQUIRE(el_AB[2].get_location() == 0.5*(element[0].get_location() + element[2].get_location()));
	}

	SECTION("One can generate new inner elements out of old element when refining the mesh") {
		Element <2, 3, Point <2, double> > inner_el = divider.get_inner_element(0, midpoint_vertices);
		REQUIRE(inner_el.get_function(0)(inner_el[0].get_location()) == 1);
		REQUIRE(inner_el.get_function(1)(inner_el[0].get_location()) == 0);
		REQUIRE(inner_el[0].get_location() == 0.5*(element[0].get_location() + element[1].get_location()));
		REQUIRE(inner_el[1].get_location() == 0.5*(element[0].get_location() + element[2].get_location()));
		REQUIRE(inner_el[2].get_location() == 0.5*(element[1].get_location() + element[2].get_location()));
	}

	SECTION( "Generating new Elements without boundary adjustments should succeed" ) {
		map< array<int, 2>, Vertex<2, Point <2, double> >* >  coms;
		int vertices_at_t0 = n_1.how_many();
		
		vector <Element <2, 3, Point <2,double> >* > els = divider.divide(element, coms, empty_edges);
		
		REQUIRE(els.size() == 4);
		REQUIRE(vertices_at_t0 + 3 == n_1.how_many() );
		empty_edges.erase(empty_edges.begin(), empty_edges.end());
		vector <Element <2, 3, Point <2, double> >* > els2 = divider.divide(el2, coms, empty_edges);
		REQUIRE(els2.size() == 4);
		REQUIRE(vertices_at_t0 + 5 == n_1.how_many() );
		coms.erase(coms.begin(), coms.end()); 
	}

	SECTION("Generating new Elements with boundary adjustments should succeed") {
		map< array<int, 2>, Vertex<2, Point <2, double> >* >  adj_coms;
		int vertices_amount = n_1.how_many();
		cout << "Empty edges size" << empty_edges.size() << endl;
		vector <Element <2, 3, Point <2, double> >* > adjusted_els = divider.divide(element, adj_coms, empty_edges);

		REQUIRE(adjusted_els.size() == 4);
		REQUIRE(adjusted_els.size() == 4);
		REQUIRE(vertices_amount + 3 == n_1.how_many());
		double volume = 0;
		for (int i = 0; i < adjusted_els.size(); i++) {
			volume += adjusted_els[i]->get_volume();
		}
		REQUIRE(abs(volume - 0.5) < 0.01);
		
	}

}*/


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

TEST_CASE("Dividing elements in 3-D Element<3,4,Vertex<2, VectorXd> > should succeed") {

	VectorXd p_loc(3);
	p_loc << 0, 0, 0;
	Vertex<3, VectorXd> v1(p_loc);
	p_loc << 1, 0, 0;
	Vertex<3, VectorXd> v2(p_loc);
	p_loc << 0, 1, 0;
	Vertex<3, VectorXd> v3(p_loc);
	p_loc << 0, 0, 1;
	Vertex<3, VectorXd> v4(p_loc);
	vector<Vertex<3, VectorXd> *> p_vertices(4, nullptr);
	p_loc << 0, 0, 0;
	p_vertices[0] = new Vertex<3, VectorXd>(p_loc);
	p_loc << 1, 0, 0;
	p_vertices[1] = new Vertex<3, VectorXd>(p_loc);
	p_loc << 0, 1, 0;
	p_vertices[2] = new Vertex<3, VectorXd>(p_loc);
	p_loc << 0, 0, 1;
	p_vertices[3] = new Vertex<3, VectorXd>(p_loc);



	vector<Vertex<3, VectorXd> *> p_mid_vertices(6, nullptr);//6 = sum(i=0, i=n) (i)  when n = 3
	p_loc << 0.5, 0, 0;
	p_mid_vertices[0] = new Vertex<3, VectorXd>(p_loc);
	p_loc << 0, 0.5, 0;
	p_mid_vertices[1] = new Vertex<3, VectorXd>(p_loc);
	p_loc << 0, 0, 0.5;
	p_mid_vertices[2] = new Vertex<3, VectorXd>(p_loc);


	p_loc << 0.5, 0.5, 0;
	p_mid_vertices[3] = new Vertex<3, VectorXd>(p_loc);
	p_loc << 0.5, 0, 0.5;
	p_mid_vertices[4] = new Vertex<3, VectorXd>(p_loc);

	p_loc << 0, 0.5, 0.5;
	p_mid_vertices[5] = new Vertex<3, VectorXd>(p_loc);

	ElementFactory<3, 4, VectorXd> p_factory;
	Element<3, 4, VectorXd> p_element = p_factory.build(p_vertices);
	p_element.set_indices(-1);
	p_element.set_index_maps();

	vector<Vertex<3, VectorXd> *> p_vertices2(4, nullptr);
	
	p_vertices2[0] = p_vertices[0];
	p_vertices2[3] = p_vertices[3];
	p_loc << 1, 0.5, 0.5;
	p_vertices2[1] = new Vertex<3, VectorXd>(p_loc);
	p_loc << 0.5, 1, 0.5;
	p_vertices2[2] = new Vertex<3, VectorXd>(p_loc);
	Element<3, 4, VectorXd> p_element2 = p_factory.build(p_vertices2);
	p_element2.set_indices(3);
	p_element2.set_index_maps();

	BoundaryConditions<VectorXd> p_boundaries = { bound_cond, bound_is_inside, bound_val, bound_normal, 0.0001 };
	ElementDivider <3, 4, VectorXd> p_divider(p_boundaries);

	map< array<int, 2>, int> P_EDGES_MAP;
	int A, B;
	int I = 0;
	for (int i = 0; i < 4; i++) {
		A = p_element2[i].get_index();
		for (int j = i + 1; j < 4; j++) {
			B = p_element2[j].get_index();
			P_EDGES_MAP.insert(pair< array<int, 2>, int>({ min(A,B), max(A,B) }, 1));
			I++;
		}
	}
	map<array<int, 2>, int> p_empty_edges;

	/*SECTION("Generating corner element should succeed") {
		Element<3, 4, VectorXd> corner0 = p_divider.get_corner_element(0, p_mid_vertices, p_element);
		REQUIRE(corner0[0] == *p_vertices[0]);
		REQUIRE(corner0[1] == *p_mid_vertices[0]);
		REQUIRE(corner0[2] == *p_mid_vertices[1]);
		REQUIRE(corner0[3] == *p_mid_vertices[2]);

		Element<3, 4, VectorXd> corner3 = p_divider.get_corner_element(3, p_mid_vertices, p_element);
		//corner3.show();
		cout << corner3.get_volume() << endl;
		cout << p_element.get_volume() << endl;
	}

	SECTION("Generating inner element should succeed") {
		Element<3, 4, VectorXd> inner0 = p_divider.get_inner_element(0, p_mid_vertices);
		//REQUIRE(inner0[0] == *p_vertices[0]);
		Element<3, 4, VectorXd> inner2 = p_divider.get_inner_element(2, p_mid_vertices);
		//REQUIRE(inner0[0] == *p_vertices[0]);
		REQUIRE(inner0[3] == inner2[3]);
	}

	SECTION("Generating mid vertices should succeed") {
		map<array<int, 2>, Vertex<3, VectorXd>* > p_common_vert;
		map <array<int, 2>, Vertex<3, VectorXd>* > mid_verts = p_divider.get_mid_vertices_map(p_element, p_common_vert, p_empty_edges);

		//REQUIRE(inner0[0] == *p_vertices[0]);
		cout << mid_verts.size() << endl;
		for (auto iter = mid_verts.begin(); iter != mid_verts.end(); iter++) {
			cout << iter->second->get_location() << endl;
			cout << endl;
		}
	}

	SECTION("Dividing element should succeed") {
		map<array<int, 2>, Vertex<3, VectorXd>* > p_comms;
		vector<Element<3, 4, VectorXd>* > p_els = p_divider.divide(p_element, p_comms, p_empty_edges);
		
		//REQUIRE(inner0[0] == *p_vertices[0]);
		double vol = 0;
		VectorXd avg(3);
		avg << 0, 0, 0;
		for (int i = 0; i < p_els.size(); i++) {
			vol += p_els[i]->get_volume();
			avg = avg + p_els[i]->get_avg_location()*p_els[i]->get_volume();
			//p_els[i]->show();
			cout << p_els[i]->get_avg_location() << endl;
			cout << endl;
		}
		avg = avg * (1 / vol);
		REQUIRE(limit_decimals(p_element.get_volume(), 4) == 0.1666);
		REQUIRE(limit_decimals(vol, 4) == 0.1666);
		REQUIRE(avg == p_element.get_avg_location());

	}*/

	SECTION("Dividing element with adjustments towards the boundary should succeed") {
		map < array<int, 2>, int> edges;
		map<array<int, 2>, Vertex<3, VectorXd>* > p_comm;
		vector<Element<3, 4, VectorXd>* > adj_p_els = p_divider.divide(p_element2, p_comm, P_EDGES_MAP);

		//REQUIRE(inner0[0] == *p_vertices[0]);
		double vol = 0;
		VectorXd avg(3);
		avg << 0, 0, 0;
		for (int i = 0; i < adj_p_els.size(); i++) {
			vol += adj_p_els[i]->get_volume();
			avg = avg + adj_p_els[i]->get_avg_location()*adj_p_els[i]->get_volume();
			//p_els[i]->show();
			cout << (*adj_p_els[i])[3].get_location() << endl;
			cout << p_boundaries.cond((*adj_p_els[i])[3].get_location()) << endl;
			cout << endl;
			//REQUIRE(p_boundaries.cond((*adj_p_els[i])[3].get_location()) == true);
		}
		avg = avg * (1 / vol);
		adj_p_els[2]->show();

	}

	

}