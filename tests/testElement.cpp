#include "../include/point.h"
#include "../include/node.h"
#include "../include/element.h"
#include <math.h>
#include <tuple>
#include <array>
#include <map>
#include <memory>

using namespace std;
using namespace Eigen;


#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"

/*
TEST_CASE( "Test Element template containing Node template initiated with 2-D double vector from Eigen lib" ) {

    VectorXd location(2);
    location << 0.0, 0.0;
    //cout << "VectorXd size: " << location.size() << endl;
    Node <2,VectorXd> node1(location);
    location << 1.0, 0.0;
    Node <2,VectorXd> node2(location);
    location << 1.0, 1.0;
    Node <2,VectorXd> node3(location);
    Node <2,VectorXd> * nodes[3];
    nodes[0] = &node1;
    nodes[1] = &node2;
    nodes[2] = &node3;
	cout << nodes[1]->get_location() << endl;
    vector <SimplexFunction <VectorXd> > funcs(3);
    VectorXd coeffs(3);
    coeffs << -1,0,1;
    funcs[0].coeff = coeffs;
    coeffs << 1,-1,0;
    funcs[1].coeff = coeffs;
    coeffs << 0,1,0;
    funcs[2].coeff = coeffs;
	node1.show();
    Element <2, 3, VectorXd> element(nodes, funcs);

	
    SECTION( "Test default constructor()" ){
        Element <2, 3, VectorXd> empty_element;
        REQUIRE( empty_element[0].get_shared_elements() == 1 );
        REQUIRE( empty_element[0].get_index() == -1 );
        REQUIRE( empty_element[2].get_shared_elements() == 1 );
        REQUIRE( empty_element[2].get_index() == -1 );
    }

	SECTION("Test get_nodes") {
		Node <2, VectorXd> **retrieved_nodes = element.get_nodes();
		retrieved_nodes[0]->show();
		retrieved_nodes[1]->show();
		retrieved_nodes[2]->show();
		REQUIRE( retrieved_nodes[2]->get_node_amount() == 6 );
		REQUIRE( retrieved_nodes[2]->get_location() == node3.get_location() );
		REQUIRE( retrieved_nodes[2]->get_index() == node3.get_index() );
		REQUIRE( retrieved_nodes[2]->get_shared_elements() == node3.get_shared_elements()+1 );
		//Node <2, VectorXd> nod_2 = *retrieved_nodes[1];
		//cout << nod_2.get_location() << endl;
		//REQUIRE(retrieved_nodes[2]->get_location() == location);
	}

    SECTION( "Test operator []" ){
        REQUIRE( element[0].get_location() == nodes[0]->get_location() );
        REQUIRE( element[0].get_location() != nodes[1]->get_location() );
    }

    SECTION( "Test show()" ){
        cout << "showing element[0]" << endl;
        element[0].show();
        cout << "showing node1" << endl;
        nodes[0]->show();
    }

    SECTION( "Test copy constructor" ){
        Element <2, 3, VectorXd> copyed_element(element);
        REQUIRE( copyed_element[0] == element[0] );
        REQUIRE( copyed_element[1] == element[1] );
        REQUIRE( copyed_element[2] == element[2] );
    }

    SECTION( "Test increase_shared_elements" ){
        int shared_els = element[0].get_shared_elements();
        element.increase_shared_elements();
        REQUIRE( shared_els+1 == element[0].get_shared_elements() );
    }

    SECTION( "Test set_indices" ){
        element.set_indices();
        REQUIRE( 0 == element[0].get_index() );
        REQUIRE( 1 == element[1].get_index() );
    }

    SECTION( "Test get_function(int)" ){
        SimplexFunction<VectorXd> first_fn = element.get_function(0);
        REQUIRE( first_fn(node1.get_location()) == 1 );
        REQUIRE( first_fn(node2.get_location()) == 0 );
        REQUIRE( first_fn(node3.get_location()) == 0 );
    }

    SECTION( "Test assignment operator" ){
        Element <2, 3, VectorXd> assigned_element = element;
        REQUIRE( assigned_element[0] == element[0] );
        REQUIRE( assigned_element[1] == element[1] );
        REQUIRE( assigned_element[2] == element[2] );
    }


    SECTION( "Test operators == and !=" ){
        Element <2, 3, VectorXd> new_element(nodes, funcs);
        REQUIRE( new_element == element );
        Node <2,VectorXd> *reflected_nodes[3];
        reflected_nodes[0] = &node3;
        reflected_nodes[1] = &node2;
        reflected_nodes[2] = &node1;
        Element <2, 3, VectorXd> reflected_element(reflected_nodes, funcs);
        REQUIRE( reflected_element != element );
    }

    SECTION( "Test get_volume()" ){
        REQUIRE(element.get_volume() == 0.5 );
    }

	SECTION("Test midpoints()") {
		vector <pair <int[2], VectorXd> > midpoints = element.get_midpoints();
		VectorXd mid_loc(2);
		mid_loc << 0.5, 0.0;
		REQUIRE( midpoints[0].second == mid_loc);
		mid_loc << 1, 0.5;
		REQUIRE(midpoints[2].second == mid_loc);
		REQUIRE(midpoints[0].first[0] == 0);
		REQUIRE(midpoints[0].first[1] == 1);
		REQUIRE(midpoints[2].first[0] == 1);
		REQUIRE(midpoints[2].first[1] == 2);
		REQUIRE(midpoints.size() == 3);
	}

	
	SECTION("Test midpoint_nodes()") {
		VectorXd loc(2);
		loc << 0.5, 0.0;
		int node_no = element[0].get_node_amount();
		vector <pair <int[2], VectorXd> > mpoints = element.get_midpoints();
		Node<2, VectorXd> ** new_nodes = element.get_midpoint_nodes(mpoints);
		new_nodes[2]->show();
		REQUIRE( element[0].get_node_amount() == node_no+3 );
		//REQUIRE(new_nodes[0]->get_location() == loc);
		//loc << 0.5, 0.5;
		//REQUIRE(new_nodes[1]->get_location() == loc);
		//loc << 0.5, 1.0;
		//REQUIRE(new_nodes[2]->get_location() == loc);
		//REQUIRE(midnodes[0].first[1] == 1);
	}

	SECTION("Test vector of elements") {
		Element <2, 3, VectorXd> e1(nodes, funcs);
		Element <2, 3, VectorXd> e2(nodes, funcs);
		vector < Element <2, 3, VectorXd> > els;
		els.push_back(e1);
		els.push_back(e2);
		els[0].show();
	}
	 
	SECTION( "Testing getting new nodes()" ) {
		map< array<int, 2>, Node<2, VectorXd> * > new_nodes = element.get_new_nodes();
		//for (int i = 0; i < 3; i++) {
			//for (int j = i + 1; j < 3; j++) {
				//REQUIRE(*new_nodes[{i, j}] == element[i + j - 1]);
				//cout << new_nodes[{i, j}]->get_location() << endl;
			//}
		//}
		cout << "size:" << new_nodes.size();
		//typedef std::shared_ptr< Node<2, VectorXd> > NodePtr;
		//shared_ptr< NodePtr > ptr = make_shared<Node<2, VectorXd> >();
		//Node<2, VectorXd> *ptr;
		//ptr = new_nodes[{0, 0}];
		//ptr->show();
		//VectorXd loc_00 = new_nodes[{0, 0}]->get_location();
		//cout << loc_00 << endl;
	}
	
/*
	SECTION("Test map of nodes") {
		vector <pair <int[2], VectorXd> > m_points = element.get_midpoints();
		Node<2, VectorXd> ** m_nodes = element.get_midpoint_nodes(m_points);
		map< array<int, 2>, Node<2, VectorXd>* > node_map;
		//Node<2, VectorXd> n6(*m_nodes[1]);
		//cout << n6.get_location() << endl;
		for (int i = 0; i < 3; i++) {
			for (int j = i + 1; j < 3; j++) {
				node_map[{i, j}] = new Node<2, VectorXd>(*m_nodes[i+j-1]);
			}
		}
		node_map[{0, 1}]->show();
		//cout << node_map[{0, 1}]->get_location() << endl;
	}

}*/





TEST_CASE( "Test Element template containing Node template initiated with 2-D Point <double> template objects" ) {

    vector <double> vec1 = {0.0, 0.0};
    vector <double> vec2 = {1.0, 0.0};
    vector <double> vec3 = {1.0, 1.0};
    Point <double> point1(vec1);
    Point <double> point2(vec2);
    Point <double> point3(vec3);
    Node <2,Point <double> > n_1(point1);
    Node <2,Point <double> > n_2(point2);
    Node <2,Point <double> > n_3(point3);

	vector < Node <2, Point <double> >* > node_vec;
	cout << n_1.how_many() << endl;
	node_vec.push_back(&n_1);
	node_vec.push_back(&n_2);
	node_vec.push_back(&n_3);

    vector <SimplexFunction <Point <double> > > fns(3);
    VectorXd coeff(3);
    coeff << -1,0,1;
    fns[0].coeff = coeff;
    coeff << 1,-1,0;
    fns[1].coeff = coeff;
    coeff << 0,1,0;
    fns[2].coeff = coeff;
    Element <2, 3, Point <double> > el(node_vec, fns);

	map< array<int, 2>, int> MIDPOINTS_MAP;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			MIDPOINTS_MAP.insert(pair< array<int, 2>, int>({ i, j }, i*3 + j - 1));
		}
	}

	SECTION("Test operator []") {
		REQUIRE(el[0].get_location() == node_vec[0]->get_location());
		REQUIRE(el[0].get_location() != node_vec[1]->get_location());
		REQUIRE(el[0].get_location()[0] == 0.0);
		REQUIRE(el[0].get_location()[1] == 0.0);
	}

    SECTION( "Test show()" ){
        cout << "showing element[0]" << endl;
        el[0].show();
        cout << "showing node1" << endl;
        node_vec[0]->show();
    }

    SECTION( "Test get_simplex_matrix(T &el)" ){
        Matrix<double, 2, 2> s_mat = el.get_simplex_matrix(el);
        REQUIRE( s_mat == MatrixXd::Identity(2,2) );
    }

    SECTION( "Test get_volume()" ){
        REQUIRE( el.get_volume() == 0.5 );
    }

    SECTION( "Test copy constructor" ){
        Element <2, 3, Point <double> > copyed_el(el);
        REQUIRE( copyed_el[0] == el[0] );
        REQUIRE( copyed_el[1] == el[1] );
        REQUIRE( copyed_el[2] == el[2] );
    }

    SECTION( "Test increase_shared_elements" ){
        int shared_els = el[0].get_shared_elements();
        el.increase_shared_elements();
        REQUIRE( shared_els+1 == el[0].get_shared_elements() );
    }

    SECTION( "Test set_indices" ){
        el.set_indices();
        REQUIRE( 0 == el[0].get_index() );
        REQUIRE( 1 == el[1].get_index() );
    }

    SECTION( "Test assignment operator" ){
        Element <2, 3, Point <double> > assigned_el = el;
        REQUIRE( assigned_el[0] == el[0] );
        REQUIRE( assigned_el[1] == el[1] );
        REQUIRE( assigned_el[2] == el[2] );
        REQUIRE( assigned_el.get_function(0) == el.get_function(0) );
    }

    SECTION( "Test get_function(int)" ){
        SimplexFunction< Point <double> > func_1 = el.get_function(0);
        REQUIRE( func_1(n_1.get_location()) == 1 );
        REQUIRE( func_1(n_2.get_location()) == 0 );
        REQUIRE( func_1(n_3.get_location()) == 0 );
    }
	
	SECTION("Test midpoint_nodes()") {
		VectorXd loc(2);
		loc << 0.5, 0.0;
		int node_no = el[0].how_many();
		vector <pair <int[2], Point <double> > > m_points = el.get_midpoints();
		vector <Node <2, Point <double> >* > m_nodes = el.get_midpoint_nodes();
		m_nodes[2]->show();
		REQUIRE(el[0].how_many() == node_no + 3);
		REQUIRE(m_nodes[0]->get_location()[0] == 0.5*(vec1[0]+ vec2[0]));
		REQUIRE(m_nodes[1]->get_location()[1] == 0.5*(vec1[1] + vec3[1]));
		REQUIRE(m_nodes[2]->get_location()[0] == 0.5*(vec2[0] + vec3[0]));
	}


	SECTION("One can generate new vertex elements out of old element while refining the mesh") {
		vector <Node <2, Point <double> >* > mid_nodes = el.get_midpoint_nodes();
		Element <2, 3, Point <double> > el_AB = el.get_vertex_element(0, mid_nodes, MIDPOINTS_MAP);
		el_AB.show();
	}
}
