#include "../include/point.h"
#include "../include/node.h"
#include "../include/element.h"

using namespace std;
using namespace Eigen;

//N-dim box's boundary
bool bound_cond(VectorXd coords) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) { return true; }
	}
	return false;
}

bool bound_cond_large_box(VectorXd coords) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 2.0)) { return true; }
	}
	return false;
}

double bound_val(VectorXd coords) {
	if (coords[0] == 1.0) { return 1; }
	return 0;
}


#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


/*TEST_CASE( "Test Element template containing Node template initiated with 2-D double vector from Eigen lib" ){

    VectorXd location(2);
    location << 0.0, 0.0;
    Node <2,VectorXd> node1(location);
    location << 1.0, 0.0;
    Node <2,VectorXd> node2(location);
    location << 1.0, 1.0;
    Node <2,VectorXd> node3(location);
    vector<Node<2,VectorXd> *> nodes(3, nullptr);
	location << 0.0, 0.0;
    nodes[0] = new Node<2, VectorXd>(location);
	location << 1.0, 0.0;
    nodes[1] = new Node<2, VectorXd>(location);
	location << 1.0, 1.0;
    nodes[2] = new Node<2, VectorXd>(location);
    vector<SimplexFunction <VectorXd> > funcs(3);
    VectorXd coeffs(3);
    coeffs << -1,0,1;
    funcs[0].coeff = coeffs;
    coeffs << 1,-1,0;
    funcs[1].coeff = coeffs;
    coeffs << 0,1,0;
    funcs[2].coeff = coeffs;
	
    Element <2, 3, VectorXd> element(nodes, funcs);

	BoundaryConditions<VectorXd> boundaries;
	boundaries.cond = bound_cond_large_box;//omega = [0,2]^2
	boundaries.val = bound_val;
	


    SECTION( "Test default constructor()" ){
        Element <2, 3, VectorXd> empty_element;
		vector<Node <2, VectorXd> *> empty_nodes = empty_element.get_nodes();
		REQUIRE(empty_nodes[2] == nullptr);
		REQUIRE(empty_nodes.size() == 3);
    }
	
	SECTION("Test get_nodes") {

		vector<Node <2, VectorXd> *> retrieved_nodes = element.get_nodes();
		REQUIRE( retrieved_nodes[2]->how_many() == 6 );
		REQUIRE( retrieved_nodes[2]->get_location() == node3.get_location() );
		REQUIRE( retrieved_nodes[2]->get_index() == node3.get_index() );
		REQUIRE( retrieved_nodes[2]->get_shared_elements() == node3.get_shared_elements()+1 );
	}

    SECTION( "Test operator []" ){
        REQUIRE( element[0].get_location() == nodes[0]->get_location() );
        REQUIRE( element[0].get_location() != nodes[1]->get_location() );
    }

    /*SECTION( "Test show()" ){
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

	SECTION( "Test increase/decrease_shared_elements" ){
        int shared_els = element[0].get_shared_elements();
        element.increase_shared_elements();
        REQUIRE( shared_els+1 == element[0].get_shared_elements() );
		element.decrease_shared_elements();
		REQUIRE(shared_els == element[1].get_shared_elements());
    }

    SECTION( "Test set_indices" ){
        element.set_indices(34);
        REQUIRE( 35 == element[0].get_index() );
        REQUIRE( 36 == element[1].get_index() );
    }

	SECTION("Test set_all_indices_to") {
		element.set_all_indices_to(-3);
		REQUIRE(element[2].get_index() == -3);
		REQUIRE(element[0].get_index() == -3);
		REQUIRE(element[1].get_index() == -3);
	}

	SECTION("Test set_inner_indices") {
		element.set_all_indices_to(-1);
		REQUIRE(element.set_inner_node_indices(-1, boundaries) == 0);
		REQUIRE(element[2].get_index() == 0);
		REQUIRE(element[0].get_index() == -1);
		REQUIRE(element[1].get_index() == -1);
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
		cout << assigned_element[0].how_many() << endl;
    }
	
    SECTION( "Test operators == and !=" ){
        Element <2, 3, VectorXd> new_element(nodes, funcs);
        REQUIRE( new_element == element );
		vector<Node <2, VectorXd> *> reflected_nodes;//(3, nullptr);
        reflected_nodes.push_back(nodes[2]);
        reflected_nodes.push_back(nodes[1]);
        reflected_nodes.push_back(nodes[0]);
		Element <2, 3, VectorXd> refl_element(reflected_nodes, funcs);
		REQUIRE(reflected_nodes[1]->how_many() == 6 );
        REQUIRE( refl_element != element );
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
		int node_no = element[0].how_many();
		vector <pair <int[2], VectorXd> > mpoints = element.get_midpoints();
		vector<Node<2, VectorXd> *> new_nodes = element.get_midpoint_nodes();
		new_nodes[2]->show();
		REQUIRE( element[0].how_many() == node_no+3 );
		loc << 0.5, 0.0;
		REQUIRE(new_nodes[0]->get_location() == loc);
		loc << 0.5, 0.5;
		REQUIRE(new_nodes[1]->get_location() == loc);
		loc << 1.0, 0.5;
		REQUIRE(new_nodes[2]->get_location() == loc);
		//REQUIRE(midnodes[0].first[1] == 1);
	}

	SECTION("Test map of nodes") {
		vector <pair <int[2], VectorXd> > m_points = element.get_midpoints();
		vector< Node<2, VectorXd> * > m_nodes = element.get_midpoint_nodes();
		map< array<int, 2>, Node<2, VectorXd>* > node_map;
		for (int i = 0; i < 3; i++) {
			for (int j = i + 1; j < 3; j++) {
				node_map[{i, j}] = new Node<2, VectorXd>(*m_nodes[i+j-1]);
			}
		}
		node_map[{0, 1}]->show();
	}
	
}*/


TEST_CASE( "Test Element template containing Node template initiated with 2-D Point <2 ,double> template objects" ) {

    vector <double> vec1 = {0.0, 0.0};
    vector <double> vec2 = {1.0, 0.0};
    vector <double> vec3 = {1.0, 1.0};
    Point <2 ,double> point1(vec1);
    Point <2 ,double> point2(vec2);
    Point <2 ,double> point3(vec3);
	vector < Node <2, Point <2 ,double> >* > node_vec;
	node_vec.push_back(new Node <2, Point <2 ,double> >(point1));//delete is utilized by ~Element afters the tests are run
	node_vec.push_back(new Node <2, Point <2 ,double> >(point2));
	node_vec.push_back(new Node <2, Point <2 ,double> >(point3));
    vector <SimplexFunction <Point <2 ,double> > > fns(3);
    VectorXd coeff(3);
    coeff << -1,0,1;
    fns[0].coeff = coeff;
    coeff << 1,-1,0;
    fns[1].coeff = coeff;
    coeff << 0,1,0;
    fns[2].coeff = coeff;
    Element <2, 3, Point <2 ,double> > el(node_vec, fns);
	VolumeCalculator <2, Point <2, double> > volume_calculator;

	map< array<int, 2>, int> MIDPOINTS_MAP;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			MIDPOINTS_MAP.insert(pair< array<int, 2>, int>({ i, j }, i*3 + j - 1));
		}
	}

	SECTION("Test operator []") {
		cout << "Showing location in the Element based on Point<.,.>" << endl;
		node_vec[0]->show();
		REQUIRE(el[0].get_location() == node_vec[0]->get_location());
		REQUIRE(el[0].get_location() != node_vec[1]->get_location());
		REQUIRE(el[0].get_location()[0] == 0.0);
		REQUIRE(el[0].get_location()[1] == 0.0);
	}

    /*SECTION( "Test show()" ){
        cout << "showing element[0]" << endl;
        el[0].show();
        cout << "showing node1" << endl;
        node_vec[0]->show();
    }*/

	SECTION("Test get_squared_dist_mat()") {
		MatrixXd D = volume_calculator.get_distance_squared_matrix(el.get_nodes());
		cout << D << endl;
	}

    SECTION( "Test get_simplex_matrix(T &el)" ){
        Matrix<double, 2, 2> s_mat = volume_calculator.get_simplex_matrix(el.get_nodes());
        REQUIRE( s_mat == MatrixXd::Identity(2,2) );
    }

    SECTION( "Test get_volume()" ){
        REQUIRE( el.get_volume() == 0.5 );
		double vol = volume_calculator.get_dist_volume(el.get_nodes());
		REQUIRE(vol == 0.5);
    }

    SECTION( "Test copy constructor" ){
        Element <2, 3, Point <2 ,double> > copyed_el(el);
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
        el.set_indices(665);
        REQUIRE( 666 == el[0].get_index() );
        REQUIRE( 667 == el[1].get_index() );
    }

    SECTION( "Test assignment operator" ){
        Element <2, 3, Point <2 ,double> > assigned_el = el;
        REQUIRE( assigned_el[0] == el[0] );
        REQUIRE( assigned_el[1] == el[1] );
        REQUIRE( assigned_el[2] == el[2] );
        REQUIRE( assigned_el.get_function(0) == el.get_function(0) );
    }

    SECTION( "Test get_function(int)" ){
        SimplexFunction< Point <2 ,double> > func_1 = el.get_function(0);
        REQUIRE( func_1(node_vec[0]->get_location()) == 1 );
        REQUIRE( func_1(node_vec[1]->get_location()) == 0 );
        REQUIRE( func_1(node_vec[2]->get_location()) == 0 );
    }
	
	SECTION("Test midpoint_nodes()") {
		VectorXd loc(2);
		loc << 0.5, 0.0;
		int node_no = el[0].how_many();
		vector <pair <int[2], Point <2 ,double> > > m_points = el.get_midpoints();
		vector <Node <2, Point <2 ,double> >* > m_nodes = el.get_midpoint_nodes();
		REQUIRE(el[0].how_many() == node_no + 3);
		REQUIRE(m_nodes[0]->get_location()[0] == 0.5*(vec1[0]+ vec2[0]));
		REQUIRE(m_nodes[1]->get_location()[1] == 0.5*(vec1[1] + vec3[1]));
		REQUIRE(m_nodes[2]->get_location()[0] == 0.5*(vec2[0] + vec3[0]));
	}

}
