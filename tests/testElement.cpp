#include "../include/Point.h"
#include "../include/Vertex.h"
#include "../include/Element.h"

using namespace std;
using namespace Eigen;

//N-dim box's boundary
bool bound_cond(VectorXd coords, double acc) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) { return true; }
	}
	return false;
}

bool bound_cond_large_box(VectorXd coords, double acc) {
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


TEST_CASE( "Test Element template containing vertex template initiated with 2-D double vector from Eigen lib" ){

    VectorXd location(2);
    location << 0.0, 0.0;
    Vertex<2,VectorXd> vertex1(location);
    location << 1.0, 0.0;
    Vertex<2,VectorXd> vertex2(location);
    location << 1.0, 1.0;
    Vertex<2,VectorXd> vertex3(location);
    vector<Vertex<2,VectorXd> *> vertices(3, nullptr);
	location << 0.0, 0.0;
    vertices[0] = new Vertex<2, VectorXd>(location);
	location << 1.0, 0.0;
    vertices[1] = new Vertex<2, VectorXd>(location);
	location << 1.0, 1.0;
    vertices[2] = new Vertex<2, VectorXd>(location);
    vector<SimplexFunction <VectorXd> > funcs(3);
    VectorXd coeffs(3);
    coeffs << -1,0,1;
    funcs[0].coeff = coeffs;
    coeffs << 1,-1,0;
    funcs[1].coeff = coeffs;
    coeffs << 0,1,0;
    funcs[2].coeff = coeffs;
	
    Element <2, 3, VectorXd> element(vertices, funcs);

	BoundaryConditions<VectorXd> boundaries;
	boundaries.cond_fn = bound_cond_large_box;//omega = [0,2]^2
	boundaries.val = bound_val;
	


    SECTION( "Test default constructor()" ){
        Element <2, 3, VectorXd> empty_element;
		vector<Vertex<2, VectorXd> *> empty_vertices = empty_element.get_vertices();
		REQUIRE(empty_vertices[2] == nullptr);
		REQUIRE(empty_vertices.size() == 3);
    }
	
	SECTION("Test get_vertices") {

		vector<Vertex<2, VectorXd> *> retrieved_vertices = element.get_vertices();
		REQUIRE( retrieved_vertices[2]->how_many() == 6 );
		REQUIRE( retrieved_vertices[2]->get_location() == vertex3.get_location() );
		REQUIRE( retrieved_vertices[2]->get_index() == vertex3.get_index() );
		REQUIRE( retrieved_vertices[2]->get_shared_elements() == vertex3.get_shared_elements()+1 );
	}

    SECTION( "Test operator []" ){
        REQUIRE( element[0].get_location() == vertices[0]->get_location() );
        REQUIRE( element[0].get_location() != vertices[1]->get_location() );
    }

    SECTION( "Test show()" ){
        cout << "showing element[0]" << endl;
        element[0].show();
        cout << "showing vertex1" << endl;
        vertices[0]->show();
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
		REQUIRE(element.set_inner_vertex_indices(-1, boundaries) == 0);
		REQUIRE(element[2].get_index() == 0);
		REQUIRE(element[0].get_index() == -1);
		REQUIRE(element[1].get_index() == -1);
	}

	SECTION("Test set index maps") {
		element.set_indices(0);
		element.set_index_maps();
		cout << "0 to global" << element.to_global(0) << endl;
		REQUIRE(element.to_global(0) == 1);
		REQUIRE(element.to_global(1) == 2);
		REQUIRE(element.to_global(2) == 3);
		REQUIRE(element.to_local(1) == 0);
		REQUIRE(element.to_local(2) == 1);
		REQUIRE(element.to_local(3) == 2);
		element.set_all_indices_to(-1);
	}

	SECTION("Test f_variations methods") {
		element.set_f_variation(0, .1);
		element.set_f_variation(2, .5);
		vector<double> f_vars = element.get_f_variations();
		REQUIRE(f_vars[0] == 0.1);
		REQUIRE(f_vars[1] == 0);
		REQUIRE(f_vars[2] == 0.5);
		REQUIRE(limit_decimals(element.get_avg_f_variation(),2) == 0.2);
	}

    SECTION( "Test get_function(int)" ){
        SimplexFunction<VectorXd> first_fn = element.get_function(0);
        REQUIRE( first_fn(vertex1.get_location()) == 1 );
        REQUIRE( first_fn(vertex2.get_location()) == 0 );
        REQUIRE( first_fn(vertex3.get_location()) == 0 );
    }

    SECTION( "Test assignment operator" ){
        Element <2, 3, VectorXd> assigned_element = element;
        REQUIRE( assigned_element[0] == element[0] );
        REQUIRE( assigned_element[1] == element[1] );
        REQUIRE( assigned_element[2] == element[2] );
		cout << assigned_element[0].how_many() << endl;
    }
	
    SECTION( "Test operators == and !=" ){
        Element <2, 3, VectorXd> new_element(vertices, funcs);
        REQUIRE( new_element == element );
		vector<Vertex<2, VectorXd> *> reflected_vertices;//(3, nullptr);
        reflected_vertices.push_back(vertices[2]);
        reflected_vertices.push_back(vertices[1]);
        reflected_vertices.push_back(vertices[0]);
		Element <2, 3, VectorXd> refl_element(reflected_vertices, funcs);
		REQUIRE(reflected_vertices[1]->how_many() == 6 );
        REQUIRE( refl_element != element );
    }
	
    SECTION( "Test get_volume()" ){
        REQUIRE(element.get_volume() == 0.5 );
    }

}


TEST_CASE( "Test Element template containing vertex template initiated with 2-D Point <2 ,double> template objects" ) {

    vector <double> vec1 = {0.0, 0.0};
    vector <double> vec2 = {1.0, 0.0};
    vector <double> vec3 = {1.0, 1.0};
    Point <2 ,double> point1(vec1);
    Point <2 ,double> point2(vec2);
    Point <2 ,double> point3(vec3);
	vector < Vertex<2, Point <2 ,double> >* > vertex_vec;
	vertex_vec.push_back(new Vertex<2, Point <2 ,double> >(point1));//delete is utilized by ~Element afters the tests are run
	vertex_vec.push_back(new Vertex<2, Point <2 ,double> >(point2));
	vertex_vec.push_back(new Vertex<2, Point <2 ,double> >(point3));
    vector <SimplexFunction <Point <2 ,double> > > fns(3);
    VectorXd coeff(3);
    coeff << -1,0,1;
    fns[0].coeff = coeff;
    coeff << 1,-1,0;
    fns[1].coeff = coeff;
    coeff << 0,1,0;
    fns[2].coeff = coeff;
    Element <2, 3, Point <2 ,double> > el(vertex_vec, fns);
	VolumeCalculator <2, Point <2, double> > volume_calculator;

	map< array<int, 2>, int> MIDPOINTS_MAP;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			MIDPOINTS_MAP.insert(pair< array<int, 2>, int>({ i, j }, i*3 + j - 1));
		}
	}

	SECTION("Test operator []") {
		cout << "Showing location in the Element based on Point<.,.>" << endl;
		vertex_vec[0]->show();
		REQUIRE(el[0].get_location() == vertex_vec[0]->get_location());
		REQUIRE(el[0].get_location() != vertex_vec[1]->get_location());
		REQUIRE(el[0].get_location()[0] == 0.0);
		REQUIRE(el[0].get_location()[1] == 0.0);
	}

    SECTION( "Test show()" ){
        cout << "showing element[0]" << endl;
        el[0].show();
        cout << "showing vertex1" << endl;
        vertex_vec[0]->show();
    }

    SECTION( "Test get_distance_squared_matrix(T &el)" ){
        MatrixXd s_mat = volume_calculator.get_distance_squared_matrix(el.get_vertices());
		MatrixXd s_mat_should_be(4, 4);
		s_mat_should_be << 0, 1, 2, 1, 1, 0, 1, 1, 2, 1, 0, 1, 1, 1, 1, 0;
		REQUIRE( s_mat == s_mat_should_be);
    }

    SECTION( "Test get_volume()" ){
        REQUIRE( el.get_volume() == 0.5 );
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
        REQUIRE( func_1(vertex_vec[0]->get_location()) == 1 );
        REQUIRE( func_1(vertex_vec[1]->get_location()) == 0 );
        REQUIRE( func_1(vertex_vec[2]->get_location()) == 0 );
    }


}
