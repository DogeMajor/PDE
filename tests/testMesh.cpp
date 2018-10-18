#include "../include/point.h"
#include "../include/element.h"
#include "../include/mesh.h"
#include "../include/baseMesh.h"

#include <math.h>

using namespace std;

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"

/*
TEST_CASE( "Test Element template containing Node template initiated with 2-D double vector from Eigen lib" ) {


    VectorXd location(2);
    location << 0.0, 0.0;
    Node <2,VectorXd> node1(location);
    location << 1.0, 0.0;
    Node <2,VectorXd> node2(location);
    location << 1.0, 1.0;
    Node <2,VectorXd> node3(location);
    Node <2,VectorXd> *nodes[3];
    nodes[0] = &node1;
    nodes[1] = &node2;
    nodes[2] = &node3;

    Element <2, 3, VectorXd> element(nodes);
    element.increase_shared_elements();
    vector <Element <2, 3, VectorXd>* > el_vec(2);
    el_vec[0] = new Element <2, 3, VectorXd>(nodes);
    el_vec[1] = new Element <2, 3, VectorXd>(nodes);
    //el_vec[0]->show();
    VectorXd new_location = (*el_vec[0])[0].get_location();
    //cout << new_location << endl;
    //cout << (*el_vec[1])[1].get_location() << endl;
    //cout << "vec_size: " << el_vec.size() << endl;
    BaseMesh <Element <2, 3, VectorXd> > mesh;
    //element.show();
    //Mesh <2, 3, VectorXd> mesh;
    //mesh.show();
    //Mesh <2, 3, VectorXd> second_mesh(element);
    //second_mesh.show();

    SECTION( "Test constructing BaseMesh()" ){
        BaseMesh <Element <2, 3, VectorXd> > b_mesh;
        //b_mesh.show();
    }

    SECTION( "Test constructing BaseMesh with int as Type" ){
        BaseMesh <int> s_mesh;
        s_mesh.set_top(1);
        s_mesh.show();
    }

    SECTION( "Test assignment operator" ){
        Mesh <2, 3, VectorXd> second_mesh(element);
        Mesh <2, 3, VectorXd> third_mesh(element);
        REQUIRE( third_mesh.get_next() == second_mesh.get_next() );
        REQUIRE( third_mesh.get_element() == second_mesh.get_element() );
    }


}*/


TEST_CASE( "Test SimpleMesh with Elements based on Points" ) {

    vector <double> vec1 = {0.0, 0.0};
    vector <double> vec2 = {1.0, 0.0};
    vector <double> vec3 = {1.0, 1.0};
    Point <double> point1(vec1);
    Point <double> point2(vec2);
    Point <double> point3(vec3);
    Node <2,Point <double> > n_1(point1);
    Node <2,Point <double> > n_2(point2);
    Node <2,Point <double> > n_3(point3);
    vector<Node <2,Point <double> > * > node_vec(3, nullptr);
    node_vec[0] = &n_1;
    node_vec[1] = &n_2;
    node_vec[2] = &n_3;
	ElementFactory <2, 3, Point <double> > factory;
    Element <2, 3, Point <double> > el = factory.build(node_vec);
	BaseMesh < Element <2, 3, Point <double> > > empty_mesh;
    //BaseMesh < Element <2, 3, Point <double> > > mesh1(el);
    //mesh1.show();


    /*SECTION( "Test assignment operator" ){
        BaseMesh < Element <2, 3, Point <double> > > mesh2 = mesh1;
        REQUIRE( mesh2.get_next() == mesh1.get_next() );
        REQUIRE( mesh2.get_element() == mesh1.get_element() );
    }

    SECTION( "Test operators == and !=" ){
        BaseMesh < Element <2, 3, Point <double> > > similar_mesh(el);
        BaseMesh < Element <2, 3, Point <double> > > mesh3;
        REQUIRE( mesh1 == similar_mesh );
        REQUIRE( mesh1 == mesh1 );
        REQUIRE( mesh3 != mesh1 );
    }

    SECTION( "Test get_element" ){
        REQUIRE( mesh1.get_element() == el );
    }

    SECTION( "Test get_next" ){
        REQUIRE( mesh1.get_next() == nullptr );
    }

    SECTION( "Test set_top" ){
        BaseMesh < Element <2, 3, Point <double> > > empty_mesh;
        empty_mesh.set_top(el);
        REQUIRE( mesh1.get_element() == el );
    }

/*
    SECTION( "Test get_simplex_matrix" ){
        Matrix<double,2,2> mat = mesh.get_simplex_matrix(el);
        cout << mat << endl;
        REQUIRE( mat == MatrixXd::Identity(2,2) );
    }

    SECTION( "Test get_element_volume" ){
        double volume = mesh.get_element_volume(el);
        REQUIRE( volume == 1.0 );
    }
*/

}


