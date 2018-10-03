#include "../include/element.h"
#include "../include/mesh.h"
#include "../include/baseMesh.h"
#include "../include/FunctionHandler.h"

#include <math.h>

using namespace std;

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


TEST_CASE( "Test FunctionHandler" ) {


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
    cout << new_location << endl;
    cout << (*el_vec[1])[1].get_location() << endl;
    cout << "vec_size: " << el_vec.size() << endl;
    //element.show();
    //Mesh <2, 3, VectorXd> mesh;
    //mesh.show();
    //Mesh <2, 3, VectorXd> second_mesh(element);
    //second_mesh.show();
    /*SECTION( "Test constructing SimpleMesh(T &t) with element as arg" ){
        BaseMesh <Element <2, 3, VectorXd> > second_mesh;
        second_mesh.set_top(element);
        second_mesh.show();
    }
*/
    SECTION( "Test constructing FunctionHandler" ){
        BaseMesh <Element <2, 3, VectorXd> > b_mesh;
        //b_mesh.show();
    }

    SECTION( "Test constructing BaseMesh with int as Type" ){
        BaseMesh <int> s_mesh;
        s_mesh.set_top(1);
        //s_mesh.show();
    }

    //SECTION( "Test constructing Mesh" ){
        //Mesh <2, 3, VectorXd> mesh(element);
        //mesh.show();
//}

}

/*
TEST_CASE( "Test SimpleMesh" ) {

    /*VectorXd location(2);
    location << 0.0, 0.0;
    Node <2,VectorXd> node1(location);
    location << 1.0, 0.0;
    Node <2,VectorXd> node2(location);
    location << 1.0, 1.0;
    Node <2,VectorXd> node3(location);

    float loc1 = 1.0;
    float loc2 = 2.0;
    float loc3 = 3.0;
    Node <1, float> node1(loc1);
    Node <1, float> node2(loc2);
    Node <1, float> node3(loc3);
    Node <1, float> *nodes[3];
    nodes[0] = &node1;
    nodes[1] = &node2;
    nodes[2] = &node3;
    Element <1,3, float> el(nodes);
    //el.show();
    //vector < Element < 2,3, Vector2d > >* els;
    //els[0] = el;
    //els->push_back(&el);
    //els[0].show();
    //Mesh <2,3, VectorXd> mesh(el);
    Node <1, float> *new_nodes[1];
    new_nodes[0] = &node1;
    Element <1, 1, float> new_el(new_nodes);
    typedef Element <1,3, float> FloatElement;
    //cout << mesh.weak_form_element(1,1) << endl;
    FloatElement float_el(nodes);
    SimpleMesh < FloatElement > simple_mesh(new_el);
    simple_mesh.show();


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


}*/



