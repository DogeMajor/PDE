#include "../include/point.h"
#include "../include/element.h"
#include "../include/mesh.h"
//#include "../include/LinkedMesh.h"

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
*/

/*TEST_CASE("Test LinkedMesh with Elements based on Points") {

	vector <double> vec1 = { 0.0, 0.0 };
	vector <double> vec2 = { 1.0, 0.0 };
	vector <double> vec3 = { 1.0, 1.0 };
	vector <double> vec4 = { 0.0, 1.0 };
	Point <double> point1(vec1);
	Point <double> point2(vec2);
	Point <double> point3(vec3);
	Point <double> point4(vec4);
	Node <2, Point <double> > n_1(point1);
	Node <2, Point <double> > n_2(point2);
	Node <2, Point <double> > n_3(point3);
	vector<Node <2, Point <double>  > * > node_vec(3, nullptr);
	node_vec[0] = new Node<2, Point <double> >(point1);
	node_vec[1] = new Node<2, Point <double> >(point2);
	node_vec[2] = new Node<2, Point <double> >(point3);
	ElementFactory<2, 3, Point <double> > factory;
	Element<2, 3, Point <double> > el1 = factory.build(node_vec);
	vector<Node <2, Point <double> > *> node_vec2;//(3, nullptr);
	node_vec2.push_back(node_vec[0]);
	node_vec2.push_back(node_vec[2]);
	node_vec2.push_back(new Node<2, Point <double> >(point4));
	Element<2, 3, Point <double> > el2 = factory.build(node_vec2);
	LinkedMesh<Element <2, 3, Point <double> > > el_mesh(el1);
	el_mesh.push(el2);
	el_mesh.show();
	//el_mesh.get_top().show();

	LinkedMesh< Node<2, Point <double> > > node_mesh(*node_vec[0]);//ok


	SECTION("Mesh can be initialized with default constructor") {
		LinkedMesh<Element <2, 3, Point <double> > > empty_mesh;
		REQUIRE(empty_mesh.how_many() == 2);
	}

	SECTION("Creating mesh of Point objects should succeed") {
		LinkedMesh< Point <double> > point_mesh(point1);
		REQUIRE(point_mesh.how_many() == 1);
	}

	SECTION("One can access top element"){
		REQUIRE(el_mesh.get_top() == el2);
		REQUIRE(node_mesh.get_top() == *node_vec[0]);
	}

	SECTION("One can access last element") {
		REQUIRE(el_mesh.get_last() == el1);
		REQUIRE(node_mesh.get_last() == *node_vec[0]);
	}

	SECTION("Pushing item to top should succeed") {
		REQUIRE(node_mesh.push(*node_vec[1]));
	}

	SECTION("Popping item from top should succeed") {
		REQUIRE(node_mesh.pop());
	}

	SECTION("Refining the mesh should succeed") {
		el_mesh.refine();
		REQUIRE(el_mesh.how_many_nodes() == 8);
		//REQUIRE(node_mesh.get_last() == *node_vec[0]);
	}

	//SECTION( "Test copy constructors" ){
		//LinkedMesh < Element <2, 3, Point <double> > > el_mesh2(el_mesh);
		//REQUIRE( el_mesh2.get_top() == el_mesh.get_top() );
		//REQUIRE( el_mesh2.get_last() == el_mesh.get_last());
	//}
	/*
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


}*/


TEST_CASE("Test the real Mesh with Elements based on Points") {

	vector <double> vec1 = { 0.0, 0.0 };
	vector <double> vec2 = { 1.0, 0.0 };
	vector <double> vec3 = { 1.0, 1.0 };
	vector <double> vec4 = { 0.0, 1.0 };
	vector <double> vec5 = { 2.0, 0.0 };
	vector <double> vec6 = { 2.0, 1.0 };
	Point <double> point1(vec1);
	Point <double> point2(vec2);
	Point <double> point3(vec3);
	Point <double> point4(vec4);
	Point <double> point5(vec5);
	Point <double> point6(vec6);
	Node <2, Point <double> > n_1(point1);
	Node <2, Point <double> > n_2(point2);
	Node <2, Point <double> > n_3(point3);
	vector<Node <2, Point <double>  > * > node_vec(3, nullptr);
	node_vec[0] = new Node<2, Point <double> >(point1);
	node_vec[1] = new Node<2, Point <double> >(point2);
	node_vec[2] = new Node<2, Point <double> >(point3);
	ElementFactory<2, 3, Point <double> > factory;
	Element<2, 3, Point <double> > el1 = factory.build(node_vec);
	vector<Node <2, Point <double> > *> node_vec2;//(3, nullptr);
	node_vec2.push_back(node_vec[0]);
	node_vec2.push_back(node_vec[2]);
	node_vec2.push_back(new Node<2, Point <double> >(point4));
	Element<2, 3, Point <double> > el2 = factory.build(node_vec2);
	Mesh<2, 3, Point <double> > el_mesh(el1);
	el_mesh.push(el2);
	//el_mesh.show();
	//el_mesh.get_top().show();

	//LinkedMesh< Node<2, Point <double> > > node_mesh(*node_vec[0]);//ok


	/*SECTION("Mesh can be initialized with default constructor") {//OK
		Mesh<2, 3, Point <double> > empty_mesh;
		REQUIRE(empty_mesh.how_many() == 2);
	}

	SECTION("One can access top element") {//OK
		REQUIRE(el_mesh.get_top() == el2);
	}

	SECTION("One can access last element") {//OK
		REQUIRE(el_mesh.get_last() == el1);
	}

	SECTION("One can pop and push at the top") {//OK
		REQUIRE(el_mesh.pop());
		REQUIRE(el_mesh.get_top() == el1);
		REQUIRE(el_mesh.push(el2));
		REQUIRE(el_mesh.get_top() == el2);
		REQUIRE(el_mesh.get_last() == el1);
	}

	SECTION("One can pop a MeshNode from anywhere in the Mesh") {//OK!!
		MeshNode<Element<2, 3, Point <double> > >* top_mesh_node = el_mesh.get_top_mesh_node();
		REQUIRE(el_mesh.pop(top_mesh_node));
		REQUIRE(el_mesh.get_last() == el2);
		top_mesh_node = el_mesh.get_top_mesh_node();
		REQUIRE(el_mesh.pop(top_mesh_node) == false);
		//Rebuild the mesh!!
		REQUIRE(el_mesh.pop());
		REQUIRE(el_mesh.push(el1));
		REQUIRE(el_mesh.push(el2));
		REQUIRE(el_mesh.how_many_nodes() == 2);
		REQUIRE(el_mesh.get_last() == el1);
	}

	SECTION("One can push a MeshNode from anywhere in the Mesh") {//OK!!
		MeshNode<Element<2, 3, Point <double> > >* top_mesh_node = el_mesh.get_top_mesh_node();
		REQUIRE(el_mesh.push(top_mesh_node, el2));//Should be pushed at the middle place
		Element<2, 3, Point <double> > new_el(el1);
		REQUIRE(el_mesh.push(nullptr, new_el));
		REQUIRE(el_mesh.how_many_nodes() == 4);
		el_mesh.show();
		REQUIRE(el_mesh.get_last() == el1);
		//Rebuild the mesh!!
		REQUIRE(el_mesh.pop() == true);
		REQUIRE(el_mesh.pop(top_mesh_node) == true);
		REQUIRE(el_mesh.how_many_nodes() == 2);
	}
	SECTION("One can push and pop a MeshNode from anywhere in the Mesh") {//OK!!
		vector<Node <2, Point <double>  > * > node_vec_D;
		node_vec_D.push_back(new Node<2, Point <double> >(point6));
		node_vec_D.push_back(node_vec[1]);
		node_vec_D.push_back(node_vec[2]);
		Element<2, 3, Point <double> > el4 = factory.build(node_vec_D);
		

		vector<Node <2, Point <double>  > * > node_vec_C;
		node_vec_C.push_back(new Node<2, Point <double> >(point5));
		node_vec_C.push_back(node_vec[1]);
		node_vec_C.push_back(node_vec_D[0]);
		Element<2, 3, Point <double> > el3 = factory.build(node_vec_C);

		Mesh<2, 3, Point <double> > large_mesh(el1);

		SECTION("Pushing and popping mesh nodes succeeds everywhere in the mesh") {
			large_mesh.push(el4);
			large_mesh.push(large_mesh.get_top_mesh_node(), el3);
			large_mesh.push(large_mesh.get_top_mesh_node()->next, el2);
			//large_mesh.show();
			REQUIRE(large_mesh.get_element(0) == el4);
			REQUIRE(large_mesh.get_element(1) == el3);
			REQUIRE(large_mesh.get_element(2) == el2);
			REQUIRE(large_mesh.get_element(3) == el1);

			cout << "Det var det!" << endl;
			//REQUIRE(large_mesh.pop());

			REQUIRE(large_mesh.pop(large_mesh.get_top_mesh_node()->next));
			REQUIRE(large_mesh.get_element(2) == el1);
			REQUIRE(large_mesh.get_element(1) == el3);
			REQUIRE(large_mesh.get_element(0) == el4);

			REQUIRE(large_mesh.pop(large_mesh.get_top_mesh_node()->next->next) == false);
			REQUIRE(large_mesh.pop(large_mesh.get_top_mesh_node()->next));
			REQUIRE(large_mesh.get_element(1) == el3);
			REQUIRE(large_mesh.get_element(0) == el4);

			REQUIRE(large_mesh.pop());
			REQUIRE(large_mesh.get_element(0) == el3);

			REQUIRE(large_mesh.pop());
			REQUIRE(large_mesh.how_many_nodes() == 0);

			REQUIRE(large_mesh.pop() == false);


		}

	}*/

	SECTION("Refining the mesh should succeed") {
		el_mesh.refine();
		REQUIRE(el_mesh.how_many_nodes() == 8);
		vector<double> new_loc_22_1 = { 0, 1.0 };
		vector<double> new_loc_22_2 = { 0.5, 1.0 };
		el_mesh.show();
		REQUIRE(el_mesh.get_element(2)[0].get_location() == new_loc_22_1);
	}

	//SECTION( "Test copy constructors" ){
		//LinkedMesh < Element <2, 3, Point <double> > > el_mesh2(el_mesh);
		//REQUIRE( el_mesh2.get_top() == el_mesh.get_top() );
		//REQUIRE( el_mesh2.get_last() == el_mesh.get_last());
	//}
	/*
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




