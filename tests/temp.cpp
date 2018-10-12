#include "../include/point.h"
#include "../include/node.h"
#include "../C++ libs/eigen/Eigen/Dense"

using namespace std;
using namespace Eigen;


#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"


TEST_CASE("Test Node template with 3-D double vector from Eigen lib") {
	VectorXd location(3);
	location << 1.0, 2.0, 3.0;
	Node <3, VectorXd> node(location);
	location << 1.0, 2.0, 3.0;
	Node <3, VectorXd> similar_node(location);

	SECTION("Test get_location") {
		VectorXd value = node.get_location();
		REQUIRE(value(0) == location(0));
	}

	SECTION("Test pointers and pointers to arrays") {
		int g[] = { 9,8 };
		int(*j)[2] = &g;

		//Dereference 'j' and access array element zero
		int n = (*j)[0];
		cout << n;
	}

	SECTION("Test pointers to arrays of templates") {
		Node <3, VectorXd> nodes[2] = { node, similar_node };//Instantiating arrays needs the objects right away, cannot be assigned later!!
		nodes[0].show();
		//nodes[1] = similar_node;
		Node <3, VectorXd>(*node_ptrs)[2] = &nodes;
		Node <3, VectorXd> deref_node = (*node_ptrs)[0];
		deref_node.show();

	}
	SECTION("Test constructing an array of object pointers") {
		VectorXd location(2);
		location << 0.0, 0.0;
		cout << "VectorXd size: " << location.size() << endl;
		Node <2, VectorXd> node1(location);
		location << 1.0, 0.0;
		Node <2, VectorXd> node2(location);
		location << 1.0, 1.0;
		Node <2, VectorXd> node3(location);
		Node <2, VectorXd> *nodes[3];
		nodes[0] = &node1;
		nodes[1] = &node2;
		nodes[2] = &node3;
		nodes[0]->show();
	}
	



	/*SECTION("Test vector of template class instances") {
		vector <Node <3, VectorXd> > nodes(2);
		nodes[0] = node;
		nodes[1] = similar_node;
		nodes[0].show();
		//Node <3, VectorXd>(*node_ptrs)[2] = &nodes;
		//Node <3, VectorXd> deref_node = (*node_ptrs)[0];
		//deref_node.show();

	}*/
	/*
	SECTION("Test pointers to arrays of templates") {
		//Node <3, VectorXd> nods[2];
		//nods[0].show();
		//nodes[1] = similar_node;
		Node <3, VectorXd>(*n_ptrs)[2];
		Node <3, VectorXd> n_1(node);
		n_ptrs[0] = &n_1;
		//(*n_ptrs)[1] = Node <3, VectorXd>(similar_node);
		//Node <3, VectorXd> der_node = (*n_ptrs)[0];
		//cout << "Hehei!" << endl;
		//cout << der_node.get_location();

	}*/
}