#include "../include/point.h"
#include "../include/node.h"
#include "../C++ libs/eigen/Eigen/Dense"

using namespace std;
using namespace Eigen;


#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../C++ libs/catch/catch.hpp"

TEST_CASE("Test NodeFactory") {
	VectorXd loc(3);
	loc << 1.0, 2.0, 3.0;
	Node <3, VectorXd> node(loc);
	VectorXd sim_loc(3);
	sim_loc << 1.0, 2.0, 3.0;
	NodeFactory <3, VectorXd> factory;
	cout << "Show size of a VectorXd(3)" << endl;
	cout << sim_loc.size() << endl;

	SECTION("Test default constructor") {
		NodeFactory <3, VectorXd> extra_factory;
	}

	SECTION("Building Node based on loation T loc should succeed") {
		Node <3, VectorXd> new_node = factory.build(sim_loc);
		new_node.show();
		REQUIRE(node == new_node);
		REQUIRE(node.get_location() == new_node.get_location());
		REQUIRE(node.get_location() == loc);
	}

	SECTION("Building nodes based on T=Point<Dim, double> should succeed") {
		vector <double> coords = { 1.0, 2.0, 3.0 };
		Point <3, double> point1(coords);

		NodeFactory <3, Point<3, double> > factory2;
		Node <3, Point<3, double> > node1 = factory2.build(point1);
		cout << "Showing Node <3, Point<3,double> >" << endl;
		Point<3, double> point_loc = node1.get_location();
		for (int i = 0; i < 3; i++) {
			cout << point_loc[i] << endl;
		}
		node1.show();
		for (int i = 0; i < point_loc.size(); i++) {
			REQUIRE(node1.get_location()[i] == coords[i]);
		}
		//REQUIRE(node1.get_location() == coords);
		REQUIRE(node1.get_index() == -1);
		REQUIRE(node1.get_shared_elements() == 0);
	}
	

}