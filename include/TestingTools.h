#ifndef TESTINGTOOLS_H
#define TESTINGTOOLS_H

#include <map>
#include <vector>
#include "../include/Point.h"
#include "../include/Vertex.h"
#include "../include/Function.h"
#include "../include/Element.h"
#include "../include/ElementFactory.h"
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"

using namespace std;
using namespace Eigen;


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


bool point_bound_cond(Point<2, double> coords, double acc) {
	for (int i = 0; i < coords.size(); i++) {
		if ((abs(coords[i] - 0.0) < acc) || (abs(coords[i] - 1.0) < acc)) { return true; }
	}
	return false;
}

bool point_bound_is_inside(Point<2, double> coords, double acc) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] < 0.0 + acc) || (coords[i] > 1.0 - acc)) { return false; }
	}
	return true;
}

double point_bound_val(Point<2, double> coords) {
	//if (coords[1] == 1.0) { return 1; }
	return 0;
}

Point<2, double> point_bound_normal(Point<2, double> coords) {
	vector<double> normal = { 0.0,0.0 };
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) {
			normal[i] = (coords[i] == 0.0) ? -1 : 1;
			return Point<2, double>(normal);
		}
	}
	return Point<2, double>(normal);
}

//sqrt(2) circle at (.5,.5) boundaries for Point template

bool p_circle_cond(Point<2, double> coords, double acc) {
	vector<double> orig = { 0.5, 0.5 };
	double dist = sqrt(dist_squared<2, vector<double> >(orig, coords.get_value()));
	if ((dist > 1 / sqrt(2) - acc) && (dist < 1 / sqrt(2) + acc)) { return true; }
	return false;
}

bool p_circle_is_inside(Point<2, double> coords, double acc) {
	vector<double> orig = { 0.5, 0.5 };
	double dist = sqrt(dist_squared<2, vector<double> >(orig, coords.get_value()));
	if (dist < 1/sqrt(2) - acc) { return true; }
	return false;
}

double p_circle_val(Point<2, double> coords) {
	//if (coords[1] == 1.0) { return 1; }
	return 0;
}


//N-dim box's boundary with VectorXd as Typename

bool bound_cond(VectorXd coords, double acc) {
	for (int i = 0; i < coords.size(); i++) {
		if ((abs(coords[i] - 0.0) < acc) || (abs(coords[i] - 1.0) < acc)) { return true; }
	}
	return false;
}

bool bound_is_inside(VectorXd coords, double acc) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] < 0.0 + acc) || (coords[i] > 1.0 - acc)) { return false; }
	}
	return true;
}

double bound_val(VectorXd coords) {
	if (coords[1] == 1.0) { return 1; }
	return 0;
}

VectorXd bound_normal(VectorXd coords) {
	int sz = coords.size();
	VectorXd result = VectorXd::Zero(sz);
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] == 0.0) || (coords[i] == 1.0)) {
			result(i) = (coords[i] == 0.0) ? -1 : 1;
			return result;
		}
	}
	return result;
}

//Unit circle domain!!
bool c_cond(VectorXd coords, double acc) {
	VectorXd m_p(2);
	m_p << 0, 0;
	double dist = sqrt(dist_squared<2, VectorXd>(m_p, coords));
	if ((dist > 1 - acc) && (dist < 1 + acc)) { return true; }
	return false;
}

bool c_is_inside(VectorXd coords, double acc) {
	VectorXd m_p(2);
	m_p << 0, 0;
	if (sqrt(dist_squared<2, VectorXd>(m_p, coords)) < 1 - acc) { return true; }
	return false;
}

double c_val(VectorXd coords) {
	double angle = atan(coords[1] / coords[0]);
	return cos(10.0 * angle);
}

VectorXd c_normal(VectorXd coords) {
	double length = sqrt(coords.transpose()*coords);
	return (1 / length) * coords;
}

#endif
