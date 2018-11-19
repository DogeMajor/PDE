#ifndef VOLUME_CALCULATOR_H
#define VOLUME_CALCULATOR_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include "Vertex.h"
#include <vector>
#include <math.h>
#include "Function.h"
#include "HelpfulTools.h"

using namespace std;
using namespace Eigen;

int factorial(int n) {
	return (n != 0) ? n * factorial(n - 1) : 1;
}

template <int Dimension, typename T>
class VolumeCalculator {
public:
	MatrixXd get_distance_squared_matrix(vector <Vertex<Dimension, T>* > el_vertices) const;
	double get_volume(vector <Vertex<Dimension, T>* > el_vertices) const;

};

template <int Dimension, typename T>
MatrixXd VolumeCalculator<Dimension, T>::get_distance_squared_matrix(vector <Vertex<Dimension, T>* > el_vertices) const {
	int cols = el_vertices.size();
	MatrixXd D = MatrixXd::Zero(cols+1, cols+1);
	for (int col = 0; col < cols; col++) {
		for (int row = 0; row < cols; row++) {
			D(row, col) = dist_squared<Dimension, T>(el_vertices[row]->get_location(), el_vertices[col]->get_location());
		}
		D(cols, col) = 1;
	}
	for (int row = 0; row < cols; row++) {
		D(row, cols) = 1;
	}
	return D;
}

template <int Dimension, typename T>
double VolumeCalculator<Dimension, T>::get_volume(vector <Vertex<Dimension, T>* > el_vertices) const {
	int simplex_order = el_vertices.size() - 1;
	MatrixXd D = get_distance_squared_matrix(el_vertices);
	double result = sqrt(abs(D.determinant()*pow(factorial(simplex_order), -2)*pow(2, -simplex_order)));
	return result;
}

#endif
