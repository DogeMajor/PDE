#ifndef VOLUME_CALCULATOR_H
#define VOLUME_CALCULATOR_H
#include <iostream>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"
#include "node.h"
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
	Matrix<double, Dimension, Dimension> get_simplex_matrix(vector <Node<Dimension, T>* > el_nodes) const;
	double get_volume(vector <Node<Dimension, T>* > el_nodes) const;

};

template <int Dimension, typename T>
Matrix<double, Dimension, Dimension> VolumeCalculator<Dimension, T>::get_simplex_matrix(vector <Node<Dimension, T>* > el_nodes) const {
	MatrixXd simplex_mat = MatrixXd::Zero(Dimension, Dimension);
	for (int col = 0; col < Dimension; col++) {
		for (int row = 0; row < Dimension; row++) {
			simplex_mat(row, col) = el_nodes[row + 1]->get_location()[col] - el_nodes[row]->get_location()[col];
		}
	}
	return simplex_mat;
}

template <int Dimension, typename T>
double VolumeCalculator<Dimension, T>::get_volume(vector <Node<Dimension, T>* > el_nodes) const {
	MatrixXd simplex_mat = get_simplex_matrix(el_nodes);
	return abs(simplex_mat.determinant() / factorial(Dimension));
}

#endif
