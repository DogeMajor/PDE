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
	MatrixXd get_distance_squared_matrix(vector <Node<Dimension, T>* > el_nodes) const;
	Matrix<double, Dimension, Dimension> get_simplex_matrix(vector <Node<Dimension, T>* > el_nodes) const;
	double get_volume(vector <Node<Dimension, T>* > el_nodes) const;
	double get_dist_volume(vector <Node<Dimension, T>* > el_nodes) const;

};

template <int Dimension, typename T>
MatrixXd VolumeCalculator<Dimension, T>::get_distance_squared_matrix(vector <Node<Dimension, T>* > el_nodes) const {
	int cols = el_nodes.size();
	MatrixXd D = MatrixXd::Zero(cols+1, cols+1);
	for (int col = 0; col < cols; col++) {
		for (int row = 0; row < cols; row++) {
			D(row, col) = dist_squared<Dimension, T>(el_nodes[row]->get_location(), el_nodes[col]->get_location());
		}
		D(cols, col) = 1;
	}
	for (int row = 0; row < cols; row++) {
		D(row, cols) = 1;
	}
	return D;
}

template <int Dimension, typename T>
double VolumeCalculator<Dimension, T>::get_dist_volume(vector <Node<Dimension, T>* > el_nodes) const {
	int simplex_order = el_nodes.size()-1;
	MatrixXd D = get_distance_squared_matrix(el_nodes);
	double result = sqrt(abs(D.determinant()*pow(factorial(simplex_order), -2)*pow(2, -simplex_order)));
	return result;
}

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
