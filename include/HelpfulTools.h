#ifndef HELPFULTOOLS_H
#define HELPFULTOOLS_H

#include <math.h>
#include <array>
#include <map>
#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Core"

using namespace Eigen;

double limit_decimals(double number, int decimals){
    double N = pow(10, decimals);
    return double(int(number * N)) / N;
}

typedef map< array<int, 2>, int> Map;
typedef map< array<int, 2>, int>::const_iterator MapIter;

void show_map(map< array<int, 2>, int> &m){
	for (MapIter iter = m.begin(); iter != m.end(); iter++){
		cout << "Key: " << iter->first[0] << ", " << iter->first[1] << endl << "Value:" << endl;
		cout << iter->second << endl;
	}
}


template <int Dim, typename T>
double dist_squared(const T &a, const T &b){
	double result = 0;
	for (int i = 0; i < Dim; i++) {
		result += pow(a[i] - b[i], 2);
	}
	return result;
}

template <typename T>
pair<int, double> find_smallest(const T &container, int size) {
	int index = 0;
	for (int i = 1; i < size; i++) {
		if (container[i] < container[index]) {s
			index = i;
		}
	}
	return pair<int, double>(index, container[index]);
}

#endif