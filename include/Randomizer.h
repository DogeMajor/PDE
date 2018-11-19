#ifndef RANDOMIZER_H
#define RANDOMIZER_H
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <random>
#include <chrono>


using namespace std;

typedef std::mt19937 Engine;
typedef std::uniform_real_distribution<> Distribution;

class Generator {

private:
	Engine eng;
	Distribution dist;

public:
	Generator() {}
	Generator(double min, double max, int seed) : dist(min, max) {
		eng.seed(seed);
	}
	Distribution::result_type gen() { return dist(eng); }
};



class Randomizer {
public:
	Randomizer(int seed) { generator = Generator(0, 1, seed); }
	~Randomizer() {}
	double prob() { return generator.gen(); }
	int random_int(int max) { return int(max*generator.gen()); }
	vector<double> randomize_items(vector<double> x);
	vector<double> get_convex_coeffs(int dim);

private:
	Generator generator;
};

vector<double> Randomizer::get_convex_coeffs(int dim) {
	vector <double> coeffs(dim);
	double sum = 0.0;
	for (int i = 0; i < dim - 1; i++) {
		coeffs[i] = (1 - sum)*prob();
		sum += coeffs[i];
	}
	coeffs[dim - 1] = 1 - sum;
	return randomize_items(coeffs);
}

vector<double> Randomizer::randomize_items(vector<double> x) {
	vector<double> items;
	int index;
	int sz = x.size();
	for (int i = 0; i < sz; i++) {
		index = random_int(sz - i);
		items.push_back(x[index]);
		x.erase(x.begin() + index);
	}
	return items;
}


class Seeder {
private:
	chrono::high_resolution_clock::time_point start_time;
public:
	Seeder() {
		start_time = chrono::high_resolution_clock::now();
	}

	int get_nanoseconds() const {
		chrono::high_resolution_clock::time_point current_time =
			chrono::high_resolution_clock::now();
		return chrono::duration_cast<chrono::nanoseconds>(current_time - start_time).count();
	}
	int get_milliseconds() const {
		chrono::high_resolution_clock::time_point current_time =
			chrono::high_resolution_clock::now();
		return chrono::duration_cast<chrono::milliseconds>(current_time - start_time).count();
	}
	int get_seconds() const {
		chrono::high_resolution_clock::time_point current_time =
			chrono::high_resolution_clock::now();
		return chrono::duration_cast<chrono::seconds>(current_time - start_time).count();
	}
};



#endif
