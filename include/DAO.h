#ifndef DAO_H
#define DAO_H

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "Mesh.h"
using namespace std;

enum OutputType { Normal, Mathematica };


class Delimiters {

private:
	map<OutputType, map<string, string> > d_map;

public:
	Delimiters();
	map<string, string> get_delimiters(OutputType output) {
		return d_map[output];
	}
};


Delimiters::Delimiters() {
	d_map[Mathematica]["item"] = ",";
	d_map[Mathematica]["start_row"] = "{";
	d_map[Mathematica]["end_row"] = "},";
	d_map[Mathematica]["start_file"] = "{";
	d_map[Mathematica]["end_file"] = "}}";

	d_map[Normal]["item"] = " ";
	d_map[Normal]["start_row"] = "";
	d_map[Normal]["end_row"] = "";
	d_map[Normal]["start_file"] = "";
	d_map[Normal]["end_file"] = "";
}


template <int Dim, int N, typename T>
class MeshDAO {

private:
	Delimiters delimiters;

public:
	MeshDAO() {}
	~MeshDAO() {}
	void save_matrix(string file_name, MatrixXd mat, OutputType type = Mathematica);
	MatrixXd get_grid_values(Mesh<Dim, N, T> * m_ptr);
	void save_grid(string file_name, Mesh<Dim, N, T> * m_ptr, OutputType type = Mathematica);

};


template <int Dim, int N, typename T>
void MeshDAO<Dim, N, T>::save_grid(string file_name, Mesh<Dim, N, T> * m_ptr, OutputType type = Mathematica) {
	MatrixXd grid_m = get_grid_values(m_ptr);
	save_matrix(file_name, grid_m, type);
}

template <int Dim, int N, typename T>
void MeshDAO<Dim, N, T>::save_matrix(string file_name, MatrixXd mat, OutputType type) {

	map<string, string> delimiter = delimiters.get_delimiters(type);
	ofstream file;

	try {
		file.open(file_name);
		file << delimiter["start_file"];
		for (int row = 0; row < mat.rows(); row++) {
			file << delimiter["start_row"];
			for (int col = 0; col < mat.cols() - 1; col++) {
				file << mat(row, col) << delimiter["item"];
			}
			if (row != mat.rows() - 1) {
				file << mat(row, mat.cols() - 1) << delimiter["end_row"];
			}
			else {
				file << mat(row, mat.cols() - 1) << delimiter["end_file"];
			}
		}
		file.close();
	}
	catch (std::exception err) {
		cout << "Saving matrix with name "<< file_name <<  " failed!" << endl;
	}
}


template <int Dim, int N, typename T>
MatrixXd MeshDAO<Dim, N, T>::get_grid_values(Mesh<Dim, N, T> * m_ptr) {
	int max_outer_index = m_ptr->get_max_outer_index();
	MatrixXd values = MatrixXd::Zero(max_outer_index + 1, Dim);
	MeshNode <Element<Dim, N, T> >* iter = m_ptr->get_top_mesh_node();
	int I;
	T loc;
	while (iter != nullptr) {
		for (int i = 0; i < N; i++) {
			I = iter->data[i].get_index();
			loc = iter->data[i].get_location();
			for (int J = 0; J < Dim; J++) {
				values(I, J) = loc[J];
			}
		}
		iter = iter->next;
	}
	return values;
}


#endif
