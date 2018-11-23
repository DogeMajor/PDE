#ifndef DAO_H
#define DAO_H

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "Mesh.h"

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
class BaseDAO {

private:
	Delimiters delimiters;

public:
	BaseDAO() {}
	~BaseDAO() {}
	void save_matrix(string file_name, MatrixXd mat, OutputType type = Mathematica);
	MatrixXd get_grid_values(Mesh<Dim, N, T> * m_ptr);
};

template <int Dim, int N, typename T>
MatrixXd BaseDAO<Dim, N, T>::get_grid_values(Mesh<Dim, N, T> * m_ptr) {
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

template <int Dim, int N, typename T>
void BaseDAO<Dim, N, T>::save_matrix(string file_name, MatrixXd mat, OutputType type) {

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
		cout << "Saving matrix with name " << file_name << " failed!" << endl;
	}
}

//------------------------------------------------------------------------------------------------------------------

template <int Dim, int N, typename T>
class MeshDAO: public BaseDAO<Dim, N, T> {

public:
	MeshDAO() {}
	~MeshDAO() {}
	void save_grid(string file_name, Mesh<Dim, N, T> * m_ptr, OutputType type = Mathematica);
};


template <int Dim, int N, typename T>
void MeshDAO<Dim, N, T>::save_grid(string file_name, Mesh<Dim, N, T> * m_ptr, OutputType type = Mathematica) {
	MatrixXd grid_m = get_grid_values(m_ptr);
	save_matrix(file_name, grid_m, type);
}

//-------------------------------------------------------------------------

template <int Dim, int N, typename T>
class SolverDAO : public BaseDAO<Dim, N, T> {

public:
	SolverDAO() {}
	~SolverDAO() {}
	MatrixXd get_solution_values(Mesh<Dim, N, T> * m_ptr, VectorXd solution);
	void save_solution_values(string file_name, Mesh<Dim, N, T> * m_ptr, VectorXd solution, OutputType type = Mathematica);
};

template <int Dim, int N, typename T>
MatrixXd SolverDAO<Dim, N, T>::get_solution_values(Mesh<Dim, N, T> * m_ptr, VectorXd solution) {
	MatrixXd grid = get_grid_values(m_ptr);
	MatrixXd values = MatrixXd::Zero(solution.size(), Dim + 1);
	int rows = solution.size();
	int cols = Dim + 1;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols-1; j++) {
			values(i, j) = grid(i, j);
		}
		values(i, cols - 1) = solution[i];
	}
	return values;
}

template <int Dim, int N, typename T>
void SolverDAO<Dim, N, T>::save_solution_values(string file_name, Mesh<Dim, N, T> * m_ptr, VectorXd solution, OutputType type = Mathematica) {
	MatrixXd solution_values = get_solution_values(m_ptr, solution);
	save_matrix(file_name, solution_values, type);
}

#endif
