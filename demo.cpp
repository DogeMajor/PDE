#include "../PDE/include/Solver.h"


using namespace std;
using namespace Eigen;


double PI = 3.1415926535;

//----------------boundary conditions (same as in TestingTools.h)--------------------------
//N-dim box's boundary with VectorXd as Typename

bool bound_cond(VectorXd coords, double acc) {
	for (int i = 0; i < coords.size(); i++) {
		if ((abs(coords[i] - 0.0) <= acc) || (abs(coords[i] - 1.0) <= acc)) { return true; }
	}
	return false;
}

bool bound_is_inside(VectorXd coords, double acc) {
	for (int i = 0; i < coords.size(); i++) {
		if ((coords[i] <= 0.0 + acc) || (coords[i] >= 1.0 - acc)) { return false; }
	}
	return true;
}

double bound_val(VectorXd coords) {
	//if (coords[1] == 1.0) { return 1; }
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

//------------------kernel functions of the weak form PDE--------------------------

double f_kern_sin(VectorXd coords) {//Particle in a N-Dim box...
	double result = sin(coords[0] * PI);
	for (int i = 1; i < coords.rows(); i++) {
		result *= sin(coords[i] * PI);
	}
	return result;
}


double analytic_sol(VectorXd coords) {//Particle in a N-Dim box...
	int dim = coords.size();
	double result = (1 / (dim*PI*PI))*sin(coords[0] * PI);
	for (int i = 1; i < coords.rows(); i++) {
		result *= sin(coords[i] * PI);
	}
	return result;
}


double error_norm(MatrixXd sol_values) {//For eventual future work with surface integrals
	double norm = 0;
	int sz = sol_values.cols() - 1;
	VectorXd loc;
	for (int i = 0; i < sol_values.rows(); i++) {
		loc = sol_values.row(i).head(sz);
		norm = norm + pow(sol_values(i, sz) - analytic_sol(loc), 2);
	}
	return norm * (1 / double(sol_values.rows()));
}

int main() { 
	cout << "Solves a poisson pde with sin(pi*x)sin(pi*y)sin(pi*z) as f, A as I(3,3) and boundary fn as zero" << endl;
	cout << "Domain is a simple unit cube centered at (.5,.5,.5)" << endl;
	
	BoundaryConditions<VectorXd> p_boundaries = { bound_cond, bound_is_inside, bound_val, bound_normal, 0.000001 };

	MeshDAO<3, 4, VectorXd> p_mesh_dao;
	SolverDAO<3, 4, VectorXd> p_solver_dao;

	Mesh<3, 4, VectorXd> p_empty_mesh;
	Mesh <3, 4, VectorXd>* p_mesh_ptr = &p_empty_mesh;

	BilinearFunction p_bl_fn;
	p_bl_fn.mat = MatrixXd::Identity(3, 3);
	PDE<3, VectorXd> p_pde(p_bl_fn, f_kern_sin);
	Solver<3, VectorXd> p_solver(p_pde, p_mesh_ptr, p_boundaries);
	Timer timer = Timer();
	
	VectorXd cube_center(3);
	cube_center << 0.5, 0.5, 0.5;
	VectorXd cube_lengths(3);
	cube_lengths << 1, 1, 1;

	p_solver.fill_mesh_covering_box(cube_center, cube_lengths);
	timer.reset();
	p_solver.refine();
	p_solver.refine();
	cout << "Refining the mesh twice took " << timer.get_milliseconds() << "ms." << endl;
	cout << "Max inner ind: " << p_mesh_ptr->get_max_inner_index() << endl;
	cout << "Max outer ind: " << p_mesh_ptr->get_max_outer_index() << endl;
	cout << "Mesh size: " << p_mesh_ptr->how_many_nodes() << endl;
	timer.reset();
	VectorXd p_sol = p_solver.solve();
	cout << "Solving the PDE took " << timer.get_milliseconds() << "ms." << endl;

	MatrixXd p_grid = p_mesh_dao.get_grid_values(p_mesh_ptr);
	MatrixXd p_sol_values = p_solver_dao.get_solution_values(p_mesh_ptr, p_sol);
	cout << "max value" << p_sol.maxCoeff() << endl;
	cout << "error norm: " << error_norm(p_sol_values) << endl;
	
	return 0; }

