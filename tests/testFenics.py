#!/usr/bin/env python
# -*- coding: utf-8 -*-

from fenics import *
from datetime import datetime
import codecs
import math
from xml.etree import ElementTree as ET

def file_to_etree(file_name):
    tree_root = ET.parse(file_name)
    return tree_root


def find_file_name(e_tree):
    data_sets = e_tree.findall(".//DataSet")
    file_name = data_sets[0].attrib['file']
    return file_name


def get_raw_grid_points(e_tree):
    point_els = e_tree.findall(".//Points/DataArray")
    points = point_els[0]
    return points.text

def get_solution(e_tree):
    sols = e_tree.findall(".//PointData/DataArray")
    sol = sols[0]
    return sol.text


def to_list(data_string, delimiter, to_last=0):
    if(to_last == 0):
        return data_string.split(delimiter)
    return data_string.split(delimiter)[:to_last]

def to_matrix(data_string, row_delimiter, item_delimiter):
    row_strings = to_list(data_string, row_delimiter, -1)
    matrix = [to_list(row, item_delimiter) for row in row_strings]
    return matrix

def to_solution_matrix(matrix, sol_vec):
    max_col = len(matrix[0])-1
    for index, row in enumerate(matrix):
        matrix[index][-1] = to_decimal_number(sol_vec[index])
    return matrix
    
def to_decimal_number(number):#seven decimals!
    return "{:.12f}".format(float(number))

def to_string(data_list, delimiter=","):
    return '{'+delimiter.join(data_list)+'}'

def to_mathematica_string(matrix, row_start, row_end):
    row_strings = [to_string(row, ',') for row in matrix]
    return '{'+','.join(row_strings)+'}'

def analytic_sol(x,y):
    x_ = float(to_decimal_number(x))*math.pi
    y_ = float(to_decimal_number(y))*math.pi
    return ((math.pi*math.pi*2)**-1)*math.sin(x_)*math.sin(y_)

def error_norm(sol_mat):
    errors = [(float(to_decimal_number(coords[2])) - analytic_sol(coords[0], coords[1]))**2 for coords in sol_mat]
    return (sum(errors))**.5/len(errors)

def get_solution_summary():#Just a script to compare pde program's results to fenics
    summary = {}
    parsed_tree = file_to_etree('poisson/solution.pvd')
    file_name = find_file_name(parsed_tree)
    sol_tree = file_to_etree('poisson/'+file_name)
    grid_points = get_raw_grid_points(sol_tree)
    sol_points = get_solution(sol_tree)
    sol_list = to_list(sol_points, "  ", -1)
    grid_mat = to_matrix(grid_points, "  ", " ")
    summary['sol_mat'] = to_solution_matrix(grid_mat, sol_list)
    summary['mathematica string'] = to_mathematica_string(summary['sol_mat'], '{', '},')
    summary['mesh size'] = len(summary['sol_mat'])
    summary['error'] = error_norm(summary['sol_mat'])
    return summary

# Create mesh and define function space
start = datetime.now()
mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
#u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
u_D = Expression('0', degree=0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
#f = Constant(-6.0)
f = Expression('sin(pi* x[0])*sin(pi* x[1])', degree=1)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

t_mesh_building = datetime.now() - start
#print(t_mesh_building.second)
start = datetime.now()

if __name__ == "__main__":
    # Compute solution
    
    u = Function(V)
    start = datetime.now()
    solve(a == L, u, bc)
    t_mesh_solving = datetime.now() - start
    # Save solution to file in VTK format
    vtkfile = File('poisson/solution.pvd')
    vtkfile << u

    # Compute error in L2 norm
    error_L2 = errornorm(u_D, u, 'L2')

    # Compute maximum error at vertices
    vertex_values_u_D = u_D.compute_vertex_values(mesh)
    vertex_values_u = u.compute_vertex_values(mesh)

    print("Time elapsed building the mesh and time elapsed solving the pde: {0}, {1}".format(t_mesh_building, t_mesh_solving))
    
    summary = get_solution_summary()
    print(summary['error'])
    print(summary['mesh size'])
    print(summary['mathematica string'])