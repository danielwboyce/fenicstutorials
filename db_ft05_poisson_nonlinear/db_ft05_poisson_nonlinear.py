#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 17:12:02 2019

@author: daniel
"""

"""
Gonna solve for the nonlinear Poisson equation

    -div(q(u) * grad(u)) = f    in the unit square
    
    u  = u_D                    on the boundary (boundary condition).
                                Chosen to be zero. Dirichlet boundary.
    u  = u_0                    at t = 0 (initial condition).
                                Chosen to be a Gaussian hill. exp(-a*x**2 - a*y**2)

"""

import fenics as fs
import sympy as sym
import numpy as np

def q(u):
    "Return nonlinear coefficient"
    return 1 + u**2

# use SymPy to compute f form the manufactured solution u
x, y = sym.symbols('x[0], x[1]')
u = 1 + x + 2*y
f = - sym.diff(q(u)*sym.diff(u, x), x) - sym.diff(q(u)*sym.diff(u,y), y)
f = sym.simplify(f)
u_code = sym.printing.ccode(u)
f_code = sym.printing.ccode(f)
print('u =', u_code)
print('f =', f_code)

# create mesh and define function space
mesh = fs.UnitSquareMesh(8,8)
V = fs.FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = fs.Expression(u_code, degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = fs.DirichletBC(V, u_D, boundary)

# Define variational problem
u = fs.Function(V)
v = fs.TestFunction(V)
f = fs.Expression(f_code, degree=2)
F = q(u)*fs.dot(fs.grad(u), fs.grad(v))*fs.dx - f*v*fs.dx
a, L = fs.lhs(F), fs.rhs(F)

# Create VTK file for saving solution
vtkfile = fs.File('solution.pvd')

#Compute Solution
fs.solve(F == 0, u, bc)

#plot solution
fs.plot(u)

# Compute maximum error at vertices. This computation illustrates
# an alternative to using compute_vertex_values as in poisson.py.
u_e = fs.interpolate(u_D, V)
error_max = np.abs(u_e.vector() - u.vector()).max()
print('error_max = {:.5e}'.format(error_max))