#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 20:58:46 2019

@author: daniel
"""

"""
FEniCS tutorial demo program: Linear elastic problem.

  -div(sigma(u)) = f
  sigma(epsilon) = lambda*tr(epsilon)*I + 2*mu*epsilon
  epsilon(u)     = (1/2) * (grad(u) + transpose(grad(u)))

The model is used to simulate an elastic beam clamped at
its left end and deformed under its own weight.
"""

import fenics as fs
import matplotlib.pyplot as plt

# scaled variables
L = 1
W = 0.2
mu = 1
rho = 1
delta = W/ L
gamma = 0.4 * delta**2
beta = 1.25
lambda_ = beta
g = gamma

# Create mesh and define function space
mesh = fs.BoxMesh(fs.Point(0, 0, 0), fs.Point(L, W, W), 10, 3, 3)
V = fs.VectorFunctionSpace(mesh, 'P', 1)

# define boundary condition
tol = 1e-14

def clamped_boundary(x, on_boundary):
    return on_boundary and x[0] < tol

bc = fs.DirichletBC(V, fs.Constant((0, 0, 0)), clamped_boundary)

# define strain and stress
def epsilon(u):
    return 0.5 * (fs.nabla_grad(u) + fs.nabla_grad(u).T)

def sigma(u):
    return lambda_*fs.nabla_grad(u)*fs.Identity(d) + 2*mu*epsilon(u)

# Define variational problem
u = fs.TrialFunction(V)
d = u.geometric_dimension()
v = fs.TestFunction(V)
f = fs.Constant((0, 0, -rho*g))
T = fs.Constant((0, 0, 0))
a = fs.inner(sigma(u), epsilon(v))*fs.dx
L = fs.dot(f, v)*fs.dx + fs.dot(T, v)*fs.ds

# Compute solution
u = fs.Function(V)
fs.solve(a == L, u, bc)

# Plot solution
plt.figure()
fs.plot(u, title='Displacement', mode='displacement')

# Plot stress
s = sigma(u) - (1./3)*fs.tr(sigma(u))*fs.Identity(d) # deviatoric stress
von_Mises = fs.sqrt(3./2*fs.inner(s, s))
V = fs.FunctionSpace(mesh, 'P', 1)
von_Mises = fs.project(von_Mises, V)
plt.figure()
fs.plot(von_Mises, title='Stress intensity')

# Compute magnitude of displacement
u_magnitude = fs.sqrt(fs.dot(u, u))
u_magnitude = fs.project(u_magnitude, V)
plt.figure()
fs.plot(u_magnitude, title='Displacement magnitude')
print('min/max u:', u_magnitude.vector().min(), u_magnitude.vector().max())

# Save solution to file in VTK format
fs.File('results/displacement.pvd') << u
fs.File('results/von_mises.pvd') << von_Mises
fs.File('results/magntiude.pvd') << u_magnitude