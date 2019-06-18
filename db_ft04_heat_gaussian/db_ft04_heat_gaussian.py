#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 19:37:24 2019

@author: daniel
"""

"""
Gonna solve for the diffusion of a Gaussian hill

    u' = Laplacian(u) + f   in a square domain [-2,-2]x[-2,2]
    
    u  = u_D                on the boundary (boundary condition).
                            Chosen to be zero. Dirichlet boundary.
    u  = u_0                at t = 0 (initial condition).
                            Chosen to be a Gaussian hill. exp(-a*x**2 - a*y**2)

"""

import fenics as fs
import matplotlib.pyplot as plt

T = 2.0            # final time
num_steps = 50     # number of time steps
dt = T / num_steps # time step size

# create mesh and define function space
nx = ny = 30
mesh = fs.RectangleMesh(fs.Point(-2,-2), fs.Point(2,2), nx, ny)
V = fs.FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = fs.Constant(0)

def boundary(x, on_boundary):
    return on_boundary

bc = fs.DirichletBC(V, u_D, boundary)

# Define initial value
u_0 = fs.Expression('exp(-a * pow(x[0], 2) - a * pow(x[1], 2))', degree=2, a=5)
u_n = fs.interpolate(u_0, V)

# Define variational problem
u = fs.TrialFunction(V)
v = fs.TestFunction(V)
f = fs.Constant(0)

F = u * v * fs.dx + dt * fs.dot(fs.grad(u), fs.grad(v)) *fs.dx - (u_n + dt * f) * v * fs.dx
a, L = fs.lhs(F), fs.rhs(F)

# Create VTK file for saving solution
vtkfile = fs.File('heat_gaussian/solution.pvd')

# Time-stepping
u = fs.Function(V)
t = 0
for n in range(num_steps):

    # Update current time
    t += dt

    # Compute solution
    fs.solve(a == L, u, bc)

    # Save to file and plot solution
    vtkfile << (u,t)
    fs.plot(u)

    # Update previous solution
    u_n.assign(u)
    
    plt.pause(0.1)

# Hold plot
#interactive()