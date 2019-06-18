#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 00:20:55 2019

@author: daniel
"""

"""
FEniCS tutorial demo program: Incompressible Navier-Stokes equations
for channel flow (Poisseuille) on the unit square using the
Incremental Pressure Correction Scheme (IPCS).

  u' + u . nabla(u)) - div(sigma(u, p)) = f
                                 div(u) = 0
"""

import fenics as fs
import numpy as np
import matplotlib.pyplot as plt

# scaled variables
T = 10.0                # final time
num_steps = 500         # number of time steps
dt = T / num_steps      # time step size
mu = 1                  # kinematic viscosity
rho = 1                 # density

# Create mesh and define function spaces SEPARATELY for pressure and velocity spaces. 
mesh = fs.UnitSquareMesh(16, 16)
V = fs.VectorFunctionSpace(mesh, 'P', 2)    # velocity space, a vector space
Q = fs.FunctionSpace(mesh, 'P', 1)          # pressure space, a scalar space

# Define trial and test functions
u = fs.TrialFunction(V)     # in velocity space
v = fs.TestFunction(V)
p = fs.TrialFunction(Q)     # in pressure space
q = fs.TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = fs.Function(V)    # u^(n)
u_  = fs.Function(V)    # u^(n+1)
p_n = fs.Function(Q)    # p^(n)
p_  = fs.Function(Q)    # p^(n+1)

# Define boundaries
inflow  = 'near(x[0], 0)'
outflow = 'near(x[0], 1)'
walls   = 'near(x[1], 0) || near (x[1], 1)'

# Define boundary bonditions
bcu_noslip  = fs.DirichletBC(V, fs.Constant((0,0)), walls)
bcp_inflow  = fs.DirichletBC(Q, fs.Constant(8), inflow)
bcp_outflow = fs.DirichletBC(Q, fs.Constant(0), outflow)
bcu = [bcu_noslip]
bcp = [bcp_inflow, bcp_outflow]

# Constants used to define variational problems
U   = 0.5*(u_n + u)
n   = fs.FacetNormal(mesh)
f   = fs.Constant((0, 0))
k   = fs.Constant(dt)
mu  = fs.Constant(mu)
rho = fs.Constant(rho)

# Define strain-rate tensor
def epsilon(u):
    return fs.sym(fs.nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*fs.Identity(len(u))

# Define variational problem for step 1 (Tentative velocity step)
F1 = rho*fs.dot((u - u_n) / k, v)*fs.dx + \
    rho*fs.dot(fs.dot(u_n, fs.nabla_grad(u_n)), v)*fs.dx \
  + fs.inner(sigma(U, p_n), epsilon(v))*fs.dx \
  + fs.dot(p_n*n, v)*fs.ds - fs.dot(mu*fs.nabla_grad(U)*n, v)*fs.ds \
  - fs.dot(f, v)*fs.dx
a1 = fs.lhs(F1)
L1 = fs.rhs(F1)

# Define variational problem for step 2 (Pressure correction step)
a2 = fs.dot(fs.nabla_grad(p),   fs.nabla_grad(q))*fs.dx
L2 = fs.dot(fs.nabla_grad(p_n), fs.nabla_grad(q))*fs.dx \
        - (1/k)*fs.div(u_)*q*fs.dx

# Define variational problem for step 3 (Velocity correction step)
a3 = fs.dot(u,  v)*fs.dx
L3 = fs.dot(u_, v)*fs.dx - k*fs.dot(fs.nabla_grad(p_ - p_n), v)*fs.dx

# Assemble matrices
A1 = fs.assemble(a1)
A2 = fs.assemble(a2)
A3 = fs.assemble(a3)

# Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# Time-stepping
t = 0
for n in range(num_steps):
    
    # Update current time
    t += dt
    
    # Step 1: Tentative velocity step
    b1 = fs.assemble(L1)
    [bc.apply(b1) for bc in bcu]
    fs.solve(A1, u_.vector(), b1)
    
    # Step 2: Pressure correction step
    b2 = fs.assemble(L2)
    [bc.apply(b2) for bc in bcp]
    fs.solve(A2, p_.vector(), b2)
    
    # Step 3: Velocity correction step
    b3 = fs.assemble(L3)
    fs.solve(A3, u_.vector(), b3)
    
    # Plot solutions
    fs.plot(u_)
    plt.pause(0.1)
    
    # Compute error
    u_e = fs.Expression(('4*x[1]*(1.0 - x[1])', '0'), degree=2)
    u_e = fs.interpolate(u_e, V)
    error = np.abs(u_e.vector() - u_.vector()).max()
    print('t = %.2f: error = %.3g' % (t, error))
    print('max u:', u_.vector().max())
    
    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)