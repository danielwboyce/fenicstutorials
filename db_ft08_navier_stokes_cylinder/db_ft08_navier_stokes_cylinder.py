#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:03:19 2019

@author: daniel
"""


"""
FEniCS tutorial demo program: Incompressible Navier-Stokes equations
for flow around a cylinder on the unit square using the
Incremental Pressure Correction Scheme (IPCS).

  u' + u . nabla(u)) - div(sigma(u, p)) = f
                                 div(u) = 0
"""

import fenics as fs
import mshr
import matplotlib.pyplot as plt

import time
start_time = time.time()

# scaled variables
T = 5.0                 # final time
num_steps = 5000       # number of time steps
dt = T / num_steps      # time step size
mu = 0.001              # kinematic viscosity
rho = 1                 # density

# Create mesh
channel = mshr.Rectangle(fs.Point(0,0), fs.Point(2.2, 0.41))
cylinder = mshr.Circle(fs.Point(0.2, 0.2), 0.05)
domain = channel - cylinder
mesh = mshr.generate_mesh(domain, 64)

# Define function spaces separately for velocity and pressure space separately
V = fs.VectorFunctionSpace(mesh, 'P', 2)    # velocity space, a vector space
Q = fs.FunctionSpace(mesh, 'P', 1)          # pressure space, a scalar space

# Define boundaries
inflow  = 'near(x[0], 0)'
outflow = 'near(x[0], 2.2)'
walls   = 'near(x[1], 0) || near (x[1], 0.41)'
cylinder = 'on_boundary && x[0]>0.1 && x[0]<0.3 && x[1]>0.1 && x[1]<0.3'

# Define inflow profile
inflow_profile = ('4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)', '0')

# Define boundary bonditions
bcu_inflow   = fs.DirichletBC(V, fs.Expression(inflow_profile, degree=2), inflow)
bcu_walls    = fs.DirichletBC(V, fs.Constant((0, 0)), walls)
bcu_cylinder = fs.DirichletBC(V, fs.Constant((0, 0)), cylinder)
bcp_outflow   = fs.DirichletBC(Q, fs.Constant(0), outflow)
bcu = [bcu_inflow, bcu_walls, bcu_cylinder]
bcp = [bcp_outflow]

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

# Create XDMF files for visualization output
xdmffile_u = fs.XDMFFile('navier_stokes_cylinder/velocity.xdmf')
xdmffile_p = fs.XDMFFile('navier_stokes_cylinder/pressure.xdmf')

# Create time series (for use in reaction_system.py)
timeseries_u = fs.TimeSeries('navier_stokes_cylinder/velocity_series')
timeseries_p = fs.TimeSeries('navier_stokes_cylinder/pressure_series')

# Save mesh to file (for use in reaction_system.py)
fs.File('navier_stokes_cylinder/cylinder.xml.gz') << mesh

# Create progress bar
progress = fs.Progress('Time-stepping', num_steps)
#fs.set_log_level(0)

# Time-stepping
t = 0
for n in range(num_steps):
    
    # Update current time
    t += dt
    
    # Step 1: Tentative velocity step
    b1 = fs.assemble(L1)
    [bc.apply(b1) for bc in bcu]
    fs.solve(A1, u_.vector(), b1, 'bicgstab', 'hypre_amg')
    
    # Step 2: Pressure correction step
    b2 = fs.assemble(L2)
    [bc.apply(b2) for bc in bcp]
    fs.solve(A2, p_.vector(), b2, 'bicgstab', 'hypre_amg')
    
    # Step 3: Velocity correction step
    b3 = fs.assemble(L3)
    fs.solve(A3, u_.vector(), b3, 'cg', 'sor')
    
    # Plot solutions
    #plt.figure(1)
    #fs.plot(u_, title='Velocity')
    #plt.figure(2)
    #fs.plot(p_, title='Pressure')
    #plt.pause(0.001)
    #plt.clf()
    
    # Save solution to file (XDMF/HDF5)
    xdmffile_u.write(u_, t)
    xdmffile_p.write(p_, t)
    
    # Save nodal values to file
    timeseries_u.store(u_.vector(), t)
    timeseries_p.store(p_.vector(), t)
    
    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)
    
    # Update progress bar
    progress += 1
    print('u max:', u_.vector().max())
    print('The program is at {:1.3f}'.format(t), ' seconds of ', T, '.')
    print('This program has been running for ', time.time() - start_time, ' seconds.')
    
