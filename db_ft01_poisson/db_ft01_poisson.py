"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.

  -Laplace(u) = f    in the unit square
            u = u_D  on the boundary

  u_D = 1 + x^2 + 2y^2
    f = -6
"""

import fenics as fs
import numpy as np
import matplotlib.pyplot as plt

# create mesh and define function space
mesh = fs.UnitSquareMesh(8, 8)
V = fs.FunctionSpace(mesh, 'P', 1)

# define boundary condition
u_D = fs.Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)

def boundary(x, on_boundary):
	return on_boundary
	
bc = fs.DirichletBC(V, u_D, boundary)

# define variational problem
u = fs.TrialFunction(V)
v = fs.TestFunction(V)
f = fs.Constant(-6.0)
a = fs.dot(fs.grad(u), fs.grad(v)) * fs.dx
L = f * v * fs.dx

# compute solution
u = fs.Function(V)
fs.solve(a == L, u, bc)

# plot solution and mesh
fs.plot(u)
fs.plot(mesh)

# save solution to file in VTK format
vtkfile = fs.File('solution.pvd')
vtkfile << u

# compute error in L2 norm
error_L2 = fs.errornorm(u_D, u, 'L2')

# #compute maximum error at vertices
vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

# print errors
print('error_L2 =', error_L2)
print('error_max = ', error_max)

plt.savefig('newtest',Format='png')