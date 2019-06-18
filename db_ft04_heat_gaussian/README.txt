--------------------------------------------------------------------------------------------------
About This Problem
--------------------------------------------------------------------------------------------------
This time we are solving the heat equation, where heat the distribution begins as a Gaussian hill
		u'  = Laplace(u) + f  in a square domain [-2,-2]x[2,2]
		u   = u_D             on the boundary
		u   = u_0             at t = 0. (The Gaussian hill)

		u_D = 0
		u_0 = exp(-a*x**2 - a*y**2)
		f   = 0

--------------------------------------------------------------------------------------------------
Method
--------------------------------------------------------------------------------------------------
This problem is very similar to the previous one, but showcases a couple of neat tricks that we'll be able to use when solving other problems.

The first harks back to how we change PDEs to their weak form, the form that FEniCS understands. Normally after multiplying by a test function and integrating, we sort the two sides of the new weak form of the PDE such that one side (that we'll call L) is independent of the trial function and that the other side (that we'll call a) has the other terms that are dependent on the trial function. However, if we have not sorted the terms and have some function F, we can use the FEniCS functions fenics.lhs(F) and fenics.rhs(F) to sort F into a and L. As a matter of convention fenics.lhs() returns a (the function dependent on the trial function) and fenics.rhs() returns L (the function independent of the trial function).

The other main trick has to do with the storage of the solution at each time step. By creating a vtkfile (fenics.File() object) variable and then writing u and t to it at each time step we yield a single .vtk file that contains the solution for *each* time step, and it can even be played as an animation in ParaView (as mentioned before, an open-source 3D renderer that can be obtained at paraview.org/download/).
