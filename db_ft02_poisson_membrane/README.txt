--------------------------------------------------------------------------------------------------
About This Problem
--------------------------------------------------------------------------------------------------
Here we're still solving Poisson's equation, but it's ultimately quite similar to the previous problem. We're solving
		-Laplace(w) = p  in the unit circle
         		  w = 0  on the boundary
to pysically model the deflection of a membrane with a load p acting on it. The load function p is a Gaussian function centered at (0, 0.6).

--------------------------------------------------------------------------------------------------
Method
--------------------------------------------------------------------------------------------------
The overall method employed here is the same as in the first tutorial, but because the geometry of the boundary condition and the form of f (or here, p) is different.

First we notice that when creating the mesh and function space, a domain variable is first defined and then used as the basis for the mesh variable. The function space after that is defined in the normal way. Defining a domain before creating a mesh is necessary whenever we are creating a function space that isn't perfectly rectilinear; here the domain is simply circular but as we move into more complex problems were holes are removed from the domain (like in a complex flow problem) where nothing is evaluated, this becomes especially important.

We also notice that even if a function has a constant value (such as the boundary condition here w_D = 0), we cannot simply code in
		w_D = 0,
it is necessary to use the fenics.Constant() command to set w_D, as done in this example

This problem also demonstrates the definition of a real-valued function p, using the fenics.Expression() funciton. This was briefly touched on in the previous problem, but here I will note that Expression() takes a string in correct C++ syntax and backend program (written in C++) to change the string into a workable FEniCS object that can be used when solving a variational problem.

Finally, this tutorial gives some good sample code for saving a plot in a file format that can be read by ParaView, an open-source graphics program (https://paraview.org/download/).
