--------------------------------------------------------------------------------------------------
About This Problem
--------------------------------------------------------------------------------------------------
Here we're solving one of the most basic PDEs there is: Poisson's equation
		-Laplace(u) = f		in the unit area
			  u = u_D	on the boundary
with boundary condition
		u_D = 1 + x^2 + 2*y^2
		  f = -6
on the unit square.

--------------------------------------------------------------------------------------------------
Method
--------------------------------------------------------------------------------------------------
From this problem and solution we glean a basic order of things in setting up a FEniCS (a.k.a. fenics) program:
	1) Create a mesh and define function space
	2) Define boundary conditions
	3) Define the variational problem (from the weak form of the PDE)
	4) Compute solution
	Optional:
	5) Plot solution and mesh
	6) Save solution to a file
	7) Compute error

Now is also a good time to take note of a few important points:
	* The tutorial imports fenics wholesale, but I think that's bad form and I don't like not knowing where my functions are being called from, so I just use
		import fenics as fs
	* There's a function in fenics for designing most meshes; you just have to say what kind of mesh it should be, its size, etc.
	* Fenics has an Expression function which takes a string and builds a numerical expression using a built-in interpreter.
	* The boundary(x, on_boundary) function seems a little weird to me. The FEniCS tutorial says that apparently the x object already knows if its on the boundary (stored in the Boolean variable on_boundary), but that is assigned in the DirichletBC function. Unfortunately this function seems a little like a magic function that's just there to make the code work, I'll just have to deal with it.
	* TrialFunction and TestFunction are distinct object types; it's important to make sure I put the right one in the right place.
