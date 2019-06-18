--------------------------------------------------------------------------------------------------
About This Problem
--------------------------------------------------------------------------------------------------
Moving into fluid dynamics, our interests move into solving PDEs for TWO different functions. The incompressible Navier-Stokes equations are dependent on both velocity u and pressure p
		nabla * (u' + u.nabla_grad(u)) - div(sigma(u, p)) = f
					    div(u) = 0
Where sigma(u, p) is the stress tensor, which for a Newtonian fluid is
		sigma(u, p) = 2*mu*epsilon(u) - p*Identity(len(u))
and where epsilon(u) is the strain-rate tensor:
		epsiolon(u) = (1/2) * (grad(u) + transpose(grad(u)) = sym(nabla_grad(u))
The parameter mu is dynamic viscosity.

We will be testing our solution scheme for the Navier-Stokes equations in this part of the tutorial by checking fluid flow through a simple channel (Poisseuille flow).

--------------------------------------------------------------------------------------------------
Method
--------------------------------------------------------------------------------------------------
To solve this problem, we have to implement an algorithm that represents a rather large departure from how we've solved the previous problems, called the Incremental Pressure Correction Scheme (IPCS). THe IPCS includes three steps:
	1) Calculate a tentative velocity u* by advancing the momentum equation nabla * (u' + u.nabla_grad(u)) - div(sigma(u, p)) = f
	2) Use u* to calculate the new pressure p_n
	3) Calculate the corrected velocity u_ using u* and p_n.

The variational formulation of each of these three steps is extremely complex, but the derivation is found in Sect. 3.4.2 (pp. 57-60) of the FEniCS tutorial (https://fenicsproject.org/tutorial/).

In the tutorial text there is also an important discussion of the difference between the FEniCS functions grad() and nabla_grad(). If you take the gradient of a scalar-valued function s, then grad(s) and nabla_grad(s) would evaluate to the same vector. However is you are working with a vector-valued function, the question of what sort of gradient you're taking becomes critical. In FEniCS, grad(u) refers to the matrix pd(u_i)/pd(x_j), but nabla_grad(u) refers to the matrix pd(u_j)/pd(x_i). The use of nabla_grad() is more common than grad() in continuum mechanics, and we will be using nabla_grad() in this problem, but it is essential that you understand the difference between grad() and nabla_grad() and which one you need to be using in your own FEniCS solutions. In general, when the basis of your PDE considers del to be a VECTOR and not an operator, you will use nabla_grad(), nabla_div(), or nabla_curl(), but if you consider del to be an operator you will use the standard grad(), div(), or curl(). A more thorough exposition of the distinction between situations when you should use grad() or nabla_grad() (or the equivalents for div and curl) can be found in the box on pp. 58-59 of the FEniCS tutorial.

Now looking at the code for this solution, there are quite a few things to point out. First is while the base mesh for the different function spaces can be the same, distinct function spaces should be defined for velocity and pressure. Here the space for velocity is a vector space and the space for pressure is a scalar space, but if solving for two different functions defining two different function spaces is good practice. Likewise the trial and test functions for each function must be defined in their respective spaces, as well as functions for solutions at previous and current time steps if your system varies in time.

The definition of boundary conditions are also different from before in that they are defined as strings in C++ syntax. The various types of boundaries are calculated using the FEniCS C++ function near, which checks if the position variable is near (within 3E-16) a given constant. Dirichlet boundary conditions are defined in arrays for easy access later, and then applied using for loops only AFTER the other parts of the variational problem are defined.

Variational problems are also defined separately, one variational problem for each step of IPCS. After finding a and L for each variational problem, we'll slightly change how we solve the variational problems for greater efficiency.

Typically, each time fenics.solve(a, L) is called, several substeps occur:
		A = fenics.assemble(a)
		b = fenics.assemble(L)
		bc.apply(A, b)
		fenics.solve(A, u.vector(), b)
where the overloaded fenics.solve() function solves the matrix equation A.u=b for u. This time, we will be explicitly be assembling our own matrices because the A matrices for each step are time-independent and assembling them in each time step by using the regular fenics.solve() function would waste computational time. The b matrices are depdendent on time, however, so they must be assembled at each time step.

Putting all this together, we have a coherent solution scheme for the incompressible Navier-Stokes equations, though the assemblage of our scheme was quite tedious.
