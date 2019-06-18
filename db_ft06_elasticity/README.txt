--------------------------------------------------------------------------------------------------
About This Problem
--------------------------------------------------------------------------------------------------
Here we wish to investigate a problem of linear elasticity. The PDEs describing this system are nonlinear, as given:
Linear elastic problem.

	  -div(sigma(u)) = f
	  sigma(epsilon) = lambda*tr(epsilon)*I + 2*mu*epsilon
	  epsilon(u)     = (1/2) * (grad(u) + transpose(grad(u)))
We are looking to model a long, rectangular brick that is clamped to the wall at one end and then bent due to the force of gravity. We'll also be scaling the variables to try and simplify the terms that will be included in the variational problem.

--------------------------------------------------------------------------------------------------
Method
--------------------------------------------------------------------------------------------------
While this problem seems rather daunting, in truth a big reason is because of the derived values that we will be calculating. We start the program with our scaled variables that will be used elsewhere, and then define the function space in largely the same way we've done before. The only difference there is that we are defining the space as a *vector* space, which contrasts slightly from what we've done before but really only means that we need to define vectors at the boundaries, etc., instead of scalars.

When looking at defining the boundary condition, we also see that we make a slight change of boundary(), changing it to clamped_boundary(). This new function returns true only if it's at the outer boundary of the domain (like normal) AND if the x-coordinate of its position is within tol away from x = 0, which together is a clever way to check if you're at the right plane or not.

Epsilon and sigma are named quantities for strain and stress, respectively, but the definition of them is the same in principle as the nonlinear coefficient we defined in the last tutorial problem, q(u).

Solving the problem is roughly the same, but then after that we are simply calculating derived values relevant to the physics of elastics, and while it may be a little intimidating, their calculation doesn't really rely on our understanding of FEniCS, but of physics and that's a little outside the scope of what we're worrying about here.
