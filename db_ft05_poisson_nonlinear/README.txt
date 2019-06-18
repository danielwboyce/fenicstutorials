--------------------------------------------------------------------------------------------------
About This Problem
--------------------------------------------------------------------------------------------------
With this problem we begin our foray into solving nonlinear PDEs. In particular here we'll be solving the nonlinear Poisson's equation:
		-div(q(u) * grad(u)) = f    in the unit square
with boundary conditions
   		u    = u_D 		    on the boundary (boundary condition).
    		u    = u_0                  at t = 0 (initial condition).
		u_0  = exp(-a*x**2 - a*y**2)
		q(u) = 1 + u**2
Because this problem is nonlinear, we need to use some different methods to define the variational problem, though after the initial problem is defined things remain largely the same. We are also manufacturing the solution to this problem, so at the end we should find
		u   = 1 + x + 2*y
		u_D = 1 + x + 2*y

--------------------------------------------------------------------------------------------------
Method
--------------------------------------------------------------------------------------------------
A first, easy, note on this problem is that our nonlinear coefficent q(u) has to be defined independently of the other parts of the problem. We see that in lines 25-27.

After that definition, we have to use the SymPy package to fully define our variational problem, which can be a little tricky. The first step in using SymPy is defining our symbolic variables. We do this here with the code
		x, y = sym.symbols('x[0], x[1]')
Intuition would have told us to define x and y
		x, y = sym.symbols('x, y')
but by using 'x[0], x[1]' instead of 'x, y' we won't have to make any changes to our code when it gets exported to C++.

Our function f is given to be equal to -div(q(u) * grad(u)), but we have to write this in SymPy syntax. The task isn't too difficult here since we're using Cartesian two-dimensional coordinates, and the code to fulfill this task is quite readable:
		f = - sym.diff(q(u)*sym.diff(u, x), x) - sym.diff(q(u)*sym.diff(u,y), y)
		f = sym.simplify(f)
We just have to note to find the divergence of the gradient for each variable one-by-one. Other problems will be expressed differently with different tricks; it's likely that some algebra will be necessary to find the correct way to express f in different problems using SymPy.

After defining f and u using SymPy, we can use the function sympy.printing.ccode(u) to transform the function u into C++ syntax, which will of course be useful in defining the problem later with fenics.Expression().

Beyond the changes made to defining the variational problem with SymPy, the solution to this problem isn't any different to the previous ones we've solved.

