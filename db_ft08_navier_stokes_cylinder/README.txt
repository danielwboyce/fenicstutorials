--------------------------------------------------------------------------------------------------
About This Problem
--------------------------------------------------------------------------------------------------
This problem is an extension of the previous one; instead of looking at incompressible fluid flow through a simple channel this time we will be looking at incompressible fluid flow through a channel that has a cylindrical obstacle in it.

Incompressible Navier-Stokes equations for flow around a cylinder on the unit square using the Incremental Pressure Correction Scheme (IPCS).

	  u' + u . nabla(u)) - div(sigma(u, p)) = f
	                                 div(u) = 0

--------------------------------------------------------------------------------------------------
Method
--------------------------------------------------------------------------------------------------
This problem shows us that once you have a problem set up in FEniCS it's an almost trivial matter to change the geometry of the problem to be more complicated. All we had to do in this problem to add the cylindrical barrier is to add about five lines of code defining the new domain on which to define a function space and to define a boundary condition around the obstacle.

The tutorial also tried to show us how to use a progress bar to show how much of the computation has been completed, but the provided code for that part of the program didn't compile. I was able to fix it, but even then the progress bar didn't appear and I wasn't terribly interested in troubleshooting that.

Running the computation for this problem took about a while, 1066.311 seconds to be exact. It ran much faster without attempting to display two graphs at every iteration, and so I saved the files into XDMF/HDF5 format, which can be opened with ParaView.
