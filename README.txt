--------------------------------------------------------------------------------------------------
About FEniCS
--------------------------------------------------------------------------------------------------
FEniCS is a convenient open-source package for Python and C++ that can be used to solve PDEs using the finite element method. (Details on how the finite element method works in general can be found at https://en.wikipedia.org/wiki/Finite_element_method.) While implementing the finite element method is typically a sordid affair, where the user has to define finite meshes and implement a solver algorithm themself, FEniCS allows the user to easily define a mesh domain, boundary conditions, and choose from a variety of pre-selected solving schemes to quickly model physical systems described by a system of PDEs.

FEniCS is developed by the FEniCS Project, which offers a great deal of exposition into using FEniCS, and the mathematical basis behind FEniCS and the finite element method. Books published by the FEniCS Project can be found at https://fenicsproject.org/book/ and https://fenicsproject.org/tutorial/. The second link provides a quick tutorial in getting started with FEniCS, whereas the first is a far more detailed text on FEniCS and the finite element method.

--------------------------------------------------------------------------------------------------
A Note on the Calculus of Variations
--------------------------------------------------------------------------------------------------
FEniCS interprets problems using the calculus of variations; that is, PDEs cannot be simply be input to a FEniCS in their standard form. To solve a PDE using FEniCS, the PDE must be converted into what is commonly called a "weak form." The oversimplified version of this process is that if you're attempting to solve for a function u (aka the "trial function"), you multiply both sides of the PDE by a so-called "test function" and integrate. This projects u into a different function space (defined by your finite element mesh), which is much easier for FEniCS to interpret and solve for. A much more mathematically rigorous explanation can be found in Chapters 1 and 2 of the FEniCS tutorial https://fenicsproject.org/tutorial/, in the FEniCS book https://fenicsproject.org/book/, or in a textbook on the calculus of variations of your choice.)

--------------------------------------------------------------------------------------------------
My Goals
--------------------------------------------------------------------------------------------------
I'm using this space to keep track of my work through the FEniCS tutorials found at https://fenicsproject.org/tutorial/. This is useful for research I'm performing with Dr. John Colton of the BYU Physics Department. We wish to use FEniCS to form a training data set of resonance characteristics of electromagnetic resonant cavities that can be later used by a neural network to predict the same characteristics of a resonant cavity of arbitrary (likely cylindrically-symmetric) geometry. While FEniCS can be a relatively fast tool (that can be potentially run on a laptop), we envision implementing our neural network as an open-source web-based app.

--------------------------------------------------------------------------------------------------
On Installation
--------------------------------------------------------------------------------------------------
Installing FEniCS ended up being a lot more tricky than I had anticipated. Conda distributions of FEniCS only supported deployment on macOS and Linux platforms, and the home computer where I was doing most of my work only had Windows 10 installed.

For a Windows 10 installation, the best option I found was using a Docker container created by the FEniCS project. Installing Docker (which, if you don't know, is a containerization platform) is a pain on Windows, but doable using the tutorial at https://docs.docker.com/docker-for-windows/install/. I didn't like using the Docker container because file management was so difficult; the Docker container is essentially a walled garden containing an independent installation of Python and FEniCS, and moving development files into the Docker was annoying. Docker installations also don't support dynamic figure generation; you have to save whatever figure you're making to a separate file, which is annoying while debugging.

Eventually using the Docker environment was so irritating that I decided to install Ubuntu for dual-boot on my home computer and just used the Conda distribution of FEniCS. I have a personal preference for coding using Spyder, so making a FEniCS environment in Anaconda and being able to use Spyder, Jupyter Notebooks, etc. was my preferred outcome.

It is also possible to use the Google Colaboratory to run FEniCS, but my personal preference is to use the Conda distribution on my Ubuntu 18.04LTS installation.

--------------------------------------------------------------------------------------------------
Technical Notes Affecting All Tutorials
--------------------------------------------------------------------------------------------------
	* Many of the tutorials compute the error of the solution found by FEniCS to a known analytic solution. The code published in the tutorial is incorrect as it references a structure that doesn't exist in FEniCS Function objects. Much of the code references variables of the type
		.vector().array()
belonging to Function objects when calculating error, but it should be just variables of the type
		.vector()
	* The Conda installation (or at least Spyder) doesn't seem to support the fenics.interactive() command, which is supposed to make figures generated by fenics.plot() more interactive. I'm not sure what's broken about it. Docker environments didn't support dynamic figure generation, so this command obviously didn't do anything there, either.
