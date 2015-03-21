To plot a lattice instantiate a plot object, the parameters include: 
	N- size of lattice
	B- strength of magnetic field (default is zero)
	start- low(cold) or high(hot) start (default is low)
	inc- size of increments in plots (default is 0.01)
	x0- starting point of plots (default is 1)
	x1- final point of plots (default is 5)
	steps-number of steps (default is 50000)
	T - Temperature (default is 1) 

Some recipes to showcase how the program works:

Plots specific heat capacity with a cold start and with increments of 0.1 for a 10x10 lattice:

	P=plots(N=10,start='Low',inc=0.1)
	P.spec_heat()
	P.show()

Plots magnetization with with a hot start increments of 0.1 for a 15x15 lattice

	P=plots(N=15,start='Low',inc=0.1)
	P.mag()
	P.show()

Plots lattice with a hot start, with and a 150x150 lattice and equilibrated at a temperature of 2

	P=plots(N=150,start='High',inc=0.1,T=2)
	P.lattice()
	P.show()

Plots lattice with a hot start and with and a 100x100 lattice, with an external magnetic field of 1, and ran for 60000 steps

	P=plots(N=100,start='High',B=1,steps=60000)
	P.lattice()
	P.show()