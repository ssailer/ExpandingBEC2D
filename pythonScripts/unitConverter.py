
#!/usr/bin/env python
import numpy as np
import math
# import h5py
# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():
	L = 100 * 10**(-6)
	N = 2048
	number = 2 * 10**5
	
	m = 87 * 1.66 * 10**(-27)
	hbar = 1.054 * 10**(-34)
	a0 = 5.29 * 10**(-11)
	# m = 1
	# hbar = 1
	g =  4 * math.pi * 106 * a0 * hbar * hbar / m
	cg = g * m / ( hbar * hbar)
	U = g / hbar
	a = 10**-6
	cnumber = a * a * number	
	xmin = 50 * 10**(-6)
	cxmin = xmin / a
	omega_x = 180 * 2 * math.pi
	omega_y = 207 * 2 * math.pi
	comega_x = omega_x * m * a * a / hbar
	comega_y = omega_y * m * a * a / hbar
	runtime = 10**(-5)	
	C_t = hbar / (m * a * a)
	C_w = m * a * a / hbar
	cruntime = runtime / C_t

	mu = math.sqrt(0.75 * m * g * omega_x * omega_y * number / 2)
	# mu = 10**-6
	Rx = math.sqrt( 2 * mu / (m * omega_x * omega_x))
	Ry = math.sqrt( 2 * mu / (m * omega_y * omega_y))

	solvedMu = 0.5 * m * omega_y * omega_y * 40 * 40 * 10**-12
	solvedG = solvedMu * solvedMu * 8 / ( 3 * m * omega_y * omega_x * number)
	usedG = solvedG * m / hbar / hbar

	Rxsolved = math.sqrt( 2 * solvedMu / (m * omega_x * omega_x))
	print Rxsolved
	print "alpha = ", 16 * number * hbar * hbar / (Rxsolved * m * m )
	print "beta = ", 4 * ( hbar * 1 / m)**2
	print "gamma = ", Rxsolved**2

	print "a =", a
	print "cxmin =", cxmin
	print "comega_x =", comega_x, "comega_y =", comega_y
	print "Number of particles =", number
	print "if runtime =",runtime, "then cruntime =", cruntime
	print "g =", g / hbar
	print "cg =", cg
	# print "cnumber =", cnumber
	print "C_t =", C_t, "C_w =", omega_x * C_w, omega_y * C_w
	print "delta t ", 10**(-1) / C_t
	print "U = ", U
	print "mu =", mu, " solvedMu =", solvedMu, "solvedG =", solvedG * m / hbar / hbar
	print "Rx in micrometers=", Rx * 10**6
	print "Ry in micrometers=", Ry * 10**6
	print "usedG =", usedG, "usedOmegaX =", comega_x, "usedOmegaY =", comega_y, "usedDeltaT", cruntime






if __name__ == '__main__':
   	main()

