
#!/usr/bin/env python
import numpy as np
import math
# import h5py
# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():
	L = 1500 * 10**(-6)
	N = 2048
	number = 2 * 10**5
	
	m = 87 * 1.66 * 10**(-27)
	hbar = 1.054 * 10**(-34)
	a0 = 5.29 * 10**(-11)
	# m = 1
	# hbar = 1
	g =  4 * math.pi * 100 * a0 * hbar * hbar / m
	cg = g * 2 * m / ( hbar * hbar)
	a = L / N
	cnumber = a * a * number	
	xmin = 750 * 10**(-6)
	cxmin = xmin / a
	omega_x = 300
	omega_y = 200
	comega_x = omega_x * m * a * a / hbar
	comega_y = omega_y * m * a * a / hbar
	runtime = 1 * 10**(-5)
	cruntime = runtime * m / ( hbar * a * a)
	C_t = hbar / (m * a * a)
	C_w = m * a * a / hbar

	print "a =", a
	print "cxmin =", cxmin
	print "comega_x =", comega_x, "comega_y =", comega_y
	print "Number of particles =", number
	print "if runtime =",runtime, "then cruntime =", cruntime
	print "cg =", cg
	print "cnumber =", cnumber
	print m / hbar
	print "C_t =", C_t
	print "C_w =", omega_x * C_w, omega_y * C_w
	print "delta t ", 10**(-1) / C_t






if __name__ == '__main__':
   	main()

