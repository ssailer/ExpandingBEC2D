
#!/usr/bin/env python
import numpy as np
import math
# import h5py
# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():
	# L = 150 * 10**-6
	# N = 2048
	# number = 2 * 10**5
	
	m = 87 * 1.66 * 10**(-27)
	hbar = 1.054 * 10**(-34)
	# a0 = 5.29 * 10**(-11)
	# # m = 1
	# # hbar = 1
	# g =  4 * math.pi * 106 * a0 * hbar * hbar / m
	# cg = g * m / ( hbar * hbar)
	# U = g / hbar
	# a = 10**-6
	# cnumber = a * a * number	
	# xmin = 75 * 10**-6
	# cxmin = xmin / a
	# omega_x = 180 * 2 * math.pi
	# omega_y = 207 * 2 * math.pi
	# comega_x = omega_x * m * a * a 
	# comega_y = omega_y * m * a * a
	# runtime = 10**(-5)	
	# C_t = 1/ (m * a * a)
	# C_w = m * a * a
	# cruntime = runtime / C_t

	# mu = math.sqrt(0.75 * m * g * omega_x * omega_y * number / 2)
	# # mu = 10**-6
	# Rx = math.sqrt( 2 * mu / (m * omega_x * omega_x))
	# Ry = math.sqrt( 2 * mu / (m * omega_y * omega_y))

	# solvedMu = 0.5 * m * omega_y * omega_y * 40 * 40 * 10**-12
	# solvedG = solvedMu * solvedMu * 8 / ( 3 * m * omega_y * omega_x * number)
	# usedG = solvedG * m / hbar / hbar

	# Rxsolved = math.sqrt( 2 * solvedMu / (m * omega_x * omega_x))
	# print Rxsolved
	# print "alpha = ", 16 * number * hbar * hbar / (Rxsolved * m * m )
	# print "beta = ", 4 * ( hbar * 1 / m)**2
	# print "gamma = ", Rxsolved**2

	# print "a =", a
	# print "cxmin =", cxmin
	# print "comega_x =", comega_x, "comega_y =", comega_y
	# print "Number of particles =", number
	# print "if runtime =",runtime, "then cruntime =", cruntime
	# print "g =", g / hbar
	# print "cg =", cg
	# # print "cnumber =", cnumber
	# print "C_t =", C_t, "C_w =", omega_x * C_w, omega_y * C_w
	# print "delta t ", 10**(-1) / C_t
	# print "U = ", U
	# print "mu =", mu, " solvedMu =", solvedMu, "solvedG =", solvedG * m / hbar / hbar
	# print "Rx in micrometers=", Rx * 10**6
	# print "Ry in micrometers=", Ry * 10**6
	# print "usedG =", usedG, "usedOmegaX =", comega_x, "usedOmegaY =", comega_y, "usedDeltaT", cruntime


	# N = 2 * 10**5
	# # omega_x = 180 * 2 * math.pi
	# # omega_y = 207 * 2 * math.pi
	# mu2d = 0.5 * m * comega_x**2 * cxmin**2 / 4 
	# # mu2d = math.sqrt( (3/8) * m * g * omega_x * omega_y * N)
	# # g2d = mu2d**2 * 8 / (3 * m * comega_x * comega_y * N)
	# g2d = (2/3) * comega_x**3 * (cxmin/2)**4 / (comega_y * N)
	# Rx = math.sqrt( 2 * mu2d / (m * comega_x**2))
	# Ry = math.sqrt( 2 * mu2d / (m * comega_y**2))
	# print N, comega_x, comega_y, mu2d
	# print Rx , Ry 
	# print g2d

	# print "---------------------------------------------------"
	# print "delta t = ", 10**-05 * C_t
	# print "omega_x = ", comega_x 
	# print "omega_y = ", comega_y 
	# print "g2d = ", g2d 


	# g2d = (2/3) * (omega_x**3 / omega_y) * Rx**4 / N

	# lengthx
	# lengthy
	# radiusx
	# radiusy

	
	# a = 10**-5
	# N = 2 * 10**5
	# cN = N / a**2
	# cN = N
	# deltaT = 1.0 * 10**-5
	# cdeltaT = deltaT / (hbar / (m * a * a))


	# # m = 1.
	# Rx = 32.12 * 10**-6
	# Ry = 66.495/2 * 10**-6
	# cRx = Rx / a
	# cRy = Ry / a
	# omega_x = 200 * 2 * math.pi
	# omega_y = 150 * 2 * math.pi
	# comega_x = omega_x * m * a * a / hbar
	# comega_y = omega_y * m * a * a / hbar
	# g2d1 = (2.0/3.0) * m * (omega_x**3) * (Rx**4) / (omega_y * N)
	# g2d2 = (2.0/3.0) * m * (omega_y**3) * (Ry**4) / (omega_x * N)
	# cg2d1 = (2.0/3.0) * m * (comega_x**3) * (cRx**4) / (comega_y * cN)
	# cg2d2 = (2.0/3.0) * m * (comega_y**3) * (cRy**4) / (comega_x * cN)
	# N0 = 2 * N / ( math.pi * cRx * cRy)
	# print N, cN, N0
	# print deltaT, cdeltaT
	# print comega_x, comega_y
	# print "#1-----------"
	# print g2d1, g2d2, Rx, Ry
	# print cg2d1, cg2d2, cRx, cRy

	# # g2d = 4 *math.pi * hbar * 5.1 * 10**-35 / m
	
	# g2d = (g2d1) #5.26777092251e-71 #15.0 #1.0e-9 
	# # cg2d =  15.0 #(cg2d2 + cg2d1)/2
	# cg2d = g2d * m / ( hbar **2 )
	# # cN = 1000
	# # cg2d = 11
	# # comega_x = 172
	# # comega_y = 178
	# # N = 1000 #200000
	# # # omega_x = 40.0
	# # # omega_y = 80.0
	# # # comega_x = omega_x / (m * a**2)
	# # # comega_y = omega_y / (m * a**2)
	# mu2d = math.sqrt( 3.0 * m * g2d * omega_x * omega_y * N / 8.0)
	# Ry = math.sqrt(2.0 * mu2d / (m * omega_y**2))
	# Rx = math.sqrt(2.0 * mu2d / (m * omega_x**2))

	# cmu2d = math.sqrt( 3.0  * cg2d * comega_x * comega_y * cN / 8.0)
	# cRy = math.sqrt(2.0 * cmu2d / ( comega_y**2))
	# cRx = math.sqrt(2.0 * cmu2d / ( comega_x**2))
	# # cRx = Rx * a;
	# # cRy = Ry * a;

	# # print N
	# # print omega_x, omega_y
	# # # print comega_x, comega_y
	# print "#2-----------"
	# print g2d, Rx , Ry
	# # print cg2d / 4 / math.pi
	# print "                 "
	# print "   NEW UNITS"
	# print "                 "
	# print "potential x", "potential y"
	# print comega_x, comega_y
	# print "g", "Rx", "Ry"
	# print cg2d, cRx , cRy
	# print "number"
	# print cN
	# print "                 " 
	# # # print 10**-5 / (m * a**2)
	# # print (2.0/3.0) * (omega_x**3 / omega_y) * 2.0**Rx / N
	# Ry = 30.0
	# Rx = 25.0
	# Nv = 1000

	# alpha = g2d * 4.0 * N / (math.pi * m *  Ry)
	# beta = 4 * hbar**2 / ( m**2 )
	# delta = Ry**2
	# print "#3-----------"
	# print "alpha =", alpha
	# print "beta =", beta
	# print "delta =", delta
	# print "             "
	# print m * a * a / hbar

	N = 2.0 * 10**5
	hbar = hbar * 10**12
	# hbar = 1
	# m = 1
	omega_z = 10 * 2 * math.pi
	As = 5.8 * 10**-3
	g2D = As * N * math.sqrt(8.0 * math.pi * omega_z * hbar**3 / m)
	
	
	Ag = 200.0 / 2048.0
	print "Rx, Ry", 676.007  * Ag, 455.29 * Ag
	OmegaG = hbar / ( m * Ag * Ag)
	       
	
	cN = N * Ag * Ag
	omega_x = 20 * 2 * math.pi
	omega_y = 25 * 2 * math.pi
	
	comega_x = omega_x / OmegaG
	comega_y = omega_y / OmegaG

	alpha = 2 * math.pi / 360;
	t = 45
	cRomega_x = comega_x * math.cos(alpha * t) + comega_y * math.sin(alpha * t);
	cRomega_y = - comega_x * math.sin(alpha * t) + comega_y * math.cos(alpha * t);

	deltaT = 1.0 * 10**-6  * OmegaG
	
	cg2D = g2D * Ag * Ag / (hbar * OmegaG)
	print "          "
	print "          "
	print "          "
	print "m / hbar ", m / hbar, 1.0/730.0
	print "Ag = ", Ag
	print "OmegaG = ", OmegaG
	print "N = ", cN
	print "Omega X, Omega Y", comega_x, comega_y
	# print "cRomega X, cRomega Y", cRomega_x, cRomega_y
	print "delta T", deltaT
	# print "g2D = ", g2D
	print "cg2D = ", cg2D
	print "g2D", (math.sqrt(8 * math.pi) * 5.8 / 200 )
	print "test", 8.0 / (math.pi), 15.0 / (4 * math.pi)
	print 2 * N * (0.145 / math.pi)



if __name__ == '__main__':
   	main()

