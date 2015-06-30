
#!/usr/bin/env python
#import h5py
import pandas as pd
import numpy as np
import matplotlib
from math import cos,sin,pi
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():
	# Ag = 50.0 / 500.0
	# m = 87.0 * 1.66e-27;
	# hbar = 1.054e-22;	
	# OmegaG = hbar / ( m * Ag * Ag);

	datafile = 'runObservables/EXP_Observables.dat'

	cols = ["Timestep","Time","X_max","Y_max","Vortexnumber","Alpha","Rx","Ry","R_Ratio","E_Major","E_Minor","E_Major_Angle","E_Minor_Angle","E_Ratio","N","V","N/V","E_kin","N_0"]

	from_data = pd.read_csv(datafile,skiprows=1,sep=',',names=cols)

	dataset1 = from_data['Time']

	# dataset2 = from_data['E_Major']
	# dataset3 = from_data['E_Minor']
	# dataset4 = from_data['E_Major_Angle']
	d1 = from_data['Rx']
	d2 = from_data['Ry']
	d3 = from_data['R_Ratio']

	# d1 = []
	# d2 = []
	# d3 = []

	# for i in range(0,len(dataset1)):
	# 	d1.append(dataset2[i] * cos(dataset4[i] * pi/180) + dataset3[i] * sin(dataset4[i] * pi/180))
	# 	d2.append(dataset2[i] * sin(dataset4[i] * pi/180) + dataset3[i] * cos(dataset4[i] * pi/180))
	# 	d3.append(d1[i]/d2[i])


	# ratio1 = from_data['E_Ratio']
	# dataset2 = from_data['Rx']
	# dataset3 = from_data['Ry']
	# ratio1 = from_data['R_Ratio']

	# print(dataset1.values)

	# dataset1.astype(float)
	# dataset2.astype(float)
	# dataset3.astype(float)
	# ratio1.astype(float)
	#dataset4 = from_data['Ratio']

				# dataset4 = from_data['D_Ratio']
				# dataset5 = from_data['D_max_Angle']
				# dataset1 = dataset1 * ( 0.0340119 / OmegaG )
				# dataset1 *= 1.33
				# dataset2 = dataset2 * Ag
				# dataset3 = dataset3 * Ag
				# value = 0;
				# timeline = []
				# for i in range(0,13):
				# 	value +=  1000.0 * 5000 * 0.0340119 / OmegaG
				# 	timeline.append(value)
				# for i in range(13,26):
				# 	value +=  1000.0 * 5000 * 0.0340119 / OmegaG
				# 	timeline.append(value)
				# for i in range(26,len(dataset1)):
				# 	value +=  1000.0 * 5000 * 0.0340119 / OmegaG 
				# 	timeline.append(value)
				# print dataset1
				# print timeline

	datafile2 = 'runObservables/hydro.dat'
	# datafile3 = 'ode_1000_Rx_Ry.dat'

	cols2 = ["Time","Rx","Ry"]

	from_data2 = pd.read_csv(datafile2,skiprows=1,sep=',',names=cols2)
	data1 = from_data2['Time']
	#data1 *= 1000.0
	data2 = from_data2['Rx']
	data3 = from_data2['Ry']

	data1.astype(float)
	data2.astype(float)
	data3.astype(float)

	# ratio1 = dataset2 / dataset3
	# ratio1 = from_data['D_max/D_min']
	# ratio1[17:] = 1/ratio1[17:]

	ratio2 = data2 / data3
	
	fig = plt.figure()
	fig.set_size_inches(10.0,12.0)

	ax1 = fig.add_subplot(311)
	ax1.plot(dataset1,d1,'.',color='r',label=r'$R_x$ $GPE$')
	ax1.plot(data1,data2,color='r',label=r'$R_x$ $Hydro$')


	name = r'$Radii$'
	name += " "
	name += r'$in$'
	name += " "
	name += r'$\mu m$'
	plt.ylabel(name)
	plt.xlabel('Time in s')
	plt.legend(loc='upper left')

	ax2 = fig.add_subplot(312)
	ax2.plot(dataset1,d2,'.',color='b',label=r'$R_y$ $GPE$')
	ax2.plot(data1,data3,color='b',label=r'$R_y$ $Hydro$')

	plt.ylabel(name)
	plt.xlabel('Time in s')
	plt.legend(loc='upper left')

	# ax1 = fig.add_subplot(312)	


	# plt.ylabel(name)
	# plt.xlabel('Time in ms')
	# plt.legend(loc='upper left',title='1 Vortex')

	ax3 = fig.add_subplot(313)

	name2 = r'$Aspect$ $Ratio$ $R_x$ $/$ $R_y$'
	ax3.plot(dataset1,d3,'.',color='black',label=r'$R_x$$/$$R_y$ $GPE$')
	ax3.plot(data1,ratio2,color='g',label=r'$R_x$$/$$R_y$ $Hydro$')
	plt.ylabel(name2)
	plt.xlabel('Time in ms')
	plt.legend(loc='upper right')


	plt.tight_layout()	
	# fig.text(.001,.001,txt)
	plt.savefig('Rx_Ry_Nv.png',dpi=150)

if __name__ == '__main__':
   	main()
