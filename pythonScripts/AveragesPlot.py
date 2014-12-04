
#!/usr/bin/env python
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():
	Ag = 150.0 / 2048.0
	m = 87 * 1.66 * 1.0e-27;
	hbar = 1.054 * 1.0e-22;	
	OmegaG = hbar / ( m * Ag * Ag);
	datafile = 'stuff.dat'

	cols = ["Timestep","Rx","Ry"]

	from_data = pd.read_csv(datafile,header=0,names=cols)

	dataset1 = from_data['Timestep']

	dataset2 = from_data['Rx']
	dataset3 = from_data['Ry']
	# dataset4 = from_data['D_Ratio']
	# dataset5 = from_data['D_max_Angle']
	dataset1 = dataset1 * ( 0.0340119 / OmegaG )
	dataset1 *= 1000.0
	dataset2 = dataset2 * Ag
	dataset3 = dataset3 * Ag
	value = 0;
	timeline = []
	for i in range(0,13):
		value +=  1000.0 * 5000 * 0.0340119 / OmegaG
		timeline.append(value)
	for i in range(13,26):
		value +=  1000.0 * 5000 * 0.0340119 / OmegaG / 2.0
		timeline.append(value)
	for i in range(26,len(dataset1)):
		value +=  1000.0 * 5000 * 0.0340119 / OmegaG / 10.0
		timeline.append(value)
	print dataset1
	print timeline

	datafile2 = 'ode_30_Rx_Ry.dat'
	# datafile3 = 'ode_1000_Rx_Ry.dat'

	from_data2 = pd.read_csv(datafile2,header=1,names=cols)
	data1 = from_data2['Timestep']
	data1 *= 1000.0
	data2 = from_data2['Rx']
	data3 = from_data2['Ry']

	ratio1 = dataset2 / dataset3
	ratio2 = data2 / data3


	# from_data3 = pd.read_csv(datafile3,header=1,names=cols)
	# dataa1 = from_data3['Timestep']
	# dataa2 = from_data3['Rx']
	# dataa3 = from_data3['Ry']

	# for i in range(0,len(dataset3)):
	# 	if dataset2[i] < 1:
	# 		dataset2[i] = 1/dataset2[i]

	# hdf5_file_name = 'test06/runData/Expanding_RunEval.h5'
	# group_name = '/800/Observables'
	
	# file = h5py.File(hdf5_file_name, 'r')
	# group = file[group_name]
	
	# dataset1 = group['KVector'][:]
	
	# dataset2 = group['OccupationNumber'][:]
	
	# # for index in range(dataset1)
	# # 	array1[index] = dataset1[index]
	# # 	array2[index] = dataset2[index]
	
	# # array1[0] = dataset1[0]
	
	# file.close()
	# # dataset1.transpose()
	# # dataset2.transpose()
	
	# # stack = np.vstack((dataset1,dataset2))
	# # data1 = np.array(dataset1)
	# # data2 = np.array(dataset2)
	# data = np.column_stack((dataset1,dataset2))
	
	
	fig = plt.figure()

	ax1 = fig.add_subplot(311)
	# ax1.plot(dataset1,dataset2,'ro',color='g',label='Rx')
	# ax1.plot(dataset1,dataset3,'ro',color='b',label='Ry')
	ax1.plot(timeline,dataset2,'.',color='b',label='Rx GPE')
	ax1.plot(data1,data2,color='r',label='Rx Hydro')

	plt.ylabel('Radii in micrometer')
	plt.xlabel('Time in ms')
	plt.legend(loc='upper left')

	ax1 = fig.add_subplot(312)	

	ax1.plot(dataset1,dataset3,'.',color='b',label='Ry GPE')
	ax1.plot(data1,data3,color='r',label='Ry Hydro')
	plt.ylabel('Radii in micrometer')
	plt.xlabel('Time in ms')
	plt.legend(loc='upper left')

	ax1 = fig.add_subplot(313)

	ax1.plot(dataset1,ratio1,'.',color='b',label='Rx/Ry GPE')
	ax1.plot(data1,ratio2,color='r',label='R Hydro')
	plt.ylabel('Aspect Ratio R_x / R_y in micrometer')
	plt.xlabel('Time in ms')
	plt.legend(loc='upper right')	

	# ratioPlot2 = ax1.plot(dataset1,dataset3,'ro')
	# ax1.set_xlim([0.01,4])
	# ax1.set_ylim([0.0001,10000000])
	# ax1.set_xscale('log')
	# ax1.set_yscale('log')
	# plt.ylabel('Rx')
	# plt.xlabel('Timestep')
	# plt.legend(loc='upper left')

	# ax2 = fig.add_subplot(212)
	# anglePlot = 
	# plt.ylabel('Ry')
	# plt.xlabel('Timestep')

	# ax3 = fig.add_subplot(223)
	# anglePlot = ax3.plot(dataset1,dataset4,'ro')
	# plt.ylabel('D_Ratio')
	# plt.xlabel('Timestep')

	# ax4 = fig.add_subplot(224)
	# anglePlot = ax4.plot(dataset1,dataset5,'ro')
	# plt.ylabel('D_max_Angle')
	# plt.xlabel('Timestep')


	# envelope_plot(dataset1,dataset2,winsize=20,ax=ax2)

	# x = np.arange(0.001,10,0.01)

	# g1 = x**(-2)
	# g2 = x**(-4.66)
	# g3 = x**(-5)

	# ax2.plot(x,g1)
	# ax2.plot(x,g2)
	# ax2.plot(x,g3)

	# ax2.set_xlim([0.01,4])
	# ax2.set_ylim([0.0001,10000000])
	# ax2.set_xscale('log')
	# ax2.set_yscale('log')
	# plt.ylabel('Occupation Number')
	# plt.xlabel('radial k-Vector')
	
	
	
	# # scatter = ax2.plot(dataset1,dataset2)
	# # ax2.set_xscale('log')
	# # ax2.set_yscale('log')
	# # plt.ylabel('Occupation Number')
	# # plt.xlabel('radial k-Vector')
	# 
	
	plt.savefig('Rx_Ry_Nv.pdf')
	plt.show()

	# plt.plot(dataset1,dataset2)
	# plt.show()

def envelope_plot(x, y, winsize, ax=None, fill='gray', color='blue'):
    if ax is None:
        ax = plt.gca()
    # Coarsely chunk the data, discarding the last window if it's not evenly
    # divisible. (Fast and memory-efficient)
    numwin = x.size // winsize
    ywin = y[:winsize * numwin].reshape(-1, winsize)
    xwin = x[:winsize * numwin].reshape(-1, winsize)
    # Find the min, max, and mean within each window 
    ymin = ywin.min(axis=1)
    ymax = ywin.max(axis=1)
    ymean = ywin.mean(axis=1)
    xmean = xwin.mean(axis=1)

    fill_artist = ax.fill_between(xmean, ymin, ymax, color=fill, 
                                  edgecolor='none', alpha=0.5)
    line, = ax.plot(xmean, ymean, color=color, linestyle='-')
    return fill_artist, line

if __name__ == '__main__':
   	main()

