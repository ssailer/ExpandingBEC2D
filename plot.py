
#!/usr/bin/env python
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():
	hdf5_file_name = 'test06/runData/Expanding_RunEval.h5'
	group_name = '/800/Observables'
	
	file = h5py.File(hdf5_file_name, 'r')
	group = file[group_name]
	
	dataset1 = group['KVector'][:]
	
	dataset2 = group['OccupationNumber'][:]
	
	# for index in range(dataset1)
	# 	array1[index] = dataset1[index]
	# 	array2[index] = dataset2[index]
	
	# array1[0] = dataset1[0]
	
	file.close()
	# dataset1.transpose()
	# dataset2.transpose()
	
	# stack = np.vstack((dataset1,dataset2))
	# data1 = np.array(dataset1)
	# data2 = np.array(dataset2)
	data = np.column_stack((dataset1,dataset2))
	
	
	fig = plt.figure()
	ax1 = fig.add_subplot(2,1,1)
	histogram = ax1.plot(dataset1,dataset2,'ro')
	ax1.set_xlim([0.01,4])
	ax1.set_ylim([0.0001,10000000])
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	plt.ylabel('Occupation Number')
	plt.xlabel('radial k-Vector')

	ax2 = fig.add_subplot(2,1,2)

	envelope_plot(dataset1,dataset2,winsize=20,ax=ax2)

	x = np.arange(0.001,10,0.01)

	g1 = x**(-2)
	g2 = x**(-4.66)
	g3 = x**(-5)

	ax2.plot(x,g1)
	ax2.plot(x,g2)
	ax2.plot(x,g3)

	ax2.set_xlim([0.01,4])
	ax2.set_ylim([0.0001,10000000])
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	plt.ylabel('Occupation Number')
	plt.xlabel('radial k-Vector')
	
	
	
	# scatter = ax2.plot(dataset1,dataset2)
	# ax2.set_xscale('log')
	# ax2.set_yscale('log')
	# plt.ylabel('Occupation Number')
	# plt.xlabel('radial k-Vector')
	
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

