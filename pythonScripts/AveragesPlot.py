
#!/usr/bin/env python
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():
	Ag = 300.0 / 2048.0
	m = 87.0 * 1.66e-27;
	hbar = 1.054e-22;	
	OmegaG = hbar / ( m * Ag * Ag);
	datafile = 'stuff2.dat'

	cols = ["Timestep","Time","X_max","Y_max","D_max","D_min","Rx","Ry","D_max/D_min","D_maxAngle","D_minAngle","Ratio","RatioAngle","N","V","N/V","E_kin","N_0"]

	from_data = pd.read_csv(datafile,header=0,sep=',',names=cols)

	dataset1 = from_data['Time']

	dataset2 = from_data['Rx']
	dataset3 = from_data['Ry']
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

	datafile2 = 'ode_Rx_Ry.dat'
	# datafile3 = 'ode_1000_Rx_Ry.dat'

	cols2 = ["Time","Rx", "Ry"]

	from_data2 = pd.read_csv(datafile2,header=1,names=cols2)
	data1 = from_data2['Time']
	data1 *= 1000.0
	data2 = from_data2['Rx']
	data3 = from_data2['Ry']

	ratio1 = dataset2 / dataset3
	# ratio1 = from_data['D_max/D_min']
	# ratio1[17:] = 1/ratio1[17:]

	ratio2 = data2 / data3
	
	fig = plt.figure()
	fig.set_size_inches(10.0,12.0)

	ax1 = fig.add_subplot(311)

	ax1.plot(dataset1,dataset2,'.',color='b',label='Rx GPE')
	ax1.plot(data1,data2,color='r',label='Rx Hydro')

	name = "Radii in "
	name += r'$\mu m$'
	plt.ylabel(name)
	plt.xlabel('Time in ms')
	plt.legend(loc='upper left',title='1 Vortex')

	ax1 = fig.add_subplot(312)	

	ax1.plot(dataset1,dataset3,'.',color='b',label='Ry GPE')
	ax1.plot(data1,data3,color='r',label='Ry Hydro')
	plt.ylabel(name)
	plt.xlabel('Time in ms')
	plt.legend(loc='upper left',title='1 Vortex')

	ax1 = fig.add_subplot(313)

	name2 = "Aspect Ratio R_x / R_y"
	ax1.plot(dataset1,ratio1,'.',color='b',label='Rx/Ry GPE')
	ax1.plot(data1,ratio2,color='r',label='Rx/Ry Hydro')
	plt.ylabel(name2)
	plt.xlabel('Time in ms')
	plt.legend(loc='upper right',title='1 Vortex')

	# txt = ''' These graphs describe a situation with 30 Vortices in the initial setup.'''


	plt.tight_layout()	
	# fig.text(.001,.001,txt)
	plt.savefig('Rx_Ry_Nv.pdf',dpi=600)
	plt.show()
	
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

