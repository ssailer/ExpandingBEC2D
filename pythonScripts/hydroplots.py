
#!/usr/bin/env python
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():
	Ag = 100.0 / 2048.0
	m = 87.0 * 1.66e-27;
	hbar = 1.054e-22;	
	OmegaG = hbar / ( m * Ag * Ag);


	datafile = ["5","12","18"]
	cols = ["Time","Rx", "Ry"]

	fig = plt.figure()
	fig.set_size_inches(10.0,12.0)
	ax1 = fig.add_subplot(111)
	# ax2 = fig.add_subplot(212)

	for i in datafile:
		from_data= pd.read_csv(i,header=1,names=cols)
	
		data1 = from_data['Time']
		data1 *= 1000.0
		data2 = from_data['Rx']
		data3 = from_data['Ry']
	
		ratio2 = data2 / data3

		ax1.plot(data1,ratio2,label=i+"Vortices")

		# ax1.plot(data1,data2,color='r',label='Rx Hydro')
		# ax1.plot(data1,data3,color='b',label='Ry Hydro')

		# name = "Radii in "
		# name += r'$\mu m$'
		# plt.ylabel(name)
		# plt.xlabel('Time in ms')
		# plt.legend(loc='upper left',title=i)

	# ax1 = fig.add_subplot(312)	


	# plt.ylabel(name)
	# plt.xlabel('Time in ms')
	# plt.legend(loc='upper left',title='1 Vortex')	

		
		
	name2 = "Aspect Ratio R_x / R_y"
	plt.ylabel(name2)
	plt.xlabel('Time in ms')
	plt.legend(loc='upper right',title='Hydro')

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

