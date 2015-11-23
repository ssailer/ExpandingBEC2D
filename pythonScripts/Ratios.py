
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

	name = r'$Radii$'
	name += " "
	name += r'$in$'
	name += " "
	name += r'$\mu m$'

	name2 = r'$Aspect$ $Ratio$ $R_x$ $/$ $R_y$'

	cols = ["Timestep","Time","X_max","Y_max","Vortexnumber","Alpha","Rx","Ry","R_Ratio","E_Major","E_Minor","E_Major_Angle","E_Minor_Angle","E_Ratio","N","V","N/V","E_kin","N_0"]

	fig = plt.figure()
	fig.set_size_inches(10.0,12.0)
	
	stringone = 'Rx'
	stringtwo = 'Ry'
	stringthree = 'R_Ratio'
	



	datafile1 = '10000_400x400_0.145_14.0/runObservables/EXP_Observables.dat'
	datafile2 = '10000_400x400_0.145_15.0/runObservables/EXP_Observables.dat'
	datafile3 = '10000_400x400_0.145_16.0/runObservables/EXP_Observables.dat'
	datafile4 = '10000_400x400_0.145_17.0/runObservables/EXP_Observables.dat'
	datafile5 = '10000_400x400_0.145_18.0/runObservables/EXP_Observables.dat'
	datafile6 = '10000_400x400_0.145_19.0/runObservables/EXP_Observables.dat'
	datafile7 = '10000_400x400_0.145_20.0/runObservables/EXP_Observables.dat'
	datafile8 = '10000_400x400_0.145_21.0/runObservables/EXP_Observables.dat'
	# datafile9 = ''
	

	from_data = pd.read_csv(datafile1,header=0,sep=',',names=cols)

	dataset11 = from_data['Time']
	dataset12 = from_data[stringone]
	dataset13 = from_data[stringtwo]
	dataset14 = from_data[stringthree]
	ratio1 = dataset12 / dataset13

	from_data = pd.read_csv(datafile2,header=0,sep=',',names=cols)

	dataset21 = from_data['Time']
	dataset22 = from_data[stringone]
	dataset23 = from_data[stringtwo]
	dataset24 = from_data[stringthree]
	ratio2 = dataset22 / dataset23

	from_data = pd.read_csv(datafile3,header=0,sep=',',names=cols)

	dataset31 = from_data['Time']
	dataset32 = from_data[stringone]
	dataset33 = from_data[stringtwo]
	dataset34 = from_data[stringthree]
	ratio3 = dataset32 / dataset33

	from_data = pd.read_csv(datafile4,header=0,sep=',',names=cols)

	dataset41 = from_data['Time']
	dataset42 = from_data[stringone]
	dataset43 = from_data[stringtwo]
	dataset44 = from_data[stringthree]
	ratio4 = dataset42 / dataset43

	from_data = pd.read_csv(datafile5,header=0,sep=',',names=cols)

	dataset51 = from_data['Time']
	dataset52 = from_data[stringone]
	dataset53 = from_data[stringtwo]
	dataset54 = from_data[stringthree]
	ratio5 = dataset52 / dataset53

	from_data = pd.read_csv(datafile6,header=0,sep=',',names=cols)

	dataset61 = from_data['Time']
	dataset62 = from_data[stringone]
	dataset63 = from_data[stringtwo]
	dataset64 = from_data[stringthree]
	ratio6 = dataset62 / dataset63

	from_data = pd.read_csv(datafile7,header=0,sep=',',names=cols)

	dataset71 = from_data['Time']
	dataset72 = from_data[stringone]
	dataset73 = from_data[stringtwo]
	dataset74 = from_data[stringthree]
	ratio7 = dataset72 / dataset73

	from_data = pd.read_csv(datafile8,header=0,sep=',',names=cols)

	dataset81 = from_data['Time']
	dataset82 = from_data[stringone]
	dataset83 = from_data[stringtwo]
	dataset84 = from_data[stringthree]
	ratio8 = dataset82 / dataset83

	# from_data = pd.read_csv(datafile9,header=0,sep=',',names=cols)

	# dataset91 = from_data['Time']
	# dataset92 = from_data['Rx']
	# dataset93 = from_data['Ry']
	# dataset94 = from_data['Ratio']
	# ratio9 = dataset92 / dataset93


	ax1 = fig.add_subplot(311)
	ax1.plot(dataset11,dataset12,color='r',label=r'$14$ $Vortices$')	
	ax1.plot(dataset21,dataset22,color='b',label=r'$15$ $Vortices$')
	ax1.plot(dataset31,dataset32,color='g',label=r'$16$ $Vortices$')	
	ax1.plot(dataset41,dataset42,color='y',label=r'$17$ $Vortices$')
	ax1.plot(dataset51,dataset52,color='c',label=r'$18$ $Vortices$')
	ax1.plot(dataset61,dataset62,color='k',label=r'$19$ $Vortices$')
	ax1.plot(dataset71,dataset72,color='m',label=r'$20$ $Vortices$')
	ax1.plot(dataset81,dataset82,color='k',label=r'$21$ $Vortices$')
	# ax1.plot(dataset91,dataset92,color='r',label=r'$98$ $Vortices$')
	# ax1.plot(data1,data2,color='r',label=r'$R_x$ $Hydro$')
	plt.ylabel(name)
	plt.xlabel('Time in ms')
	plt.legend(loc='upper left',title=r'$R_x$')
	
	ax2 = fig.add_subplot(312)
	ax2.plot(dataset11,dataset13,color='r',label=r'$14$ $Vortices$')
	ax2.plot(dataset21,dataset23,color='b',label=r'$15$ $Vortices$')	
	ax2.plot(dataset31,dataset33,color='g',label=r'$16$ $Vortices$')	
	ax2.plot(dataset41,dataset43,color='y',label=r'$17$ $Vortices$')
	ax2.plot(dataset51,dataset53,color='c',label=r'$18$ $Vortices$')
	ax2.plot(dataset61,dataset63,color='k',label=r'$19$ $Vortices$')
	ax2.plot(dataset71,dataset73,color='m',label=r'$20$ $Vortices$')	
	ax2.plot(dataset81,dataset83,color='k',label=r'$21$ $Vortices$')
	# ax2.plot(dataset91,dataset93,color='r',label=r'$98$ $Vortices$')
	# ax2.plot(data1,data3,color='b',label=r'$R_y$ $Hydro$')
	plt.ylabel(name)
	plt.xlabel('Time in ms')
	plt.legend(loc='upper left',title=r'$R_y$')

	ax3 = fig.add_subplot(313)
	ax3.plot(dataset11,ratio1,color='r',label=r'$14$ $Vortices$')	
	ax3.plot(dataset21,ratio2,color='b',label=r'$15$ $Vortices$')
	ax3.plot(dataset31,ratio3,color='g',label=r'$16$ $Vortices$')	
	ax3.plot(dataset41,ratio4,color='y',label=r'$17$ $Vortices$')
	ax3.plot(dataset51,ratio5,color='c',label=r'$18$ $Vortices$')
	ax3.plot(dataset61,ratio6,color='k',label=r'$19$ $Vortices$')
	ax3.plot(dataset71,ratio7,color='m',label=r'$20$ $Vortices$')	
	ax3.plot(dataset81,ratio8,color='k',label=r'$21$ $Vortices$')
	# ax3.plot(dataset91,ratio9,color='r',label=r'$98$ $Vortices$')
	# ax3.plot(data1,ratio2,color='g',label=r'$R_x$$/$$R_y$ $Hydro$')
	plt.ylabel(name2)
	plt.xlabel('Time in ms')
	plt.legend(loc='upper right',title=r'$R_x$$/$$R_y$')



	plt.tight_layout()	
	# fig.text(.001,.001,txt)
	plt.savefig('GPE_Results.pdf',dpi=600)
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

