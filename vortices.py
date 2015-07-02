
#!/usr/bin/env python
import h5py
import pandas as pd
import numpy as np
from os import getcwd
from os.path import join
from glob import glob
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():
	cols = ["Timestep","Time","X_max","Y_max","Vortexnumber","Alpha","Rx","Ry","R_Ratio","E_Major","E_Minor","E_Major_Angle","E_Minor_Angle","E_Ratio","N","V","N/V","E_kin","N_0"]

	fig = plt.figure()
	fig.set_size_inches(10.0,12.0)



	filenames = glob(join(getcwd(), '*', 'runObservables'))

	filenames = [s + '/ROT_Observables.dat' for s in filenames]
	
	rot = []
	before = "145_"
	end = "/runObservables"
	for s in filenames:
		name = s[s.find(before) + len(before):s.find(end)]
		rot.append(name)

	ratio = []
	for r in rot:
		ratio.append(float(r)/26.0)
	ax1 = fig.add_subplot(111)
	length = len(filenames)
	# from_data = [None] * length
	
	for i in range(0,length):
		dataset = []
		from_data = pd.read_csv(filenames[i],header=0,sep=',',names=cols)
		dataset1 = from_data['Timestep']
		dataset2 = from_data['Vortexnumber']
		dataset3 = from_data['N']
		plot1 = []
		plot2 = []
		average_length = 20
		avl2 = int(average_length / 2)
		for j in range(avl2,len(dataset1)-avl2):
			av = 0
			for k in range(-avl2,avl2):
				av += ( dataset1[j + k] )
			av /= average_length
			plot1.append(av)
			av = 0
			for k in range(-avl2,avl2):
				av += ( dataset2[j + k] )
			av /= average_length
			plot2.append(av)

		ax1.plot(plot1,plot2,label=str(dataset3[0]) + " " + rot[i] + " " + str(ratio[i]))

	plt.ylabel('Vortexnumber')
	plt.xlabel('Time in ms')
	plt.legend(loc='upper right')

	plt.tight_layout()	
	# fig.text(.001,.001,txt)
	plt.savefig('Vortex.pdf',dpi=600)
	plt.show()

if __name__ == '__main__':
   	main()

