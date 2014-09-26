
#!/usr/bin/env python
import sys
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():



	# datafile = 'linma1_run07_4x12samples_no_vortices/Combined_Observables.dat'
	datafile = sys.argv[1] + "/CombinedRunObservables/Expanding_Ratios.csv"

	# cols = ["Timestep","X_max","Y_max","D_max","D_min","D_Ratio","D_max_Angle","D_min_Angle","Ratio","RatioAngle","N","V","Density","E_kin"]

	from_data = pd.read_csv(datafile)

	dataset1 = from_data['Timestep']
	dataset2 = from_data['1']
	dataset3 = from_data['30']
	dataset4 = from_data['60']
	dataset5 = from_data['88']

	fig = plt.figure()
	ax1 = fig.add_subplot(411)
	ratioPlot2 = ax1.plot(dataset1,dataset2,'ro')
	plt.ylabel('Angle 0')
	plt.xlabel('Timestep')

	ax2 = fig.add_subplot(412)
	ratioPlot3 = ax2.plot(dataset1,dataset3,'ro')
	plt.ylabel('Angle 30')
	plt.xlabel('Timestep')

	ax3 = fig.add_subplot(413)
	ratioPlot4 = ax3.plot(dataset1,dataset4,'ro')
	plt.ylabel('Angle 60')
	plt.xlabel('Timestep')

	ax4 = fig.add_subplot(414)
	ratioPlot5 = ax4.plot(dataset1,dataset5,'ro')
	plt.ylabel('Angle 88')
	plt.xlabel('Timestep')
	
	plt.savefig(sys.argv[1] + "/Combined_Plots.png")
	plt.show()

	# plt.plot(dataset1,dataset2)
	# plt.show()

if __name__ == '__main__':
   	main()

