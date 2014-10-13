
#!/usr/bin/env python
import sys
# import h5py
# import pandas as pd
# import csv
import numpy as np
import matplotlib.pyplot as plt
# from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():



	# datafile = 'linma1_run07_4x12samples_no_vortices/Combined_Observables.dat'
	datafile1 = sys.argv[1] + "/CombinedRunObservables/Expanding_Ratios.csv"
	datafile2 = sys.argv[2] + "/CombinedRunObservables/Expanding_Ratios.csv"
	angle = int(sys.argv[3])

	# cols = ["Timestep","X_max","Y_max","D_max","D_min","D_Ratio","D_max_Angle","D_min_Angle","Ratio","RatioAngle","N","V","Density","E_kin"]

	# from_data1 = pd.read_csv(datafile1)
	# from_data2 = pd.read_csv(datafile2)

	dataset01=np.loadtxt(datafile1,dtype=float,delimiter=',',skiprows=1,usecols=(1,))
	dataset02=np.loadtxt(datafile1,dtype=float,delimiter=',',skiprows=1,usecols=(angle+2,))

	dataset11=np.loadtxt(datafile2,dtype=float,delimiter=',',skiprows=1,usecols=(1,))
	dataset12=np.loadtxt(datafile2,dtype=float,delimiter=',',skiprows=1,usecols=(angle+2,))

	# from_data1 = open(datafile1,"rb")
	# from_data2 = open(datafile2,"rb")

	# for row in csv.reader(from_data1)

	# dataset01 = from_data['Timestep']
	# dataset02 = from_data['1']
	# dataset03 = from_data['30']
	# dataset04 = from_data['60']
	# dataset05 = from_data['88']

	# dataset11 = from_data['Timestep']
	# dataset12 = from_data['1']
	# dataset13 = from_data['30']
	# dataset14 = from_data['60']
	# dataset15 = from_data['88']

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ratioPlot2 = ax1.plot(dataset01,dataset02, c='b', marker='o', label='first')
	ratioPlot3 = ax1.plot(dataset11,dataset12, c='r', marker='o', label='second')
	plt.ylabel('Angle 0')
	plt.xlabel('Timestep')
	plt.legend(loc='upper right')

	# ax2 = fig.add_subplot(412)
	# ratioPlot3 = ax2.plot(dataset1,dataset3,'ro')
	# plt.ylabel('Angle 30')
	# plt.xlabel('Timestep')

	# ax3 = fig.add_subplot(413)
	# ratioPlot4 = ax3.plot(dataset1,dataset4,'ro')
	# plt.ylabel('Angle 60')
	# plt.xlabel('Timestep')

	# ax4 = fig.add_subplot(414)
	# ratioPlot5 = ax4.plot(dataset1,dataset5,'ro')
	# plt.ylabel('Angle 88')
	# plt.xlabel('Timestep')
	
	plt.savefig(sys.argv[1] + "/Combined_Plots.png")
	plt.show()

	# plt.plot(dataset1,dataset2)
	# plt.show()

if __name__ == '__main__':
   	main()

