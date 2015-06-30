
#!/usr/bin/env python
import h5py
import pandas as pd
import numpy as np
import matplotlib
import sys
import os
import statsmodels.api as smf
from scipy import interpolate
from math import cos,sin,pi
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from scipy.cluster.vq import kmeans, kmeans2, whiten

from matplotlib.pylab import gca, figure, plot, subplot, title, xlabel, ylabel, xlim, ylim, show, savefig, close, yscale, xscale
from matplotlib.lines import Line2D
import segment
import fit

def draw_plot(bins,data,plot_title):
	edata = np.exp(data)
	ebins = np.exp(bins)
	# plot(range(len(data)),data,alpha=0.8,color='red')
	plot(ebins,edata,alpha=0.8,marker='o',linestyle=' ',color='red')
	title(plot_title)
	xlabel("k in 1/m")
	ylabel("n OccupationNumber")
	xlim(ebins.min(),ebins.max())
	ylim(edata.min(),edata.max())

def draw_segments(bins,segments):
	ax = gca()
	for segment in segments:
		deltaX = (bins[segment[2]] - bins[segment[0]])
		deltaY = (segment[3] - segment[1])
		# print(segment[3],segment[2])
		steigung = deltaY / deltaX
		if abs(segment[2] - segment[0]) > 2:
			if abs(steigung) < 10:
				line = Line2D((np.exp(bins[segment[0]]),np.exp(bins[segment[2]])),(np.exp(segment[1]),np.exp(segment[3])))
				ax.add_line(line) 
				ax.annotate(str(steigung)[:5],xy=(np.exp(bins[segment[0]] + deltaX/2),np.exp(segment[1] + deltaY/2)),xytext=(5,10), textcoords='offset points',)

def window(size):
	return np.ones(size)/float(size)

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

def main():

	h5fileName = sys.argv[1]
	h5file = h5py.File(h5fileName,'r')
	snapshots = h5file.attrs["SnapshotTimes"]
	directory = "runRegFit"

	if not os.path.exists(directory):
		os.makedirs(directory)

	for times in snapshots:
		print("Fitting " + str(times))
		kSetName = str(times) + "/Observables/KVector"
		kvector = h5file[kSetName]
		nSetname = str(times) + "/Observables/OccupationNumber"
		nvector = h5file[nSetname]
	
		kLog = np.log(kvector)
		kLog = np.delete(kLog,0)
		nLog = np.log(nvector)
		nLog = np.delete(nLog,0)
	
		total_bins = 100
		bins = np.linspace(kLog.min(),kLog.max(),total_bins)
		delta = bins[1] - bins[0]
		idx = np.digitize(kLog,bins)
		running_median = [np.median(nLog[idx == k]) for k in range(total_bins)]
	
		data = running_median
		data = np.array(data)
		# data = np.nan_to_num(data)
		nans, x= nan_helper(data)
		data[nans]= np.interp(x(nans), x(~nans), data[~nans])
		data = data.tolist()
	
		max_error = 0.05
		#sliding window with regression
		figure()
		yscale('log')
		xscale('log')
		# subplot()
		# segments = segment.slidingwindowsegment(data, fit.regression, fit.sumsquared_error, max_error)

		# segments = segment.bottomupsegment(data, fit.regression, fit.sumsquared_error, max_error)
		
		# segments = segment.topdownsegment(data, fit.regression, fit.sumsquared_error, max_error)
		
		segments = segment.slidingwindowsegment(data, fit.interpolate, fit.sumsquared_error, max_error)
		
		# segments = segment.bottomupsegment(data, fit.interpolate, fit.sumsquared_error, max_error)
		
		# segments = segment.topdownsegment(data, fit.interpolate, fit.sumsquared_error, max_error)



		data = np.array(data)
		# data = np.exp(data)
		# bins = np.exp(bins)
		# segments = np.exp(segments)
		draw_plot(bins,data,"Sliding window with regression")
		draw_segments(bins,segments)
	
	
		# fig = plt.figure()
		# fig.set_size_inches(10.0,12.0)
	
		# ax1 = fig.add_subplot(211)
		# subplot()
		# set_yscale('log')
		# set_xscale('log')
		# plot(kvector,nvector,'.',color='r',label=r'$R_x$ $GPE$')
	
	
	
		# ax2 = fig.add_subplot(212)
		# ax2.plot(bins-delta/2,running_median,'.')
		# plot(kLog,nLog)
		# # ax2.plot(kLog,interpolate.splev(kLog,tck,der=0))
	
		# plt.tight_layout()  
		# show()
		# # fig.text(.001,.001,txt)
		savefig(str(directory)+'/'+str(sys.argv[1])+'_'+str(times)+'_regfit.png',dpi=150)
		close()


if __name__ == '__main__':
	main()
