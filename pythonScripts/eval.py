
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

def draw_plot(bins,data,plot_title,borders):
	edata = np.exp(data)
	ebins = np.exp(bins)
	title(plot_title)
	xlabel("log(k)")
	ylabel("log(n)")
	xlim(borders[0],borders[1]*10.0)
	ylim(borders[2],borders[3]*10.0)
	# plot(range(len(data)),data,alpha=0.8,color='red')
	plot(ebins,edata,alpha=0.8,marker='o',linestyle=' ',color='red')


def draw_segments(bins,segments):
	ax = plt.gca()
	for segment in segments:
		deltaX = (bins[segment[2]] - bins[segment[0]])
		deltaY = (segment[3] - segment[1])
		# print(segment[3],segment[2])
		steigung = deltaY / deltaX
		if abs(segment[2] - segment[0]) > 2:
			if abs(steigung) < 15:
				if abs(steigung) > 1:
					line = Line2D((np.exp(bins[segment[0]]),np.exp(bins[segment[2]])),(np.exp(segment[1]),np.exp(segment[3])))
					ax.add_line(line)
					ax.annotate(str(steigung)[:5],rotation=30,xy=(np.exp(bins[segment[0]] + deltaX/2),np.exp(segment[1] + deltaY/2)),size='8',xytext=(5,30), textcoords='offset points',)
					# +"$\pm$"+str(abs(steigung*0.1))[:5]

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

	h5fileName = ["expdata.h5","rotdata.h5","trapdata.h5"]

	tmp1 = "10000_400x400_0.145_"
	tmp2 = ".0"
	directory = [tmp1+str(x)+tmp2 for x in range(14,22)]
	print(directory)
	for d in directory:
		for f in h5fileName:
			eval(d,f)

def eval(directory,h5fileName):

	h5file = h5py.File(directory+'/'+h5fileName,'r')
	snapshots = h5file.attrs["SnapshotTimes"]

	if not os.path.exists(directory+"/runRegFit"):
		os.makedirs(directory+"/runRegFit")

	borders = []
	borders.append([])
	borders.append([])
	borders.append([])
	borders.append([])
	for times in snapshots:
		kSetName = str(times) + "/Observables/KVector"
		kvector = np.array(h5file[kSetName])
		nSetname = str(times) + "/Observables/OccupationNumber"
		nvector = np.array(h5file[nSetname])
		kvector = np.delete(kvector,0)
		nvector = np.delete(nvector,0)
		borders[0].append(kvector.min())
		borders[1].append(kvector.max())
		borders[2].append(nvector.min())
		borders[3].append(nvector.max())
	borders[0] = np.array(borders[0]).min()
	borders[1] = np.array(borders[1]).max()
	borders[2] = np.array(borders[2]).min()
	borders[3] = np.array(borders[3]).max()


	for times in snapshots:
		print("Fitting " + str(times))
		kSetName = str(times) + "/Observables/KVector"
		kvector = h5file[kSetName]
		nSetname = str(times) + "/Observables/OccupationNumber"
		nvector = h5file[nSetname]
		obs = h5file[str(times) + "/Observables"]
		meta = obs.attrs["Meta"]
		the_time = meta[0] * 1000 #in ms
	
		kLog = np.log(kvector)
		kLog = np.delete(kLog,0)
		nLog = np.log(nvector)
		nLog = np.delete(nLog,0)
	
		total_bins = 200
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
	
		max_error = 0.1
		#sliding window with regression
		figure()
		yscale('log')
		xscale('log')
		# subplot()
		# segments = segment.slidingwindowsegment(data, fit.regression, fit.sumsquared_error, max_error)

		segments = segment.bottomupsegment(data, fit.regression, fit.sumsquared_error, max_error)
		
		# segments = segment.topdownsegment(data, fit.regression, fit.sumsquared_error, max_error)
		
		# segments = segment.slidingwindowsegment(data, fit.interpolate, fit.sumsquared_error, max_error)
		
		# segments = segment.bottomupsegment(data, fit.interpolate, fit.sumsquared_error, max_error)
		
		# segments = segment.topdownsegment(data, fit.interpolate, fit.sumsquared_error, max_error)



		data = np.array(data)
		# data = np.exp(data)
		# bins = np.exp(bins)
		# segments = np.exp(segments)
		draw_plot(bins,data,"Regression Fit of the Spectrum at " + str(the_time) + " ms",borders)
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
		savefig(str(directory)+'/runRegFit/'+str(h5fileName)+'_'+str(times)+'_regfit.png',dpi=150)
		close()


if __name__ == '__main__':
	main()
