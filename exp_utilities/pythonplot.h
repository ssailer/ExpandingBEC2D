#ifndef PYTHONPL0T_H__ //Python Plot
#define PYTHONPL0T_H__
	


void pyplot()
{

for(int i=0;i<snapshot_times.size(); i++)
  {

	ofstream fs;
	fs.open((dir+string("data.py")).c_str(), ios_base::trunc | ios_base::out);
	  

    fs << "#!/usr/bin/python" << endl;
	fs << "# -*- coding: utf-8 -*-" << endl;

	fs << "from matplotlib import use" << endl;
	fs << "use('Agg')" << endl;
	fs << "from matplotlib import rc" << endl;
    fs << "import scipy, pylab" << endl;
    fs << "import scipy.optimize" << endl;
	fs << "import sys" << endl;
    fs << "import matplotlib.colors " << endl;
	fs << "import pylab as p "<< endl;
	fs << "import numpy as np" << endl;
	fs << "import matplotlib.pyplot as plt" << endl;
    fs << "from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, MultipleLocator" << endl;
    fs << "from matplotlib.ticker import FixedFormatter" << endl;
	fs << "from mpl_toolkits.axes_grid1 import make_axes_locatable" << endl;
	fs << "import math as m" << endl;
	fs << "def main():" << endl;
        
	
		fs << "\tpath = '" << dirname << "Spectrum" << "_Path_"  << snapshot_times[i] <<"time"<<".png'" << endl;	
		fs << "\tfig = plt.figure(figsize=(8.3,5.7), dpi=100)" << endl;	
		fs << "\tdata = np.loadtxt('"<< (dir + string("radial_avgs.dat")).c_str() << "', delimiter='\\t', unpack=True, usecols = (0,1,2,3,4,5))" << endl;
		fs << "\ttime=" << snapshot_times[i] << endl;
	    fs << "\tdim ="  << meansFD[i].ares_1.number.size() << endl;
		fs << "\tsnapshot=" << i << endl;

	    fs << "\tstart=" << ((meansFD[i].ares_1.number.size())*i)+1 << endl;
		fs << "\tend=" << (meansFD[i].ares_1.number.size())*(i+1) << endl;

	    fs << "\tstart2=" << ((meansFD[i].ares_1.k.size())*i)+1 << endl;
		fs << "\tend2=" << (meansFD[i].ares_1.k.size())*(i+1) << endl;

	    fs << "\tHealing_k=" << sqrt((opt.N/(opt.grid[0]*opt.grid[1]*opt.grid[2]))*opt.U) << endl;
		

		fs << "\tfig.subplots_adjust(hspace = 0.10, wspace = 0.40, right  = 0.85, left=0.15, bottom = 0.10, top = 0.90)" << endl;

	    fs << "\txfit=np.arange(0.03,0.25,0.01)" << endl;
		fs << "\tyfit=(pow(xfit,(-6))*0.1)  "<< endl;

	  //fs << "\txfit2=np.arange(0.25,2.8,0.01)" << endl;
	  //fs << "\tyfit2=(pow(xfit2,(-2))*20)  "<< endl;

	  // fs << "\txfit3=np.arange(0.03,2.8,0.01)" << endl;
	  // fs << "\tyfit3=(pow(xfit3,(-4))*20)  "<< endl;
       

          //fs << "\txfit4=np.arange(0.03,0.3,0.01)" << endl;
	  //fs << "\tyfit4=(pow(xfit4,(-4.6666666666666666666666666666666666666666))*1)  "<< endl;
       


	//fs << "\tn1=data[2]" << endl;
        //fs << "\tn2=data[1]" << endl;


	    fs << "\tn1=data[1]" << endl;
		fs << "\tn2=data[2]" << endl;
	        
		fs << "\tn_q=data[3]" << endl;
		fs << "\tn_i=data[4]" << endl;
		fs << "\tn_c=data[5]" << endl;

	      





	       
		fs << "\tn4=n2[start:end]" << endl;
		fs << "\tn3=n1[start2:end2]" << endl;

		fs << "\tn_qneu=n_q[start:end]" << endl;
		fs << "\tn_ineu=n_i[start:end]" << endl;
		fs << "\tn_cneu=n_c[start:end]" << endl;	
		
		fs << "\tax = fig.add_subplot(111)" << endl;	
		fs << "\tim=plt.plot(n3, n4,'b.',n3,n_ineu,'y.',n3, n_qneu,'r.',n3,n_cneu,'g.')" << endl;

		  

		fs << "\tplt.setp(im, 'markersize', 3)" << endl; 
		fs << "\tplt.xlim(0.02,3.2)" << endl;


	           
		fs << "\tax.set_title('Time: $"<< snapshot_times[i]<<","<<"Particles:" <<meansFD[i].ares_1.particle_count/opt.N<<"$ ') " << endl;
		fs << "\tax.set_yscale('log')" << endl;
		fs << "\tax.set_xscale('log')" << endl; 
           
	    fs << "\tax.set_xticks([0.1, 0.3, 1, 2, 3, Healing_k])" << endl;
	    fs << "\tax.set_xticklabels(['0.1', '0.3','1', '2', '3','H'])" << endl;
		fs << "\tax.set_xlabel('$k$')" << endl;
		fs << "\tax.set_ylabel('$n(k)$')" << endl;
		fs << "\tax.yaxis.label.set_size(15) " << endl;
		fs << "\tax.xaxis.label.set_size(15) " << endl;
	     
		fs << "\tfig.savefig(path, dpi=100)" << endl;
	    fs << "\tsys.exit(3)" << endl;
    fs << "if __name__=='__main__':" << endl;

	    fs << "\tmain()" << endl;

  

    fs.close();
	
    system((string("python ") + dir + string("Spectrum.py")).c_str());
    }	
}

#endif // PYTHONPL0T_H__