#ifndef SAVING_H__
#define SAVING_H__

#include <fstream>

using namespace std;
 

int openDataFiles_obdm(int iterate, double time) //This function can be used for sequential file generation, which is useful when producing movies 
{	
	ofstream fileobdm;  
        char ding[1024]; 
        if(iterate<10)
        snprintf(ding,1024,"OBDM0000%d,%f.dat",iterate,time);
        if(iterate>9 && iterate<100)
        snprintf(ding,1024,"OBDM000%d,%f.dat",iterate,time);
        if(iterate>99 && iterate<1000)
        snprintf(ding,1024,"OBDM00%d,%f.dat",iterate,time);
        if(iterate>999 && iterate<10000)
        snprintf(ding,1024,"OBDM0%d,%f.dat",iterate,time);
        if(iterate>9999)
        snprintf(ding,1024,"OBDM%d,%f.dat",iterate,time);
  
        fileobdm.open(ding);
        fileobdm.setf(ios_base::scientific);
}

int closeDataFiles_obdm(void)
{	
	ofstream fileobdm;  
        fileobdm.close();
        return 0;
}

int save_obdm(double x_axis, double y_axis, double density, double phase) //Used in the save_2D() function
{	
	ofstream fileobdm;  
        fileobdm.precision(16);
        fileobdm.width(24);
        fileobdm<<x_axis;
        fileobdm.precision(16);
        fileobdm.width(24);
        fileobdm<<y_axis;
        fileobdm.precision(16);
        fileobdm.width(24);
        fileobdm<<density;
        fileobdm.precision(16);
        fileobdm.width(24);
        fileobdm<<phase;
 
        fileobdm<<endl;
        return 0;  
}

int save_r(int k, double ratio) //Used in the Aspect_Ratio() function 
{	
	ofstream fileobdm;  
        fileobdm.precision(16);
        fileobdm.width(24);
        fileobdm<<k;  
        fileobdm.precision(16);
        fileobdm.width(24);
        fileobdm<<ratio;
 
        fileobdm<<endl;
        return 0;  
}

int save_track (double x_axis, double y_axis, double timestuff, double density) //Used in track() function
{	
	ofstream fileobdm;  
        fileobdm.precision(16);
        fileobdm.width(24);
        fileobdm << x_axis;
        fileobdm.precision(16);
        fileobdm.width(24);
        fileobdm << y_axis;
        fileobdm.precision(16);
        fileobdm.width(24);
        fileobdm << timestuff; 
        fileobdm.precision(16);
        fileobdm.width(24);
        fileobdm << density; 

        fileobdm << std::endl;
        return 0;  
}

int blank_line(void) //Note that Gnuplot blocks need to be separated by two blank lines
{	
	ofstream fileobdm;  
        fileobdm<<endl;
        return 0;
}

int two_blank_lines(void) //Note that Gnuplot blocks need to be separated by two blank lines
{	
	ofstream fileobdm;  
        fileobdm<<endl;
        fileobdm<<endl;
        return 0;
}

#endif // SAVING_H__
