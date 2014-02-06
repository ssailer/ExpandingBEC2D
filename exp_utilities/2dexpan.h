#ifndef SAVING_H__
#define SAVING_H__

#include <fstream>

using namespace std; 
ofstream fileobdm;  

int openDataFiles_obdm(int iterate, int time) //This function can be used for sequential file generation, which is useful when producing movies 
{	

        char ding[1024]; 
        if(iterate<10)
        snprintf(ding,1024,"OBDM0000%d,%d.dat",iterate,time);
        if(iterate>9 && iterate<100)
        snprintf(ding,1024,"OBDM000%d,%d.dat",iterate,time);
        if(iterate>99 && iterate<1000)
        snprintf(ding,1024,"OBDM00%d,%d.dat",iterate,time);
        if(iterate>999 && iterate<10000)
        snprintf(ding,1024,"OBDM0%d,%d.dat",iterate,time);
        if(iterate>9999)
        snprintf(ding,1024,"OBDM%d,%d.dat",iterate,time);
  
        fileobdm.open(ding);
        fileobdm.setf(ios_base::scientific);
}

int closeDataFiles_obdm(void)
{	

        fileobdm.close();
        return 0;
}

int save_obdm(double x_axis,double y_axis, double density, double phase) //Used in the save_2D() function
{	

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

        fileobdm<<endl;
        return 0;
}

int two_blank_lines(void) //Note that Gnuplot blocks need to be separated by two blank lines
{	
 
        fileobdm<<endl;
        fileobdm<<endl;
        return 0;
}

#endif // SAVING_H__
