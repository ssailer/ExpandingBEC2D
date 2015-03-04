#include <iostream>
#include <complex>
#include <string>
#include <mgl2/mgl.h>

#define SUCCESS 0; 
#define ERROR_IN_COMMAND_LINE 1;
#define ERROR_IN_CONFIG_FILE 2;
#define ERROR_UNHANDLED_EXCEPTION 3;

using namespace std;


int main( int argc, char** argv){
try{
	cout << "Plotting Spectrum data file: " << argv[1] << endl;

	string name = argv[1];

	mglData dat(argv[1]);

	mglGraph gr;

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(name.c_str());
	gr.SetRange('x',0.01,4);
	gr.SetRange('y',0.1,10000000);
	gr.SetCoor(11); // log-log-coordinates

	// gr.SubPlot(2,1,0);
	// gr.Axis();
	// gr.Plot(k,number);
	// gr.SubPlot(2,1,1);

	gr.Axis();

	// gr.SetFunc("lg(x)","lg(y)");
	gr.FPlot("x^(-2)");
	gr.FPlot("x^(-4.666)");
	// gr.FPlot("x^(-5)");

	gr.Plot(dat.SubData(1),dat.SubData(2)," .");

	name = name + ".png";

	gr.WritePNG(name.c_str(),"Spectrum",false);
}
catch(std::exception& e) 
{ 
	std::cerr << "Unhandled Exception reached the top of main: " 
			<< e.what() << ", application will now exit" << std::endl; 
	return ERROR_UNHANDLED_EXCEPTION; 
}  
return SUCCESS;  
}