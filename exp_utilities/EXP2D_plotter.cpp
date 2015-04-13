#include <EXP2D_plotter.hpp>

#define COLOURBAR_MAX_VALUE 10000
#define IMAGE_SIZE 2500

Plotter::Plotter(Eval &e, Options &o){
	eval = e;
	opt = o;
	const int n = eval.data.meta.grid[0];
	const int m = eval.data.meta.grid[1];
	int k;

	density = mglData(n,m);
	phase = mglData(n,m);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			k = i + n * j;
			density.a[k] = abs2(eval.data.wavefunction[0](i,j));
			phase.a[k] = arg(eval.data.wavefunction[0](i,j));
		}
	}
	stepsString = to_string(eval.data.meta.steps);

	xrange = eval.data.meta.coord[0];
	yrange = eval.data.meta.coord[1];

	dirname = "runPlots";
    struct stat st;
    	if(stat(dirname.c_str(),&st) != 0){
        mkdir(dirname.c_str(),0755);
    }
}

Plotter::~Plotter(){}

void Plotter::plotEval(){
	control();
	spectrum();
}


void Plotter::control(){

	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	string filename = dirname + "/Control_" + stepsString + ".png";

	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(3.0);
	gr.SetQuality(3);
	// gr.Title(title.c_str());
	// gr.Alpha(true);



	// data.use_abs=false;
	gr.SetRange('x',-xrange,xrange);
	gr.SetRange('y',-yrange,yrange);
	gr.SetRange('z',phase);
	gr.SetRange('c',phase);

	gr.SubPlot(2,2,0);

	gr.Rotate(40,40);
	gr.Box();
	gr.Axis();
	gr.Surf(phase);


	gr.SubPlot(2,2,2);
	gr.Axis();
	gr.Colorbar("_");
	gr.Dens(phase);


	// data.use_abs=true;
	gr.SetRange('x',-xrange,xrange);
	gr.SetRange('y',-yrange,yrange);
	gr.SetRange('z',density);
	gr.SetRange('c',density);

	gr.SubPlot(2,2,1);

	// gr.Light(true);
	gr.Rotate(40,40);
	gr.Box();
	gr.Axis();

	gr.Surf(density);

	gr.SubPlot(2,2,3);
	gr.Axis();
	gr.Colorbar("_");
	gr.Dens(density);

	gr.WritePNG(filename.c_str(),"Control",false);

}

void Plotter::spectrum(){
    vector<double> kval;
	vector<double> numberval;
	
    for (int r = 0; r < eval.totalResult.number.size(); r++)             
	{	
		if(eval.totalResult.k(r) != 0.0){
				kval.push_back(eval.totalResult.k(r));
				numberval.push_back(eval.totalResult.number(r));
        }
	}



	int n = kval.size();//-1; // don't plot the zero mode! (why? because it looks like shit)
	mglData k(n);
	mglData number(n);



	for(int i = 0; i < n; i++){
		k.a[i] = kval[i];
		number.a[i] = numberval[i];
	}

	// cout << "copied" << endl;  

	mglGraph gr;

	gr.SetMarkSize(0.7);
	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(3.0);
	gr.SetQuality(3);
	// gr.Title(title.c_str());
	gr.SetRange('x',k);
	gr.SetRange('y',0.5,1.0e10);
	gr.SetCoor(11); // log-log-coordinates

	// gr.SubPlot(2,1,0);
	// gr.Axis();
	// gr.Plot(k,number);
	// gr.SubPlot(2,1,1);

	gr.Axis();

	// gr.SetFunc("lg(x)","lg(y)");
	gr.FPlot("x^(-2)");
	// gr.FPlot("x^(-4.66)");
	gr.FPlot("x^(-5)");

	// gr.Stem(healing_length);
	gr.Plot(k,number," .");

	string name = dirname + "/Spectrum_" + stepsString + ".png";

	gr.WritePNG(name.c_str(),"Spectrum",false);
}