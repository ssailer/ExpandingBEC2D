#include <plot_with_mgl.h>

#define COLOURBAR_MAX_VALUE 10000

using namespace std;


void plotspectrum(string name,Observables &ares){

	ofstream plotfile;
    vector<double> kval;
	vector<double> numberval;
	
	plotfile.open((name + ".dat").c_str(), ios::out | ios::trunc);
    for (int r = 0; r < ares.number.size(); r++)             
	{	
		if(ares.k(r) != 0.0){
			plotfile << r <<"\t"<< ares.k(r) <<"\t" << ares.number(r) <<"\t";
			plotfile << endl;

			kval.push_back(ares.k(r));
			numberval.push_back(ares.number(r));
        }
	}
	plotfile << endl << endl;	
	plotfile.close();


	int n = kval.size();//-1; // don't plot the zero mode! (why? because it looks like shit)

	mglData k(n);
	mglData number(n);
	// mglData healing_length(2);
	// healing_length.a[0] = ares.healing_length;
	// healing_length.a[1] = 10000000;

	for(int i = 0; i < n; i++){
		k.a[i] = kval[i];
		number.a[i] = numberval[i];
	}

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

	// gr.Stem(healing_length);
	gr.Plot(k,number," .");

	name = name + ".png";

	gr.WritePNG(name.c_str(),"Spectrum",false);
}

void plotVortexLocationMap(string name,RealGrid &VortexLocationMap){

	int n = VortexLocationMap.width();
	int m = VortexLocationMap.height();

	mglData data(n,m);
	int i,j,k;

	for(i=0;i<n;i++) for(j=0;j<m;j++)
	{
		k = i+n*j;
		data.a[k] = VortexLocationMap(0,i,j,0);		
	}

	mglGraph gr;

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(name.c_str());

	gr.SetRange('x',-2,2);
	gr.SetRange('y',-2,2);
	gr.SetRange('z',data);
	gr.SetRange('c',data);

	// gr.Axis();
	// gr.Colorbar("_");
	// gr.Dens(data);

	gr.Rotate(40,40);
	gr.Box();
	gr.Axis();
	gr.Surf(data);

	gr.WritePNG(name.c_str(),"ExpandingVortexGas2D",false);
}

void plotdatatopng(string filename,ComplexGrid* &g,Options &opt)
{
	

	int n = opt.grid[1];
	int m = opt.grid[2];

	// mglComplex data(n,m);
	mglData density(n,m);
	mglData phase(n,m);

	int i,j,k;

	// data.Create(n,m);

	// complex<double> data1;

	for(i=0;i<n;i++) for(j=0;j<m;j++)
	{	
		k = i+n*j;
		density.a[k] = abs2(g->at(0,i,j,0));
		phase.a[k] = arg(g->at(0,i,j,0));

		// data.a[k] = abs2(g(0,i,j,0));
	}

	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	filename = filename + ".png";

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(filename.c_str());
	// gr.Alpha(true);



	// data.use_abs=false;
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',-opt.min_y,opt.min_y);
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
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',-opt.min_y,opt.min_y);
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

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}

void plotdatatopng(string filename,ComplexGrid &g,Options &opt)
{
	

	int n = opt.grid[1];
	int m = opt.grid[2];

	// mglComplex data(n,m);
	mglData density(n,m);
	mglData phase(n,m);

	int i,j,k;

	// data.Create(n,m);

	// complex<double> data1;

	for(i=0;i<n;i++) for(j=0;j<m;j++)
	{	
		k = i+n*j;
		density.a[k] = abs2(g(0,i,j,0));
		phase.a[k] = arg(g(0,i,j,0));

		// data.a[k] = abs2(g(0,i,j,0));
	}

	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	filename = filename + ".png";

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(filename.c_str());
	// gr.Alpha(true);



	// data.use_abs=false;
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',-opt.min_y,opt.min_y);
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
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',-opt.min_y,opt.min_y);
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

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}

void plotdatatopng(string filename,RealGrid &g,Options &opt){

	int n = opt.grid[1];
	int m = opt.grid[2];

	// mglComplex data(n,m);
	mglData grid(n,m);


	int i,j,k;

	// data.Create(n,m);

	// complex<double> data1;

	for(i=0;i<n;i++) for(j=0;j<m;j++)
	{	
		k = i+n*j;
		grid.a[k] = g(0,i,j,0);
		// data.a[k] = abs2(g(0,i,j,0));
	}

	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	filename = filename + ".png";

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(filename.c_str());
	// gr.Alpha(true);

	// data.use_abs=true;
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',-opt.min_y,opt.min_y);
	gr.SetRange('z',grid);
	gr.SetRange('c',grid);

	gr.Axis();
	gr.Colorbar("_");
	gr.Dens(grid);

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}

void plotdatatopngEigen(Eigen::MatrixXcd& wavefct,Options &opt)
{
	

	int n = opt.grid[1];
	int m = opt.grid[2];

	mglComplex data(n,m);

	int i,j,k;

	// data.Create(n,m);

	// complex<double> data1;

	for(i=0;i<n;i++) for(j=0;j<m;j++)
	{	
		k = i+n*j;
		// data1 = g->at(0,i,j,0);
		data.a[k] = abs2(wavefct(i,j));
	}

	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',-opt.min_y,opt.min_y);
	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(opt.name.c_str());
	// gr.Alpha(true);



	// data.use_abs=false;
	// string filename = "PHASE-" + opt.name + ".png";

	// gr.SetRange('z',data);
	// gr.SetRange('c',data);

	// // gr.SubPlot(1,2,0);

	// // gr.Rotate(40,40);
	// // gr.Box();
	// // gr.Axis();
	// // gr.Surf(data);


	// // gr.SubPlot(2,2,2);
	// gr.Axis();
	// gr.Colorbar("_");
	// gr.Dens(data);
	// gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);


	data.use_abs=true;
	string filename = opt.name + "-Density-Expanding.png";
	gr.SetRange('z',data);
	// gr.SetRange('c',data);
	gr.SetRange('c',data);

	// gr.SubPlot(1,2,1);

	// // gr.Light(true);
	// gr.Rotate(40,40);
	// gr.Box();
	// gr.Axis();

	// gr.Surf(data);

	// gr.SubPlot(2,2,3);
	gr.Axis();
	gr.Colorbar("_");
	gr.Dens(data);

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}

void plotdatatopngEigenExpanding(Eigen::MatrixXcd& mPsi,vector<double> &ranges,Eigen::VectorXd &Xexpanding,Eigen::VectorXd &Yexpanding,Options &opt)
{
	

	int n = opt.grid[1];
	int m = opt.grid[2];



	mglData data(n,m);
	mglData xaxis(n);
	mglData yaxis(m);

	int i,j,k;

	// data.Create(n,m);

	for(i=0;i<n;i++) for(j=0;j<m;j++)
	{	
		k = i+n*j;

		data.a[k] = abs2(mPsi(i,j));		
	}

	for( i = 0; i < n; i++){ xaxis.a[i] = Xexpanding(i); }
	for( j = 0; j < m; j++){ yaxis.a[j] = Yexpanding(j); }



	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	string filename = opt.name + "-Density-Expanding.png";

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(opt.name.c_str());
	gr.SetRange('x',xaxis);
	gr.SetRange('y',yaxis);
	// gr.Alpha(true);


	// data.use_abs=false;

	// gr.SetRange('z',data);
	// gr.SetRange('c',data);

	// gr.SubPlot(2,2,0);

	// gr.Rotate(40,40);
	// gr.Box();
	// gr.Axis();
	// gr.Surf(xaxis,yaxis,data);


	// gr.SubPlot(2,2,2);
	// gr.Axis();
	// gr.Colorbar("_");
	// gr.Dens(xaxis,yaxis,data);


	// data.use_abs=true;
	gr.SetRange('z',data);
	// gr.SetRange('c',data);
	gr.SetRange('c',data);

	// gr.SubPlot(2,2,1);

	// // gr.Light(true);
	// gr.Rotate(40,40);
	// gr.Box();
	// gr.Axis();

	// gr.Surf(xaxis,yaxis,data);

	// gr.SubPlot(2,2,3);
	gr.Axis();
	gr.Colorbar("_");
	gr.Dens(xaxis,yaxis,data);

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}


void plotVector(string filename,vector<double> v,vector<double> w,Options &opt){

	int x = v.size();
	int y = w.size();

	mglData v_data(x);
	mglData w_data(y);

	for(int i=0;i < x;i++){	
		v_data.a[i] = v[i];
	}
	for(int j = 0; j < y; j++){
		w_data.a[j] = w[j];
	}

	mglGraph gr;
		
	filename = filename + ".png";

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(filename.c_str());
	// gr.Alpha(true);

	gr.SubPlot(2,1,0);
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',v_data);
	gr.Axis();
	gr.Plot(v_data);

	gr.SubPlot(2,1,1);
	gr.SetRange('x',-opt.min_y,opt.min_y);
	gr.SetRange('y',w_data);
	gr.Axis();
	gr.Plot(w_data);

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}

void plotVector(string filename,vector<double> v,Options &opt){

	int x = v.size();

	mglData data(x);

	for(int i=0;i < x;i++){	
		data.a[i] = v[i];
	}

	mglGraph gr;
		
	filename = filename + ".png";

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(filename.c_str());

	gr.SetRange('x',0,v.size());
	gr.SetRange('y',data);
	gr.Axis();
	gr.Plot(data);

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}

void plotVector(string filename,ArrayXd v,Options &opt){

	int x = v.size();

	mglData data(x);

	for(int i=0;i < x;i++){	
		data.a[i] = v[i];
	}

	mglGraph gr;
		
	filename = filename + ".png";

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(filename.c_str());

	gr.SetRange('x',0,v.size());
	gr.SetRange('y',data);
	gr.Axis();
	gr.Plot(data);

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}

void plotAngularDensity(string filename,vector<double> phi,vector<double> density,Options &opt){

	int x = phi.size();
	int y = density.size();

	mglData phi_data(x);
	mglData density_data(y);

	for(int i=0;i < x;i++){	
		phi_data.a[i] = phi[i];
	}
	for(int j = 0; j < y; j++){
		density_data.a[j] = density[j];
	}

	mglGraph gr;
		
	filename = filename + ".png";

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(filename.c_str());

	gr.SetRange('x',phi_data);
	gr.SetRange('y',density_data);
	gr.Axis();
	gr.Plot(phi_data,density_data);

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}