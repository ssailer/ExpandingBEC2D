#include <EXP2D_plotter.hpp>

#include <typeinfo>

#define COLOURBAR_MAX_VALUE 10000
#define IMAGE_SIZE 2500

Plotter::Plotter(Eval &e, Options &o) : eval(e), opt(o){
	// eval = e;
	// opt = o;
	const int n = eval.data.meta.grid[0];
	const int m = eval.data.meta.grid[1];
	int k;

	density = mglData(n,m);
	phase = mglData(n,m);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			k = i + n * j;
			density.a[k] = abs2(eval.data.wavefunction[0](i,j));
			if(density.a[k] != 0.0){
				phase.a[k] = arg(eval.data.wavefunction[0](i,j));
			} else {
				phase.a[k] = 0.0;
			}
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
	densityMap();
	contour();
	vortices();

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
			if(eval.totalResult.number(r) != 0.0){
				kval.push_back(eval.totalResult.k(r));
				numberval.push_back(eval.totalResult.number(r));
			}
        }
	}

	auto k_minmax = std::minmax_element(kval.begin(),kval.end());
	auto number_minmax = std::minmax_element(numberval.begin(),numberval.end());

	cout << "k_min: " << *k_minmax.first << " k_max: " << *k_minmax.second << endl;
	cout << "number_min: " << *number_minmax.first << " number_max: " << *number_minmax.second << endl;



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
	gr.Title(to_string(eval.data.meta.time).c_str());
	gr.SetRange('x',*k_minmax.first,*k_minmax.second);
	// gr.SetRange('x',*number_minmax.first,*number_minmax.second);
	gr.SetRange('y',*number_minmax.first,*number_minmax.second);
	gr.SetCoor(11); // log-log-coordinates

	// gr.SubPlot(2,1,0);
	// gr.Axis();
	// gr.Plot(k,number);
	// gr.SubPlot(2,1,1);

	gr.Axis();

	// gr.SetFunc("lg(x)","lg(y)");
	gr.FPlot("x^(-2)","k");
	gr.AddLegend("k^(-2)","k");
	gr.FPlot("x^(-4.66)","r");
	gr.AddLegend("k^(-4.66)","r");
	gr.FPlot("x^(-4)","b");
	gr.AddLegend("k^(-4)","b");

	// gr.Stem(healing_length);
	gr.Plot(k,number," .");
	gr.Legend();

	string name = dirname + "/Spectrum_" + stepsString + ".png";

	gr.WritePNG(name.c_str(),"Spectrum",false);
}

void Plotter::contour(){

	int size = eval.contour[0].size();

	mglData v_x(size);
	mglData v_y(size);

	// float tmp_x[size];
	// float tmp_y[size];

	int l = 0;
	for(std::unordered_set<Coordinate<int32_t>,Hash>::const_iterator it = eval.contour[0].begin(); it != eval.contour[0].end(); ++it){
		v_x.a[l] = it->x();
		v_y.a[l] = it->y();
		// int tmp1 = it->x();
		// int tmp2 = it->y();
		// cout << "tmp1 " << tmp1 << " tmp2 " << tmp2 << endl;
		// v_x.a[l] = tmp1;
		// v_y.a[l] = tmp2;
		
		// cout << "v_x " << v_x.a[l] << " " << "v_y " << v_y.a[l] << " type " << typeid(tmp_x[l]).name() << endl;
		l++;
	}

	// mglData v_x(size,tmp_x);
	// mglData v_y(size,tmp_y);

	mglGraph gr;

	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(3.0);
	gr.SetQuality(3);
	// gr.Title(title.c_str());

	// gr.SetRange('x',-opt.min_x,opt.min_x);
	// gr.SetRange('y',-opt.min_y,opt.min_y);
	gr.SetRange('z',density);
	gr.SetRange('c',density);
	gr.SetRange('x',0,eval.data.meta.grid[0]);
	gr.SetRange('y',0,eval.data.meta.grid[1]);

	gr.Axis();
	gr.Colorbar();
	gr.Dens(density);
	gr.Plot(v_x,v_y," .w");

	string name = dirname + "/Contour_" + stepsString + ".png";

	gr.WritePNG(name.c_str(),"Contour",false);

}

void Plotter::densityMap(){


	int n = eval.data.meta.grid[0];
	int m = eval.data.meta.grid[1];

	// mglComplex data(n,m);
	mglData density(n,m);

	int i,j,k;

	// data.Create(n,m);

	// complex<double> data1;

	for(i=0;i<n;i++) for(j=0;j<m;j++)
	{	
		k = i+n*j;
		density.a[k] = eval.densityLocationMap[0](i,j);
		// if(density.a[k] > 1){
		// 	cout << "PLOTTER" << density.a[k] << " " << i << " " << j << endl;
		// }

		// data.a[k] = abs2(g(0,i,j,0));
	}

	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	string filename = dirname + "/Density_" + stepsString + ".png";

	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(3.0);
	gr.SetQuality(3);
	// gr.Title(title.c_str());
	// gr.Alpha(true);

	gr.SetRange('x',0,eval.data.meta.grid[0]);
	gr.SetRange('y',0,eval.data.meta.grid[1]);
	gr.SetRange('z',density);
	gr.SetRange('c',density);

	// gr.SubPlot(2,2,1);

	// gr.Light(true);
	// gr.Rotate(40,40);
	// gr.Box();
	// gr.Axis();

	// gr.Surf(density);

	// gr.SubPlot(2,2,3);
	gr.Axis();
	gr.Colorbar("_");
	gr.Dens(density);

	gr.WritePNG(filename.c_str(),"Density",false);
}

void Plotter::vortices(){

	int size = eval.pres[0].vlist.size();

	mglData v_x(size);
	mglData v_y(size);

	int l = 0;
	for(list<VortexData>::const_iterator it = eval.pres[0].vlist.begin(); it != eval.pres[0].vlist.end(); ++it){
		v_x.a[l] = it->x.x();
		v_y.a[l] = it->x.y();
		l++;
	}

	mglGraph gr;

	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(3.0);
	gr.SetQuality(3);

	gr.SetRange('x',0,eval.data.meta.grid[0]);
	gr.SetRange('y',0,eval.data.meta.grid[1]);

	gr.Axis();
	gr.Colorbar();
	gr.Dens(phase);
	// gr.Plot(v_x,v_y," #xw");
	gr.Plot(v_x,v_y," .w");

	string name = dirname + "/Vortices_" + stepsString + ".png";

	gr.WritePNG(name.c_str(),"Vortices",false);
}