#include <EXP2D_plotter.hpp>

#include <typeinfo>

#define COLOURBAR_MAX_VALUE 10000
#define IMAGE_SIZE 2500
#define FONT_SIZE 2.0

Plotter::Plotter(Eval &e, Options &o) : eval(e), opt(o){
	// eval = e;
	// opt = o;

	prepareData();

	stepsString = to_string(eval.data.meta.steps);

	xrange = eval.data.meta.coord[0];
	yrange = eval.data.meta.coord[1];

	dirname = "runPlots";
    struct stat st;
    	if(stat(dirname.c_str(),&st) != 0){
        mkdir(dirname.c_str(),0755);
    }
    double currentTime = eval.data.meta.time;
	std::ostringstream out;
	if(currentTime <= 0.1){
		currentTime *= 1000;
    	out << std::setprecision(5) << currentTime;
    	title = "t = " + out.str() + " ms";
    } else {
    	out << std::setprecision(5) << currentTime;
    	title = "t = " + out.str() + " s";
    }
}

void Plotter::prepareData(){

	const int n = eval.data.meta.grid[0];
	const int m = eval.data.meta.grid[1];
	const int contour_size = eval.contour[0].size();
	const int vortex_size = eval.pres[0].vlist.size();
	int index;

	// Phase and Density
	density = mglData(n,m);
	phase = mglData(n,m);
	densitymap = mglData(n,m);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			index = i + n * j;
			density.a[index] = abs2(eval.data.wavefunction[0](i,j));
			if(density.a[index] != 0.0){
				phase.a[index] = arg(eval.data.wavefunction[0](i,j));
			} else {
				phase.a[index] = 0.0;
			}
			densitymap.a[index] = eval.densityLocationMap[0](i,j);
		}
	}

	// contour vectors
	contour_x = mglData(contour_size); 
	contour_y = mglData(contour_size);

	index = 0;
	for(std::unordered_set<Coordinate<int32_t>,Hash>::const_iterator it = eval.contour[0].begin(); it != eval.contour[0].end(); ++it){
		contour_x.a[index] = it->x();
		contour_y.a[index] = it->y();
		index++;
	}

	// vortex vectors
	vortex_x = mglData(vortex_size);
	vortex_y = mglData(vortex_size);

	index = 0;
	for(list<VortexData>::const_iterator it = eval.pres[0].vlist.begin(); it != eval.pres[0].vlist.end(); ++it){
		vortex_x.a[index] = it->x.x();
		vortex_y.a[index] = it->x.y();
		index++;
	}


	// spectrum
	vector<double> kval;
	vector<double> numberval;
	map<double,double> spectrum;
    pair<map<double,double>::iterator,bool> ret;
		
    for (int r = 0; r < eval.totalResult.number.size(); r++){
		if(eval.totalResult.k(r) != 0.0){
			if(eval.totalResult.number(r) != 0.0){
				ret = spectrum.insert(map<double,double>::value_type(eval.totalResult.k(r),eval.totalResult.number(r)));
				if(ret.second==false){
					cout << "Binning of spectrum failed, double value inserted." << endl;
				}
				// kval.push_back(k_int);
				// numberval.push_back(eval.totalResult.number(r));
			}
        }
	}

	int binsize = 200;
	int c1 = 0;
	// int c2 = 0;
	double ksum = 0;
	double nsum = 0;
	for(map<double,double>::const_iterator it = spectrum.begin(); it != spectrum.end(); ++it){
		c1++;
		ksum += it->first;
		nsum += it->second;
		if(c1 == binsize){
			c1 = 0;
			ksum /= binsize;
			nsum /= binsize;
			kval.push_back(ksum);
			numberval.push_back(nsum);
		}
	}

	// int max_k = *std::max_element(kval.begin(),kval.end());
	// vector<int> kk(max_k);
	// vector<double> nn(max_k);
	// vector<int> divisor(max_k);
	// for(int r = 0; r < kval.size(); r++){
	// 	nn[kval[r]] += numberval[r];

	// }

	auto k_minmax = std::minmax_element(kval.begin(),kval.end());
	auto number_minmax = std::minmax_element(numberval.begin(),numberval.end());

	kmin = *k_minmax.first;
	kmax = *k_minmax.second;
	numbermin = *number_minmax.first;
	numbermax = *number_minmax.second;

	int num_bins = 100;
	float bin_width = (kmax - kmin) / num_bins;

	k = mglData(kval.size());
	number = mglData(kval.size());

	for(int i = 0; i < kval.size(); i++){
		k.a[i] = kval[i];
		number.a[i] = numberval[i];
	}
}

void Plotter::plotEval(){
	combinedControl();
	// control();
	// spectrum();
	// densityMap();
	// contour();
	// vortices();
	

}

void Plotter::combinedControl(){
	mglGraph gr;

	string filename = dirname + "/control_" + stepsString + ".png";

	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(FONT_SIZE);
	gr.SetQuality(3);
	gr.Title(title.c_str());


	// SPECTRUM
	gr.SubPlot(2,2,0);
	// gr.Title("Spectrum");
	gr.SetRange('x',kmin,kmax);
	gr.SetRange('y',numbermin,numbermax);
	gr.SetCoor(11); // log-log-coordinates

	gr.Axis();
	gr.Label('x',"|k|",0); gr.Label('y',"n",0);

	gr.FPlot("x^(-2)","k");
	gr.AddLegend("k^(-2)","k");
	gr.FPlot("x^(-4.66)","r");
	gr.AddLegend("k^(-4.66)","r");
	gr.FPlot("x^(-4)","b");
	gr.AddLegend("k^(-4)","b");

	gr.Plot(k,number," .");
	gr.Legend();
	gr.SetCoor(0);


	// PHASE + VORTICES
	gr.SubPlot(2,2,2);
	// gr.Title("Phase and Vortices");
	gr.SetRange('x',0,eval.data.meta.grid[0]);
	gr.SetRange('y',0,eval.data.meta.grid[1]);
	gr.SetRange('z',phase);
	gr.SetRange('c',phase);

	// gr.SetTicks('c',M_PI,0,0,"\\pi");
	double val[]={-M_PI, -M_PI/2, 0, M_PI/2, M_PI};
	gr.SetTicksVal('c', mglData(5,val), "-\\pi\n-\\pi/2\n0\n\\pi/2\n\\pi");
	gr.Axis();
	gr.Colorbar(">");
	gr.Label('x',"x [points]",0); gr.Label('y',"y [points]",0);
	gr.Dens(phase);
	gr.Plot(vortex_x,vortex_y," .w");
	// gr.Adjust("c");
	gr.SetTicks('c');

	// DENSITY SURFACE
	gr.SubPlot(2,2,1);
	// gr.Title("Density Surface");
	gr.SetRange('x',-xrange,xrange);
	gr.SetRange('y',-yrange,yrange);
	gr.SetRange('z',density);
	gr.SetRange('c',density);
	gr.Rotate(40,40);
	gr.Box();
	gr.Axis();
	gr.Colorbar(">");
	gr.Label('x',"x [\\mu m]",0); gr.Label('y',"y [\\mu m]",0);
	
	gr.Surf(density);

	// DENSITY + CONTOUR
	gr.SubPlot(2,2,3);
	// gr.Title("Density and Contour");
	gr.SetRange('x',0,eval.data.meta.grid[0]);
	gr.SetRange('y',0,eval.data.meta.grid[1]);
	gr.SetRange('z',density);
	gr.SetRange('c',density);
	gr.Axis();
	gr.Label('x',"x [points]",0); gr.Label('y',"y [points]",0);
	gr.Colorbar(">");
	gr.Dens(density);
	gr.Plot(contour_x,contour_y," .w");

	// gr.ShowImage("eog",true);

	gr.WritePNG(filename.c_str(),"Control",false);

}


void Plotter::control(){

	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	string filename = dirname + "/Control_" + stepsString + ".png";

	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(FONT_SIZE);
	gr.SetQuality(3);
	gr.Title(title.c_str());
	// gr.Title(title.c_str());
	// gr.Alpha(true);



	// data.use_abs=false;
	gr.SetRange('x',-xrange,xrange);
	gr.SetRange('y',-yrange,yrange);
	gr.SetRange('z',phase);
	gr.SetRange('c',phase);

	gr.SubPlot(2,2,0,"^_<>");
	gr.Rotate(40,40);
	gr.Box();
	gr.Axis();
	gr.Label('x',"x",0); gr.Label('y',"y",0);
	gr.Surf(phase);

	gr.SubPlot(2,2,2,"_<>^");
	gr.Axis();
	gr.Colorbar("_");
	gr.Label('x',"x",0); gr.Label('y',"y",0);
	gr.Dens(phase);


	// data.use_abs=true;
	gr.SetRange('x',-xrange,xrange);
	gr.SetRange('y',-yrange,yrange);
	gr.SetRange('z',density);
	gr.SetRange('c',density);

	gr.SubPlot(2,2,1,"^_<>");
	// gr.Light(true);
	gr.Rotate(40,40);
	gr.Box();
	gr.Axis();
	gr.Label('x',"x",0); gr.Label('y',"y",0);
	gr.Colorbar(">");
	gr.Surf(density);

	gr.SubPlot(2,2,3,"^_<>");
	gr.Axis();
	gr.Label('x',"x",0); gr.Label('y',"y",0);
	gr.Colorbar(">");
	gr.Dens(density);

	// gr.ShowImage("eog",true);

	gr.WritePNG(filename.c_str(),"Control",false);

}

void Plotter::spectrum(){



	mglGraph gr;

	gr.SetMarkSize(0.7);
	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(FONT_SIZE);
	gr.SetQuality(3);
	gr.SubPlot(1,1,0,"_");
	gr.Title(title.c_str());
	gr.SetRange('x',kmin,kmax);
	// gr.SetRange('x',*number_minmax.first,*number_minmax.second);
	gr.SetRange('y',numbermin,numbermax);
	gr.SetCoor(11); // log-log-coordinates

	gr.Axis();
	gr.Label('x',"|k|",0); gr.Label('y',"n",0);

	// gr.SetFunc("lg(x)","lg(y)");
	gr.FPlot("x^(-2)","k");
	gr.AddLegend("k^(-2)","k");
	gr.FPlot("x^(-4.66)","r");
	gr.AddLegend("k^(-4.66)","r");
	gr.FPlot("x^(-4)","b");
	gr.AddLegend("k^(-4)","b");

	gr.Plot(k,number," .");
	gr.Legend();

	string name = dirname + "/Spectrum_" + stepsString + ".png";

	gr.WritePNG(name.c_str(),"Spectrum",false);
}

void Plotter::contour(){

	mglGraph gr;

	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(FONT_SIZE);
	gr.SetQuality(3);
	gr.SubPlot(1,1,0,"<>_");
	gr.Title(title.c_str());
	// gr.Title(title.c_str());

	// gr.SetRange('x',-opt.min_x,opt.min_x);
	// gr.SetRange('y',-opt.min_y,opt.min_y);
	gr.SetRange('z',density);
	gr.SetRange('c',density);
	gr.SetRange('x',0,eval.data.meta.grid[0]);
	gr.SetRange('y',0,eval.data.meta.grid[1]);

	gr.Axis();
	gr.Label('x',"x [points]",0); gr.Label('y',"y [points]",0);
	gr.Colorbar(">I");
	gr.Dens(density);
	gr.Plot(contour_x,contour_y," .w");

	string name = dirname + "/Contour_" + stepsString + ".png";

	gr.WritePNG(name.c_str(),"Contour",false);

}

void Plotter::densityMap(){

	mglGraph gr;

	string filename = dirname + "/Density_" + stepsString + ".png";

	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(FONT_SIZE);
	gr.SetQuality(3);
	gr.SubPlot(1,1,0,"_");
	gr.Title(title.c_str());

	gr.SetRange('x',0,eval.data.meta.grid[0]);
	gr.SetRange('y',0,eval.data.meta.grid[1]);
	gr.SetRange('z',densitymap);
	gr.SetRange('c',densitymap);

	gr.Axis();
	gr.Label('x',"x",0); gr.Label('y',"y",0);
	gr.Dens(densitymap);

	gr.WritePNG(filename.c_str(),"Density",false);
}

void Plotter::vortices(){

	mglGraph gr;

	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(FONT_SIZE);
	gr.SetQuality(3);
	gr.SubPlot(1,1,0,"<_");
	gr.Title(title.c_str());

	gr.SetRange('x',0,eval.data.meta.grid[0]);
	gr.SetRange('y',0,eval.data.meta.grid[1]);
	gr.SetRange('c',phase);

	gr.Axis();
	gr.Label('x',"x",0); gr.Label('y',"y",0);
	gr.Colorbar(">F");
	gr.Dens(phase);
	// gr.Plot(v_x,v_y," #xw");
	gr.Plot(vortex_x,vortex_y," .w");

	string name = dirname + "/Vortices_" + stepsString + ".png";

	gr.WritePNG(name.c_str(),"Vortices",false);
}

