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
	const int vortex_size = eval.vlist[0].size();
	int index;





	// // ellipse cover for phase

	// vector<int> c_x;
	// vector<int> c_y;
	// 
	// for(int i = 0; i < eval.data.meta.grid[0]; i++ ){
	// 	for(int j = 0; j < eval.data.meta.grid[1]; j++ ){
	// 		int x = i ;
	// 		int y = j ;
	// 		float cond = eval.ellipse.coef(0) * x * x + eval.ellipse.coef(1) * x * y + eval.ellipse.coef(2) * y * y + eval.ellipse.coef(3) * x + eval.ellipse.coef(4) * y + eval.ellipse.coef(5);
			
	// 		if(cond >= thresh /*&& cond >= -thresh*/){
	// 			// cerr << "x " << x << " y " << endl;
	// 			// c_x.push_back(x);
	// 			// c_y.push_back(y);
	// 			phase(x,y) = 0.0;
	// 		}
	// 	}
	// }
	// cover_x = mglData(c_x.size());
	// cover_y = mglData(c_y.size());
	// for(int i = 0; i < c_x.size(); i++){
	// 	cover_x.a[i] = c_x[i];
	// 	cover_y.a[i] = c_y[i];
	// }

		// Phase (with elliptic cover) and Density
	density = mglData(n,m);
	phase = mglData(n,m);
	densitymap = mglData(n,m);

	for(int x = 0; x < n; x++){
		for(int y = 0; y < m; y++){
			index = x + n * y;
			float cond = eval.ellipse.coef(0) * x * x + eval.ellipse.coef(1) * x * y + eval.ellipse.coef(2) * y * y + eval.ellipse.coef(3) * x + eval.ellipse.coef(4) * y + eval.ellipse.coef(5);
			density.a[index] = abs2(eval.data.wavefunction[0](x,y));
			if(cond <= 0.0){
				phase.a[index] = arg(eval.data.wavefunction[0](x,y));
			} else {
				phase.a[index] = 0.0;
			}
			densitymap.a[index] = eval.densityLocationMap[0](x,y);
		}
	}

		// // contour vectors
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
	for(list<VortexData>::const_iterator it = eval.vlist[0].begin(); it != eval.vlist[0].end(); ++it){
		vortex_x.a[index] = it->c.x();
		vortex_y.a[index] = it->c.y();
		index++;
	}


	// spectrum
	vector<double> kval;
	vector<double> numberval;
	map<double,double> spectrum;
    pair<map<double,double>::iterator,bool> ret;
    vector<double> tmpKval;
		
    for (int r = 0; r < eval.totalResult.number.size(); r++){
		if(eval.totalResult.k(r) != 0.0){
			if(eval.totalResult.number(r) != 0.0){
				ret = spectrum.insert(map<double,double>::value_type(eval.totalResult.k(r),eval.totalResult.number(r)));
				tmpKval.push_back(eval.totalResult.k(r));
				if(ret.second==false){
					cout << "Binning of spectrum failed, double value inserted." << endl;
				}

				// kval.push_back(k_int);
				// numberval.push_back(eval.totalResult.number(r));
			}
        }
	}

	auto tmpMinMax = std::minmax_element(tmpKval.begin(),tmpKval.end());



	int c = 1;
	double nsum = 0;	
	double min_value = *tmpMinMax.first;
	double max_value = *tmpMinMax.second;
	double min_log = log(min_value);
	double max_log = log(max_value);
	int nBins = sqrt(eval.data.meta.grid[0] * eval.data.meta.grid[0] + eval.data.meta.grid[1] * eval.data.meta.grid[1])/4;
	// double binMax = exp(min_value);
	double log_increment = (max_value - min_value) / nBins;
	
	double log_value = min_log + log_increment;
	double value = exp(log_value);
	
	vector<double> median;
	for(map<double,double>::const_iterator it = spectrum.begin(); it != spectrum.end(); ++it){

		if(it->first <= value){
			// nsum += it->second;
			// c++;
			median.push_back(it->second);
		} else {
			sort(median.begin(),median.end());
			int size = median.size();
			if(size%2 == 0){
				nsum = median[size/2];
			} else {
				nsum = (median[size/2] + median[size/2 +1]) /2;
			}
			// nsum /= c;
			numberval.push_back(nsum);
			double tmp_log = log_value;
			log_value += log_increment;
			double k = exp((tmp_log + log_value)/2);
			kval.push_back(k);
			value = exp(log_value);
			nsum = it->second;
			c = 1;
			median.clear();
		}
	}

	// ksum /= c;
	// nsum /= c;
	// kval.push_back(ksum);
	// numberval.push_back(nsum);

	// binsize = c1 - lastbin;
	// ksum /= binsize;
	// nsum /= binsize;
	// kval.push_back(ksum);
	// numberval.push_back(nsum);

	// int max_k = *std::max_element(kval.begin(),kval.end());
	// vector<int> kk(max_k);
	// vector<double> nn(max_k);
	// vector<int> divisor(max_k);
	// for(int r = 0; r < kval.size(); r++){
	// 	nn[kval[r]] += numberval[r];

	// }

	auto k_minmax = std::minmax_element(kval.begin(),kval.end());
	auto number_minmax = std::minmax_element(numberval.begin(),numberval.end());

	// kmin = *k_minmax.first;
	kmin = 0.5;
	kmax = *k_minmax.second;
	numbermin = *number_minmax.first;
	numbermax = *number_minmax.second;

	k = mglData(kval.size());
	number = mglData(kval.size());

	for(int i = 0; i < kval.size(); i++){
		if(kval[i] >= kmin){
			k.a[i] = kval[i];
			number.a[i] = numberval[i];
		}
	}

	// ellipse axis
	origin = mglPoint(eval.ellipse.center[0],eval.ellipse.center[1]);
	major_1 = mglPoint(eval.ellipse.center[0] + eval.ellipse.major * cos(eval.ellipse.angle),eval.ellipse.center[1] + eval.ellipse.major * sin(eval.ellipse.angle));
	minor_1 = mglPoint(eval.ellipse.center[0] - eval.ellipse.minor * sin(eval.ellipse.angle),eval.ellipse.center[1] + eval.ellipse.minor * cos(eval.ellipse.angle));

	major_2 = mglPoint(eval.ellipse.center[0] - eval.ellipse.major * cos(eval.ellipse.angle),eval.ellipse.center[1] - eval.ellipse.major * sin(eval.ellipse.angle));
	minor_2 = mglPoint(eval.ellipse.center[0] + eval.ellipse.minor * sin(eval.ellipse.angle),eval.ellipse.center[1] - eval.ellipse.minor * cos(eval.ellipse.angle));

	reg_1 = mglPoint(eval.punkte[0],eval.punkte[1]);
	reg_2 = mglPoint(eval.punkte[2],eval.punkte[3]);

	cerr << endl << eval.punkte[0] << " " << eval.punkte[1] << " " << eval.punkte[2] << " " << eval.punkte[3] << endl;

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

	// https://groups.google.com/forum/#!topic/mathgl/68VMLblLd7U
	// check this for plot margins, set by hand
	mglGraph gr;

	string filename = dirname + "/control_" + stepsString + ".png";

	gr.SetSize(IMAGE_SIZE,IMAGE_SIZE);
	gr.SetFontSize(FONT_SIZE);
	gr.SetQuality(3);
	gr.Title(title.c_str());


	// SPECTRUM
	gr.SubPlot(2,2,0);
	// gr.Title("Spectrum");
	gr.SetRange('x',kmin,kmax+2);
	gr.SetRange('y',numbermin,numbermax+2);
	gr.SetCoor(11); // log-log-coordinates

	gr.Axis();
	gr.Label('x',"|k|",0); gr.Label('y',"n",0);

	gr.FPlot("x^(-2)","k");
	gr.AddLegend("k^(-2)","k");
	// gr.FPlot("x^(-4)","b");
	// gr.AddLegend("k^(-4)","b");
	// gr.FPlot("x^(-4.66)","r");
	// gr.AddLegend("k^(-4.66)","r");

	string leg_s = "\\alpha = " + to_string(eval.steigung) + " \\pm " + to_string(eval.fehler);
	gr.Plot(k,number," .2b");
	gr.Line(reg_1,reg_2,"r2");
	gr.AddLegend(leg_s.c_str(),"r");
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
	// gr.SetMarkSize(0.01);
	// gr.Plot(cover_x,cover_y," .w");
	gr.Plot(vortex_x,vortex_y," #<w");
	gr.Plot(vortex_x,vortex_y," <k");
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
	gr.Line(major_2,major_1,"W2");
	gr.Line(minor_2,minor_1,"H2");
	// gr.Line(origin-major_2,origin,"W2");
	// gr.Line(origin-minor_2,origin,"H2");

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

