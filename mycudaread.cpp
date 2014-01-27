#include<iostream>
#include<fstream>
#include<cmath>
#include<complex>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <omp.h>

#include<bh3binaryfile.h>
#include<complexgrid.h>
#include<realgrid.h>
#include<bh3propagator.h>
#include<averageclass.h>
#include <wrapped_cuda_functions.h>
//#include<bh3save.h>

#include<coordinate.h>

#define NUM_SINGLE_PLOTS 999999
#define PNG_WIDTH 640 
#define PNG_HEIGHT 480
#define FIELD_WIDTH 6



using namespace std;


void single_shot (const string &dirname, ComplexGrid grid1, ComplexGrid grid2, ComplexGrid test, const double time);

int main (int argc, char **argv)
{
  uint plots[2] = {1,1}; // plot snapshot 1arg from file 2arg
	list<string> files;
	list<string> arguments;
	fstream file;
	PathOptions options;
	uint num_snapshots;
	string dirname;
	//init_cuda_device(0);
	
	if (argc == 1)
	{
		cout << "No snapshot-files given!" << endl;
		return 0;
		//arguments.push_back("bh3cuda_TS0.2_G4096_1_1_N40000_U0.005_KL2_2_2Run3");
	}

	for(int i = 1; i < argc; i++)
	{
		arguments.push_back(argv[i]);
	}
	
	files = get_file_list(arguments);
	dirname = arguments.front();
	
	//initialize_evaluation_dir(dirname, "_test", files);


	Bh3BinaryFile *bfile = new Bh3BinaryFile(files.front(), options, Bh3BinaryFile::in);
	options = bfile->get_options();
	num_snapshots = (uint) bfile->get_num_snapshots()/2;
	//Bh3BinaryFile test_bin ("jackfile.h5", options, Bh3BinaryFile::in);


	
	vector <double> times(num_snapshots, -1.0);
	for(int i = 0; i < num_snapshots; i++)
	{
		double test_time;
		vector<ComplexGrid> test_grid;
		
		bfile->get_snapshot(test_time, test_grid);
		times[i] = test_time;
		bfile->get_snapshot(test_time, test_grid); //second time is Omega_t
		
		
	}
	delete bfile;

	vector<string> f(files.size());
	// Dateiliste in Vektor kopieren
	int findex = 0;
	for(list<string>::const_iterator it = files.begin(); it != files.end(); ++it)
	{
		f[findex] = *it;
		findex++;
	}
	
	file.open("data/readout.dat", ios::out);
	//file.open((dirname + string("readout.dat")).c_str(), ios::out);
	
	cout << "# snapshots: " << num_snapshots << endl;
	
	int num_valid_paths = 0;
	
	AverageClass<vector<Bh3Evaluation::Averages> > av1;
	AverageClass<vector<Bh3Evaluation::Averages> > av2;
	
	#pragma omp parallel for reduction(+:num_valid_paths) schedule(dynamic,1)
	for(int path = 1; path <= f.size(); path++)
	{
		cout << "now processing path " << path << "...." << endl;
		// Snapshot-data (Zeit und Grids)
		double time;
		double Omega_t;
		vector<ComplexGrid> k1;
		vector<ComplexGrid> k2;
		
		/*test section for grid operations 
		RealGrid test(2, 6,10,1);
		RealGrid test2(1, 6,10,1);
		CudaRealGrid cuTest (2, 6,10,1);
		for (int i=0; i<6; i++)
		{
			for (int j = 0; j< 10; j++)
				test2(0,i,j,0) = sin(2*M_PI*j/10.);
		}
		
		test.set_component(test2, 0);
		test.set_component(test2, 1);
		
		
		copyHostToDevice2D(cuTest, test); 
		
		for (int i=0; i<6; i++)
		{
			for (int j = 0; j< 10; j++)
				cout << test(0,i,j,0) << " ";
			cout << endl;
		}
		cout<< endl;

		//CudaRealGrid::fft(cuTest, cuTest, CUFFT_FORWARD);
		//CudaRealGrid::fft(cuTest, cuTest, CUFFT_INVERSE);
		//RealGrid::fft_unnormalized(test,test,1);
		//RealGrid::fft_unnormalized(test,test,0);
		//test.at_fft(0,1,1,0) = complex<double> (0.66, 0.66);
		//test = cuTest;
		//RealGrid::fft_unnormalized(test,test,0);
		for (int i=0; i<6; i++)
		{
			for (int j = 0; j< 10; j++)
			{
				Coordinate <int32_t> co = test.make_coord(i,j,0);
				cout << test.at_fft(0, co) << " ";
			}
			cout << endl;
		}
		cout<< endl;

		
		//cout << test(0, 3,1,0) <<"   " << test(0,4,1,0) <<  endl;*/
		
		Bh3BinaryFile b(f[path-1], options, Bh3BinaryFile::in);
		
		if((options != b.get_options()) || (num_snapshots != b.get_num_snapshots()/2))
		{
			cout << "File parameters of file '" << f[path-1] << "' doesn't fit to parameters of first file!" << endl << "Ignoring this file!" << endl;
		}
		else
		{
			num_valid_paths++;
			vector<Bh3Evaluation::Averages> path_res1(num_snapshots);
			vector<Bh3Evaluation::Averages> path_res2(num_snapshots);
			vector<Bh3EvaluationFD::Averages> path_resFD(num_snapshots);
			
			for(int i = 0; i < num_snapshots; i++)
			{
			  
				b.get_snapshot(time, k1);
				b.get_snapshot(Omega_t, k2);
				//test_bin.append_snapshot(time, k1);
				//test_bin.append_snapshot(Omega_t, k2);
	
				if(times[i] != time)
				{
					cout << "Snapshots have different times as those of the first files!" << endl
					     << "Averages will be invalid!" << endl;
				}
				
				
				Bh3EvaluationFD ev(options,true);
				
				ev.setTime(time);
				ev.setOmega(Omega_t);
				ev.setData(k1,k2,Bh3EvaluationFD::RSpace);
				
				ev.calc_radial_averages();
				//ev.calc_g1();

				
				path_resFD[i] = ev.get_averageable_results();
				path_res1[i] = path_resFD[i].ares_1;
				path_res2[i] = path_resFD[i].ares_2;

			
				
				//path_res[i] = ev.get_averages(); 
				//path_res[i].time = time;

				
				/*vector <double> tempgrid;
				Grid evgrid = ev.get_density_field(1);
				for (int x = 0; x < evgrid.width(); x++)
					tempgrid.push_back(evgrid(x,0,0));
				path_res[i].spectrum1 = tempgrid;
				
				/*tempgrid.resize(0);
				evgrid = ev.get_density_field(2);
				for (int x = 0; x < evgrid.width(); x++)
					tempgrid.push_back(evgrid(x,0,0));
				path_res[i].spectrum2 = tempgrid;*/
				
				
				
				if (i == plots[0] && path == plots[1])
				{
					vector <double> tempgrid;
					//Grid evgrid = ev.get_density_field(1);
					//for (int x = 0; x < evgrid.width(); x++)
					//	tempgrid.push_back(evgrid(x,0,0));
					single_shot (string("data/single_shot.dat"), k1[0], k2[0], ev.get_grid1()[0], time); // Wichtig hier werden die Daten von bh3mycuda.cpp eingelesen
				}
				
			}
                        #pragma omp critical
			av1.average(path_res1);
			av2.average(path_res2);
			avFD.average(path_resFD);
			
		}
		
		
	}
//############################################### evaluation region #############################################	

      


//################# write in file ########################################################## 
	
	vector <Bh3Evaluation::Averages> means1 = av1.av();
	vector <Bh3Evaluation::Averages> means2 = av2.av();
	vector <Bh3EvaluationFD::Averages> meansFD = avFD.av();
	
	
	for(int t = 0; t < times.size(); t++)
	{
		stringstream str;
		//plot (string("/remote/lin-22/karl/plots/spectra"), options, t, means[t], times[t]);
		str << times[t] << "\t";  
		//str << means[t].E_int << "\t";
		str << meansFD[t].Sz_mean << "\t" << meansFD[t].qw_mean << endl;
		cout << str.str();
	}
	
	for(int t = 0; t < num_snapshots; t++)
	{
		//cout << means2D[t].num_vortices << endl;
		for(int x = 0; x < meansFD[t].veff2.size(); x++)
		{
			file << x << "\t" << meansFD[t].ares_1.number[x] << "\t" << meansFD[t].veff_i2[x] << "\t" << meansFD[t].veff_c2[x] << "\t" << meansFD[t].veff_i_wos2[x] << endl; 
		}
		file << endl << endl;
	}
    
	file.close();



    return 0;
  
  
}


void single_shot (const string &dirname, ComplexGrid grid1, ComplexGrid grid2, ComplexGrid test, const double time)
{
	cout << "single shot at time " << time << endl; 
	ofstream plotfile;
	plotfile.open(dirname.c_str(), ios::out);
	
	//ComplexGrid::fft(grid1, grid1, true);
	//ComplexGrid::fft(grid2, grid2, true);
	stringstream str;
	//cout << grid1(1,1,0) << grid2(1,1,0) << endl;
	for(int x = 0; x < grid1.width(); x++)
	{
		for(int y = 0; y < grid1.height(); y++)
		{
			for(int z = 0; z < grid1.depth(); z++)
			{
// 				double n1 = norm(grid1(x,y,0));  
// 				double n2 = norm(grid2(x,y,0));
// 				double phi1 = arg(grid1(x,y,0));
// 				double phi2 = arg(grid2(x,y,0)); 
			  //	str << x << "\t" << y << "\t";
				/*str << (cos(atan(n1/n2)))*(cos(atan(n1/n2))) << "\t" << phi1 - phi2 << "\t";
				//str<< norm(grid1(x,y,0)) << "\t" << norm(grid2(x,y,0)) << "\t";
				str<< n1+n2 << "\t" << phi1 + phi2 << "\t";*/
			  str << x << "\t" << y << "\t"<< real(grid1(x,y,z)) << "\t" << imag(grid1(x,y,z)) << "\t";
				//str<< norm(grid1(x,y,0)) << "\t" << norm(grid2(x,y,0)) << "\t";
				//str<< real(grid2(x,y,z)) << "\t" << imag(grid2(x,y,z))<< "\t" << real(test(x,y,z)) << "\t" << imag(test(x,y,z));
				str<< endl;
			}
		}
		//str<<endl;
	}
	
	plotfile << str.str();
	plotfile << endl << endl;
	plotfile.close();
}
  

