#include <EXP2D_observables.h>

using namespace std;
using namespace Eigen;


 Averages::Averages() {};
 Averages::~Averages() {};

void Averages::saveData(vector<MatrixXcd> &wavefctVec,Options &externalopt,int &external_snapshot_time){
		opt = externalopt;
		snapshot_time = external_snapshot_time;
		PsiVec.resize(wavefctVec.size());

		#pragma parallel for
		for(int k = 0; k < wavefctVec.size(); k++){
			PsiVec[k] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

			for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){		
			PsiVec[k](0,i,j,0) = wavefctVec[k](i,j);}}
		}
}

void Averages::saveData(MatrixXcd &wavefct,Options &externalopt,int &external_snapshot_time){
		opt = externalopt;
		snapshot_time = external_snapshot_time;
		PsiVec.resize(1);

		PsiVec[0] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

		for(int i = 0; i < opt.grid[1]; i++){
			for(int j = 0; j < opt.grid[2]; j++){		
				PsiVec[0](0,i,j,0) = wavefct(i,j);
			}
		}		
}

void Averages::evaluateData(){
	vector<Evaluation> avResult(PsiVec.size());
		
	for(int k = 0; k < PsiVec.size(); k++){
		avResult[k] = Evaluation(3*opt.grid[1]);
		avResult[k] = evaluate(PsiVec[k]);
		// avResult + evaluate(PsiVec[k]);
	}


	totalResult = Evaluation(3*opt.grid[1]);
		for(int k = 0; k < PsiVec.size(); k++){
		totalResult.k += avResult[k].k;
		totalResult.number += avResult[k].number;
		totalResult.Ekin += avResult[k].Ekin;
		totalResult.particle_count += avResult[k].particle_count;
		totalResult.healing_length += avResult[k].healing_length;
		}
	
	totalResult.k /= PsiVec.size();
	totalResult.number /= PsiVec.size();
	totalResult.Ekin /= PsiVec.size();
	totalResult.particle_count /= PsiVec.size();
	totalResult.healing_length /= PsiVec.size();
	// if(opt.samplesize == PsiVec.size()){
	// 	// for(int i = 0; totalResult.number.size();i++){
	// 	// 	totalResult.number(i) /= (double)opt.samplesize;
	// 	// }
	// totalResult /= (double)opt.samplesize;
	// }else{cout << "Error in evaluationData() member."<<endl;}

	// string pathnumber = "0";
	// plot(pathnumber,snapshot_time,avResult[0]);
	// pathnumber = "1";
	// plot(pathnumber,snapshot_time,avResult[1]);
}

void Averages::plotTotalResult(){
	plot(snapshot_time,totalResult);
}

void Averages::plot(const int &snapshot_time,Evaluation &eval)
{  
	// stringstream d;
	// d << dirname << "/" << "spectrum" << "/";
	// string dir = d.str();
	// system((string("mkdir -p ") + dir).c_str());
	string filename = "Spectrum-Step-" + to_string(snapshot_time); 
	plotspectrum(filename,eval);

	// cout << "Particle count: " << eval.particle_count << "  " << "Total energy: " << eval.Ekin << endl;
}

Evaluation Averages::evaluate(ComplexGrid &data){

Evaluation ares = Evaluation(3*opt.grid[1]);
// R-Space
double h_x = 2.*opt.min_x/opt.grid[1]; // FIXME not expanding! 
double h_y = 2.*opt.min_y/opt.grid[2]; 

double threshold = abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;


double volume = 0.0;  
for(int i=opt.grid[1];i<opt.grid[1];i++){
    for(int j=opt.grid[2];j<opt.grid[2];j++){
    	if(abs2(data(0,i,j,0))>threshold){
      		volume += h_x*h_y*(abs2(data(0,i,j,0))+abs2(data(0,i+1,j,0))+abs2(data(0,i,j+1,0))+abs2(data(0,i+1,j+1,0)))/4.0;
      	}
    }
}


// K-Space
ComplexGrid::fft(data, data);

ArrayXd divisor(ares.number.size());
divisor.setZero();

vector<vector<double>> kspace;

kspace.resize(2);
for(int d = 0; d < 2; d++)
	{
		// set k-space
		kspace[d].resize(opt.grid[d+1]);
		for (int i=0; i<opt.grid[d+1]/2; i++)
			kspace[d][i] = opt.klength[d]*sin( M_PI*((double)i)/((double)opt.grid[d+1]) );

		for (int i=opt.grid[d]/2; i<opt.grid[d+1]; i++)
			kspace[d][i] = opt.klength[d]*sin( M_PI*((double)(-opt.grid[d+1]+i))/((double)opt.grid[d+1]) );
	}

double kwidth2[3];

for(int i = 0; i < 3; i++)
	kwidth2[i] = (opt.grid[i+1] == 1) ? 0 : opt.klength[i]*opt.klength[i];

double index_factor = (ares.number.size() - 1) / sqrt(kwidth2[0] + kwidth2[1] + kwidth2[2]);

for(int x = 0; x < data.width(); x++)
{
	for (int y = 0; y < data.height(); y++)
	{
		for (int z = 0; z < data.depth(); z++)
		{
		double k = sqrt(kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]);
		Coordinate<int32_t> c = data.make_coord(x,y,z);
		int index = index_factor * k;
		ares.k(index) += k;
		divisor(index)++;
		double number = abs2(data(0,c));
		ares.number(index) += number;
		ares.particle_count += number;
		ares.Ekin += number * k * k;
		}
	}
}

ares.healing_length = 1 / sqrt(ares.particle_count * opt.g / volume);

#pragma omp parallel for schedule(guided,1)
for(int l = 0; l < ares.number.size(); l++)
{
	if(divisor[l] == 0)
		divisor[l] = 1;
	    ares.number[l] /= divisor[l];
   		ares.k[l] /= divisor[l];
}


return ares;
}

