#include <inttypes.h>
#include <stdio.h>

#include <iostream>
#include <unistd.h>
#include <cstdlib>
// #include <cstring>
#include <string>
#include <cmath>
#include <complex>
#include <omp.h>
#include <sys/stat.h>
#include <dirent.h>

#include <bh3defaultgrid.h>

#include <EXP2D_MatrixData.h>
#include <main.h>
#include <EXP2D_tools.h>
#include <EXP2D_binaryfile.h>
#include <EXP2D_rk4.hpp>
#include <EXP2D_runner.hpp>
#include <EXP2D_evaluation.h>
#include <plot_with_mgl.h>
#include <EXP2D_startgrids.h>


using namespace std;

int main(int argc, char *argv[])
{	
	InitMain initMain(argc,argv);	
	initMain.printInitVar();
	Options opt = initMain.getOptions();
	MatrixData* data = new MatrixData(initMain.getMeta());
	setGridToSinus(data,initMain.getOptions());

	Eval* initEval = new Eval(*data,opt);
	initEval->process();
	initEval->save();

	Plotter* initPlot = new Plotter(*initEval,opt);
	initPlot->plotEval();
	delete initPlot, initEval, data;

    return 0;
}