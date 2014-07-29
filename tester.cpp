#include <EXP2D_MatrixData.h>

int main( int argc, char** argv) 
{	

	MatrixData::MetaData meta;
	meta.grid[0] = 512;
	meta.grid[1] = 512;
	meta.samplesize = 8;
	meta.coord[0] = 5;
	meta.coord[1] = 5;
	meta.time = 0;
	meta.steps = 0;
	meta.spacing[0] = meta.coord[0] * 2 / meta.grid[0];
	meta.spacing[0] = meta.coord[1] * 2 / meta.grid[1];
	meta.dataToArray();

	MatrixData data(meta);
	return 0;
}