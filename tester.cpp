#include <EXP2D_MatrixData.h>

int main( int argc, char** argv) 
{	

	MetaData meta;
	meta.grid[0] = 512;
	meta.grid[1] = 512;
	meta.samplesize = 48;
	meta.coordinateBoundaries[0] = 5;
	meta.coordinateBoundaries[1] = 5;
	meta.timeState = 0;
	meta.stepState = 0;
	meta.spacing[0] = meta.coordinateBoundaries[0] * 2 / meta.grid[0];
	meta.spacing[0] = meta.coordinateBoundaries[1] * 2 / meta.grid[1];

	MatrixData data;
	data.setMeta(meta);
	data.meta.dataToArray();
	cout << data.meta.data()[3] << endl;

	return 0;
}