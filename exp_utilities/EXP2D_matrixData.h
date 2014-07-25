#ifndef EXP2D_MATRIXDATA_H__
#define EXP2D_MATRIXDATA_H__


class matrixData {
public:

    vector<MatrixXcd> wavefunction;
    double timeState;
    vector<double> coordinateBoundaries;

    matrixData() : wavefunction(0), timeState(0), coordinateBoundaries(2) {}
    matrixData(const int &samplesize,const int &gridx, const int &gridy,const int &tmpTime, const int &xsize, const int &ysize);
};

matrixData::matrixData(const int &samplesize,const int &gridx, const int &gridy,const double &extTime, const int &xsize, const int &ysize) {
        
        wavefunction.resize(samplesize);
        for(int i = 0; wavefunction.size(); ++i)
            wavefunction[i] = MatrixXcd(gridx,gridy);

        coordinateBoundaries.resize(2);
        coordinateBoundaries[0] = xsize;
        coordinateBoundaries[1] = ysize;

        timeState = extTime;
    }

matrixData::setMatrix(const vector<MatrixXcd> &extWavefct){
       wavefunction = extWavefct;
}

matrixData::setTime(const double &extTime)

#endif // EXP2D_MATRIXDATA_H__