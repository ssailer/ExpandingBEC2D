#ifndef EXP2D_MATRIXDATA_H__
#define EXP2D_MATRIXDATA_H__

#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class MetaData {
public:
    MetaData() : stepState(0), timeState(0), coordinateBoundaries(2), grid(2), spacing(2), samplesize(0) {
        for(int i = 0; i < 9; i++)
            array[i] = 0;
    }
    inline double* data();
    inline void dataToArray();
    inline void arrayToData();

    double timeState;
    int stepState, samplesize;
    vector<double> coordinateBoundaries;
    vector<int> grid;
    vector<double> spacing;

    double array[9];
};

inline double* MetaData::data(){
    return array;
}

inline void MetaData::dataToArray(){
    array[0] = timeState;
    array[1] = (double)stepState;
    array[2] = (double)samplesize;
    array[3] = (double)grid[0];
    array[4] = (double)grid[1];
    array[5] = coordinateBoundaries[0];
    array[6] = coordinateBoundaries[1];
    array[7] = spacing[0];
    array[8] = spacing[1];
}

inline void MetaData::arrayToData(){
    timeState = array[0];
    stepState = (int)array[1];
    samplesize = (int)array[2];
    grid[0] = (int)array[3];
    grid[1] = (int)array[4];
    coordinateBoundaries[0] = array[5];
    coordinateBoundaries[1] = array[6];
    spacing[0] = array[7];
    spacing[1] = array[8];
}

class MatrixData {
public:
    vector<MatrixXcd> wavefunction;
    MetaData meta;

    MatrixData() : wavefunction(0), meta() {}
    inline MatrixData(const MetaData &extMeta);
    inline MatrixData(const int &samplesize,const int &gridx, const int &gridy,const double &extTime, const int &extStep, const double &xsize, const double &ysize);

    inline void setMatrix(const vector<MatrixXcd> &extWavefct);
    inline void setTime(const double &extTime);
    inline void setStep(const int &extStep);
    inline void setCoord(const vector<double> &extCoord);
    inline void setMeta(const MetaData &extMeta);

    inline double getGridXSpacing();
    inline double getGrixYSpacing();
    inline int getStep();
    inline vector<MatrixXcd> getMatrix();    
    inline MetaData getMeta();
};



inline MatrixData::MatrixData(const int &samplesize,const int &gridx, const int &gridy,const double &extTime, const int &extStep, const double &xsize, const double &ysize) {
        
    wavefunction.resize(samplesize);
    for(int i = 0; wavefunction.size(); ++i)
        wavefunction[i] = MatrixXcd(gridx,gridy);

    meta.grid.resize(2);
    meta.grid[0] = gridx;
    meta.grid[1] = gridy;

    meta.coordinateBoundaries.resize(2);
    meta.coordinateBoundaries[0] = xsize;
    meta.coordinateBoundaries[1] = ysize;

    meta.spacing.resize(2);
    meta.spacing[0] = xsize * 2 / gridx;
    meta.spacing[1] = ysize * 2 / gridy;

    meta.timeState = extTime;
    meta.stepState = extStep;


}

inline MatrixData::MatrixData(const MetaData &extMeta){

    meta = extMeta;        
    wavefunction.resize(meta.samplesize);
    for(int i = 0; wavefunction.size(); ++i)
        wavefunction[i] = MatrixXcd(meta.grid[0],meta.grid[1]);
}

inline void MatrixData::setMatrix(const vector<MatrixXcd> &extWavefct){
    wavefunction = extWavefct;
}

inline void MatrixData::setTime(const double &extTime){
    meta.timeState = extTime;
}

inline void MatrixData::setStep(const int &extStep){
    meta.stepState = extStep;
}

inline void MatrixData::setCoord(const vector<double> &extCoord){
    if( extCoord.size() == 2){
       meta.coordinateBoundaries = extCoord;
    } else {
        cerr << "Error in MatrixData::setCoord, wrong vector size";
    }
}

inline double MatrixData::getGridXSpacing(){
    return meta.spacing[0];
}

inline double MatrixData::getGrixYSpacing(){
    return meta.spacing[1];
}

inline int MatrixData::getStep(){
    return meta.stepState;
}

inline vector<MatrixXcd> MatrixData::getMatrix(){
    return wavefunction;
}

inline void MatrixData::setMeta(const MetaData &extMeta){
    meta = extMeta;
}

inline MetaData MatrixData::getMeta(){
    return meta;
}



#endif // EXP2D_MATRIXDATA_H__