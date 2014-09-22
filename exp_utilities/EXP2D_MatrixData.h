#ifndef EXP2D_MATRIXDATA_H__
#define EXP2D_MATRIXDATA_H__

#include <iostream>
#include <complex>
#include <EXP2D_tools.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class MatrixData {
    public:

    class MetaData {
        public:
        MetaData() : steps(0), time(0), coord(2), grid(2), spacing(2), samplesize(0) { for(int i = 0; i < 9; i++){array[i] = 0;}}
        // MetaData( const MetaData &m) : coord(2), grid(2), spacing(2) {
        //     steps = m.steps;
        //     time = m.time;
        //     coord[0] = m.coord[0];
        //     coord[1] = m.coord[1];
        //     grid[0] = m.grid[0];
        //     grid[1] = m.grid[1];
        //     spacing[0] = m.grid[0];
        //     spacing[1] = m.spacing[1];
        //     samplesize = m.samplesize;
        //     arrayToData();
        // }
        inline double* data();
        inline void dataToArray();
        inline void arrayToData();
    
        double time;
        int steps, samplesize;
        vector<double> coord;
        vector<int> grid;
        vector<double> spacing;
    
        double array[9];    
    };

    vector<MatrixXcd> wavefunction;
    MetaData meta;

    MatrixData() {}
    // MatrixData(int size) : wavefunction(size), meta() {}
    // MatrixData(int x, int y) {for(int i = 0; wavefunction.size(); i++){ wavefunction[i] = MatrixXcd::Zero(x,y);}}
    // MatrixData(int size, int x, int y) : MatrixData(size), MatrixData(x,y) {}
    inline MatrixData(const MetaData &extMeta);
    inline MatrixData(const int &samplesize,const int &gridx, const int &gridy,const double &extTime, const int &extStep, const double &xsize, const double &ysize);

    inline void update(const double &extTime, const int &extSteps,const vector<double> &coordFactor);

    inline void setMatrix(const vector<MatrixXcd> &extWavefct);
    inline void setTime(const double &extTime);
    inline void setStep(const int &extStep);
    inline void setCoord(const vector<double> &extCoord);
    inline void setMeta(const MetaData &extMeta);

    inline double hX();
    inline double hY();
    inline int getStep();
    inline vector<MatrixXcd> getMatrix();    
    inline MetaData getMeta();
};

inline MatrixData::MatrixData(const MetaData &m) : meta(m), wavefunction(m.samplesize) {    
        for(int i = 0; i < wavefunction.size(); i++){
            wavefunction[i] = MatrixXcd::Zero(meta.grid[0],meta.grid[1]);
        }
    }   


inline MatrixData::MatrixData(const int &samplesize,const int &gridx, const int &gridy,const double &extTime, const int &extStep, const double &xsize, const double &ysize) {
        
    wavefunction.resize(samplesize);
    for(int i = 0; i < wavefunction.size(); i++)
        wavefunction[i] = MatrixXcd::Zero(gridx,gridy);

    meta.grid.resize(2);
    meta.grid[0] = gridx;
    meta.grid[1] = gridy;

    meta.coord.resize(2);
    meta.coord[0] = xsize;
    meta.coord[1] = ysize;

    meta.spacing.resize(2);
    meta.spacing[0] = xsize * 2 / gridx;
    meta.spacing[1] = ysize * 2 / gridy;

    meta.time = extTime;
    meta.steps = extStep;

    meta.dataToArray();
}


inline void MatrixData::update(const double &extTime,const int &extSteps,const vector<double> &coordFactor){
    meta.time = extTime;
    meta.steps = extSteps;
    meta.coord[0] *= coordFactor[0];
    meta.coord[1] *= coordFactor[1];
    meta.spacing[0] = meta.coord[0] * 2 / meta.grid[0];
    meta.spacing[1] = meta.coord[1] * 2 / meta.grid[1];
    meta.dataToArray();
}

inline void MatrixData::setMatrix(const vector<MatrixXcd> &extWavefct){
    wavefunction = extWavefct;
}

inline void MatrixData::setTime(const double &extTime){
    meta.time = extTime;
    meta.dataToArray();
}

inline void MatrixData::setStep(const int &extStep){
    meta.steps = extStep;
    meta.dataToArray();
}

inline void MatrixData::setCoord(const vector<double> &extCoord){
    if( extCoord.size() == 2){
       meta.coord = extCoord;
       meta.dataToArray();
    } else {
        cerr << "Error in MatrixData::setCoord, wrong vector size";
    }
}

inline double MatrixData::hX(){
    return meta.spacing[0];
}

inline double MatrixData::hY(){
    return meta.spacing[1];
}

inline int MatrixData::getStep(){
    return meta.steps;
}

inline vector<MatrixXcd> MatrixData::getMatrix(){
    return wavefunction;
}

inline void MatrixData::setMeta(const MetaData &extMeta){
    meta = MetaData(extMeta);
}

inline MatrixData::MetaData MatrixData::getMeta(){
    return meta;
}



inline double* MatrixData::MetaData::data(){
    return array;
}

inline void MatrixData::MetaData::dataToArray(){
    array[0] = time;
    array[1] = (double)steps;
    array[2] = (double)samplesize;
    array[3] = (double)grid[0];
    array[4] = (double)grid[1];
    array[5] = coord[0];
    array[6] = coord[1];
    array[7] = spacing[0];
    array[8] = spacing[1];
}

inline void MatrixData::MetaData::arrayToData(){
    time = array[0];
    steps = (int)array[1];
    samplesize = (int)array[2];
    grid[0] = (int)array[3];
    grid[1] = (int)array[4];
    coord[0] = array[5];
    coord[1] = array[6];
    spacing[0] = array[7];
    spacing[1] = array[8];
}



#endif // EXP2D_MATRIXDATA_H__