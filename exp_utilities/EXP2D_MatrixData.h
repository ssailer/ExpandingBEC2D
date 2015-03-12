#ifndef EXP2D_MATRIXDATA_H__
#define EXP2D_MATRIXDATA_H__

#define EIGEN_FFTW_DEFAULT

#include <iostream>
#include <complex>
#include <EXP2D_tools.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/FFT>

using namespace std;
using namespace Eigen;

class MatrixData {
    public:

    class MetaData {
        public:
        MetaData() : steps(0), time(0), coord(2), grid(2), spacing(2), initCoord(2), initSpacing(2), samplesize(0) { for(int i = 0; i < 9; i++){array[i] = 0;}}
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
        inline void convertToDimensionless();
        inline void convertFromDimensionless();
    
        double time;
        int steps, samplesize;

        double Ag, OmegaG;
               
        vector<int> grid;

        vector<double> coord;
        vector<double> initCoord;

        vector<double> initSpacing;
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
    inline void increment(const double extTime,const double factorX,const double factorY);
    inline void setStep(const int &extStep);
    inline void setCoord(const vector<double> &extCoord);
    inline void setMeta(const MetaData &extMeta);

    inline void fftForward();
    inline void fftForward(MatrixXcd &DATA);

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


inline void MatrixData::update(const double &extTime,const int &extSteps,const vector<double> &coord){
    meta.time = extTime;
    meta.steps = extSteps;
    meta.coord[0] = coord[0];
    meta.coord[1] = coord[1];
    meta.spacing[0] = meta.coord[0] * 2 / meta.grid[0];
    meta.spacing[1] = meta.coord[1] * 2 / meta.grid[1];
    meta.dataToArray();
}

inline void MatrixData::increment(const double extTime,const double factorX,const double factorY){
    meta.time += extTime;
    meta.steps++;
    meta.coord[0] = factorX * meta.initCoord[0];
    meta.coord[1] = factorY * meta.initCoord[1];
    meta.spacing[0] = meta.coord[0] * 2 / meta.grid[0];
    meta.spacing[1] = meta.coord[1] * 2 / meta.grid[1];
    meta.dataToArray();
}

inline void MatrixData::MetaData::convertToDimensionless(){

    const double m = 87 * 1.66 * 1.0e-27;
    const double hbar = 1.054 * 1.0e-22;    
    Ag = 2 * initCoord[0] / grid[0];
    OmegaG = hbar / ( m * Ag * Ag);

    initCoord[0] /= Ag;
    initCoord[1] /= Ag;
    dataToArray();
    // opt.ITP_step *= opt.OmegaG;
    // opt.RTE_step *= opt.OmegaG;
    // opt.omega_x *= 2.0 * M_PI / opt.OmegaG;
    // opt.omega_y *= 2.0 * M_PI / opt.OmegaG;
    // opt.dispersion_x *= 2.0 * M_PI / opt.OmegaG;
    // opt.dispersion_y *= 2.0 * M_PI / opt.OmegaG;
}

inline void MatrixData::MetaData::convertFromDimensionless(){
    initCoord[0] *= Ag;
    initCoord[1] *= Ag;
    time /= OmegaG;
    time *= 1000.0; // conversion to ms
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

inline void MatrixData::fftForward(){

    for(int j = 0; j < wavefunction.size(); ++j){
        FFT<double> fft;
        MatrixXcd in = wavefunction[j];
        MatrixXcd out = MatrixXcd(meta.grid[0],meta.grid[1]);
        
        for (int k = 0; k < in.rows(); k++) {
            RowVectorXcd tmpOut = RowVectorXcd(meta.grid[0]);
            fft.fwd(tmpOut, in.row(k));
            out.row(k) = tmpOut;
        }
        
        for (int k = 0; k < in.cols(); k++) {
            VectorXcd tmpOut = VectorXcd(meta.grid[1]);
            fft.fwd(tmpOut, out.col(k));
            out.col(k) = tmpOut;
        }
        wavefunction[j] = out;
    }
}

inline void MatrixData::fftForward(MatrixXcd &DATA){

    // for(int j = 0; j < wavefunction.size(); ++j){
        FFT<double> fft;
        // MatrixXcd in = wavefunction[j];
        MatrixXcd out_temp = MatrixXcd(DATA.cols(),DATA.rows());
        
        for (int k = 0; k < DATA.rows(); k++) {
            RowVectorXcd tmpOut = RowVectorXcd(DATA.cols());
            fft.fwd(tmpOut, DATA.row(k));
            out_temp.row(k) = tmpOut;
        }
        
        for (int k = 0; k < DATA.cols(); k++) {
            VectorXcd tmpOut = VectorXcd(DATA.rows());
            fft.fwd(tmpOut, out_temp.col(k));
            out_temp.col(k) = tmpOut;
        }
        // wavefunction[j] = out;
        DATA = out_temp;
    // }
}



#endif // EXP2D_MATRIXDATA_H__