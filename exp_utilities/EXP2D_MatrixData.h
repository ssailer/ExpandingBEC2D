#ifndef EXP2D_MATRIXDATA_H__
#define EXP2D_MATRIXDATA_H__

#define EIGEN_FFTW_DEFAULT
#define EIGEN_VECTORIZE
#define EIGEN_DONT_PARALLELIZE

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
        MetaData() : steps(0), time(0), coord(2), grid(2), spacing(2), initCoord(2), initSpacing(2), samplesize(0), array(14), isDimensionless(false), Ag(0), OmegaG(0) {}
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
        inline int size();
        inline void print();
        inline void dataToArray();
        inline void arrayToData();
        inline void convertToDimensionless();
        inline void convertFromDimensionless();
    
        double time;
        int steps, samplesize;
        bool isDimensionless;

        double Ag, OmegaG;
               
        vector<int> grid;

        vector<double> coord;
        vector<double> initCoord;

        vector<double> initSpacing;
        vector<double> spacing;
    
        vector<double> array;   
    };

    class FFT_internal {
    public:
        FFT_internal(vector<MatrixXcd>* r_Wavefunction, MetaData* r_meta) : wavefunction(r_Wavefunction), meta(r_meta) {
            eigenFFT = new FFT<double>;
            eigenFFT->SetFlag(Eigen::FFT<double>::Unscaled);
        }
        ~FFT_internal(){
            delete eigenFFT;
        }
        inline void Forward();
        inline void Forward_X();
        inline void Forward_Y();
        inline void Backward();
        inline void Backward_X();
        inline void Backward_Y();
        inline void Forward(MatrixXcd &DATA);

        FFT<double>* eigenFFT;
        vector<MatrixXcd>* wavefunction;
        MetaData* meta; 
    };

    vector<MatrixXcd> wavefunction;
    MetaData meta;
    FFT_internal fft;

    MatrixData() : fft(&wavefunction,&meta) {}
    MatrixData(const MatrixData &cData) : wavefunction(cData.wavefunction), meta(cData.meta), fft(&wavefunction,&meta) {}
    // MatrixData(int size) : wavefunction(size), meta() {}
    // MatrixData(int x, int y) {for(int i = 0; wavefunction.size(); i++){ wavefunction[i] = MatrixXcd::Zero(x,y);}}
    // MatrixData(int size, int x, int y) : MatrixData(size), MatrixData(x,y) {}
    inline MatrixData(const MetaData &extMeta);
    inline MatrixData(const int &samplesize,const int &gridx, const int &gridy,const double &extTime, const int &extStep, const double &xsize, const double &ysize);

    inline void updateTo(const double &extTime, const int &extSteps,const vector<double> &coordFactor);
    inline void convertToDimensionlessUnits();
    inline void convertToPhysicalUnits();
    inline void setMatrix(const vector<MatrixXcd> &extWavefct);
    inline void setTime(const double &extTime);
    inline void increment(const double extTime,const double factorX,const double factorY);
    inline void setStep(const int &extStep);
    inline void setCoord(const vector<double> &extCoord);
    inline void setMeta(const MetaData &extMeta);
    inline void checkedCopy(MatrixData* &e);
    inline void resize();



    inline double hX();
    inline double hY();
    inline int getStep();
    inline vector<MatrixXcd> getMatrix();    
    inline MetaData getMeta();
};

inline MatrixData::MatrixData(const MetaData &m) : meta(m), wavefunction(m.samplesize), fft(&wavefunction,&meta) {    
        for(int i = 0; i < wavefunction.size(); i++){
            wavefunction[i] = MatrixXcd::Zero(meta.grid[0],meta.grid[1]);
        }
    }   


inline MatrixData::MatrixData(const int &samplesize,const int &gridx, const int &gridy,const double &extTime, const int &extStep, const double &xsize, const double &ysize) : fft(&wavefunction,&meta) {
        
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

inline void MatrixData::checkedCopy(MatrixData* &e){
    if(meta.grid[0] > e->meta.grid[0] || meta.grid[1] > e->meta.grid[1]){
        wavefunction.resize(e->wavefunction.size());
        for(int i = 0; i < wavefunction.size(); i++){
            wavefunction[i] = MatrixXcd::Zero(meta.grid[0],meta.grid[1]);
            wavefunction[i].block(meta.grid[0]/2 - e->meta.grid[0]/2, meta.grid[1]/2 - e->meta.grid[1]/2 , e->meta.grid[0], e->meta.grid[1]) = e->wavefunction[i];
        }
        meta.spacing[0] = e->meta.spacing[0];
        meta.spacing[1] = e->meta.spacing[1];
        // meta.initSpacing[0] = e->meta.initSpacing[0];
        // meta.initSpacing[1] = e->meta.initSpacing[1];
        meta.time = e->meta.time;
        meta.steps = e->meta.steps;
        meta.coord[0] = e->meta.coord[0] * meta.grid[0] / e->meta.grid[0];
        meta.coord[1] = e->meta.coord[1] * meta.grid[1] / e->meta.grid[1];
        // meta.initCoord[0] = e->meta.initCoord[0] * meta.grid[0] / e->meta.grid[0];
        // meta.initCoord[1] = e->meta.initCoord[1] * meta.grid[1] / e->meta.grid[1];
        meta.dataToArray();

        cout << "Checked Copy Control e" << endl
            << e->meta.spacing[0] << " " << endl
            << e->meta.initSpacing[0] << " " << endl
            << e->meta.time << " " << endl
            << e->meta.steps << " " << endl
            << e->meta.coord[0] << " " << endl
            << e->meta.initCoord[0] << " " << endl
            << endl;

        cout << "Checked Copy Control" << endl
            << meta.spacing[0] << " " << endl
            << meta.initSpacing[0] << " " << endl
            << meta.time << " " << endl
            << meta.steps << " " << endl
            << meta.coord[0] << " " << endl
            << meta.initCoord[0] << " " << endl
            << endl;
    } else {
        meta = e->meta;
        wavefunction = e->wavefunction;
    }
}

inline void MatrixData::resize(){
    int edge = 400;
    int grid0 = meta.grid[0] + edge;
    int grid1 = meta.grid[1] + edge;
    double ratio0 = (double)grid0 / (double)meta.grid[0];
    double ratio1 = (double)grid0 / (double)meta.grid[1];


    // for(int i = 0; i < wavefunction.size(); i++){
    //     MatrixXcd tmp_wavefunction = MatrixXcd::Zero(grid0,grid1);
    //     tmp_wavefunction.block(edge/2+1, edge/2+1, meta.grid[0],meta.grid[1]) = wavefunction[i];
    //     wavefunction[i] = tmp_wavefunction;
    // }
    MatrixXcd tmp_wavefunction = MatrixXcd::Zero(grid0,grid1);
    for(int i = 0; i < meta.grid[0]; i++){
        for(int j = 0; j < meta.grid[1]; j++){
            tmp_wavefunction(i+edge/2,j+edge/2) = wavefunction[0](i,j);
        }
    }
    wavefunction[0] = tmp_wavefunction;

    meta.initCoord[0] *= ratio0;
    meta.initCoord[1] *= ratio1;
    meta.coord[0] *= ratio0;
    meta.coord[1] *= ratio1;
    meta.grid[0] = grid0;
    meta.grid[1] = grid1;

    meta.dataToArray();
}


inline void MatrixData::updateTo(const double &extTime,const int &extSteps,const vector<double> &coord){
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

    if(Ag != 0.0 && OmegaG != 0.0 && isDimensionless == false){
        initCoord[0] /= Ag;
        initCoord[1] /= Ag;
        initSpacing[0] /= Ag;
        initSpacing[1] /= Ag;
        coord[0] /= Ag;
        coord[1] /= Ag;
        spacing[0] /= Ag;
        spacing[1] /= Ag;
        time *= OmegaG / 1000.0;
        isDimensionless = true;
        dataToArray();    
    } else {
        cout << "Cannot convert, MetaData is already dimensionless." << endl;
    }
}

inline void MatrixData::MetaData::convertFromDimensionless(){
    if(Ag != 0.0 && OmegaG != 0.0 && isDimensionless == true){
        initCoord[0] *= Ag;
        initCoord[1] *= Ag;
        initSpacing[0] *= Ag;
        initSpacing[1] *= Ag;
        coord[0] *= Ag;
        coord[1] *= Ag;
        spacing[0] *= Ag;
        spacing[1] *= Ag;
        time /= OmegaG;
        time *= 1000.0; // conversion to ms
        isDimensionless = false;
        dataToArray();
    } else {
        cout << "Cannot convert, MetaData is not dimensionless." << endl;
    }
}

inline void MatrixData::convertToDimensionlessUnits(){
    if(meta.Ag != 0.0 && meta.OmegaG != 0.0){
        meta.convertToDimensionless();
        for(int k = 0; k < wavefunction.size(); k++){  
            wavefunction[k] *= complex<double>(meta.Ag,0.0);
        }
    } else {
        cout << "Cannot convert MatrixData's units, Ag and OmegaG are not set!" << endl;
    }
}

inline void MatrixData::convertToPhysicalUnits(){
    if (meta.Ag != 0.0 && meta.OmegaG != 0.0){
        meta.convertFromDimensionless();
        for(int k = 0; k < wavefunction.size(); k++){  
            wavefunction[k] /= complex<double>(meta.Ag,0.0);
        }
    } else {
        cout << "Cannot convert MatrixData's units, Ag and OmegaG are not set!" << endl;
    }
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



inline double * MatrixData::MetaData::data(){
    return array.data();
}

inline int MatrixData::MetaData::size(){
    return array.size();
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
    array[9] = initCoord[0];
    array[10] = initCoord[1];
    array[11] = initSpacing[0];
    array[12] = initSpacing[1];
    array[13] = isDimensionless;
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
    initCoord[0] = array[9] ;
    initCoord[1] = array[10];
    initSpacing[0] = array[11];
    initSpacing[1] = array[12];
    isDimensionless = array[13];
}

inline void MatrixData::FFT_internal::Forward(){

    for(int j = 0; j < wavefunction->size(); ++j){
        // FFT<double> fft;
        // fft.SetFlag(Eigen::FFT<double>::Unscaled);
        // MatrixXcd in = wavefunction[j];
        MatrixXcd out = MatrixXcd(meta->grid[0],meta->grid[1]);
        
        for (int k = 0; k < wavefunction->at(j).rows(); k++) {
            RowVectorXcd tmpOut = RowVectorXcd(meta->grid[0]);
            eigenFFT->fwd(tmpOut, wavefunction->at(j).row(k));
            out.row(k) = tmpOut;
        }
        
        for (int k = 0; k < wavefunction->at(j).cols(); k++) {
            VectorXcd tmpOut = VectorXcd(meta->grid[1]);
            eigenFFT->fwd(tmpOut, out.col(k));
            wavefunction->at(j).col(k) = tmpOut;
        }
        // wavefunction[j] = out;
    }
}

inline void MatrixData::FFT_internal::Forward_X(){

    for(int j = 0; j < wavefunction->size(); ++j){

        // MatrixXcd in = wavefunction[j];
        // MatrixXcd out = MatrixXcd(meta.grid[0],meta.grid[1]);
        
        #pragma omp parallel for
        for (int k = 0; k < wavefunction->at(j).rows(); k++) {
            RowVectorXcd tmpOut = RowVectorXcd(meta->grid[0]);
            eigenFFT->fwd(tmpOut, wavefunction->at(j).row(k));
            wavefunction->at(j).row(k) = tmpOut;
        }
        
        // wavefunction[j] = out;
    }
}

inline void MatrixData::FFT_internal::Forward_Y(){

    for(int j = 0; j < wavefunction->size(); ++j){
        // FFT<double> fft;
        // fft.SetFlag(Eigen::FFT<double>::Unscaled);
        // MatrixXcd in = wavefunction[j];
        // MatrixXcd out = MatrixXcd(meta.grid[0],meta.grid[1]);
        
        #pragma omp parallel for
        for (int k = 0; k < wavefunction->at(j).cols(); k++) {
            VectorXcd tmpOut = VectorXcd(meta->grid[1]);
            eigenFFT->fwd(tmpOut, wavefunction->at(j).col(k));
            wavefunction->at(j).col(k) = tmpOut;
        }
        // wavefunction[j] = out;
    }
}

inline void MatrixData::FFT_internal::Backward(){

    for(int j = 0; j < wavefunction->size(); ++j){
        // FFT<double> fft;
        // fft.SetFlag(Eigen::FFT<double>::Unscaled);
        // MatrixXcd in = wavefunction[j];
        MatrixXcd out = MatrixXcd(meta->grid[0],meta->grid[1]);
        
        for (int k = 0; k < wavefunction->at(j).rows(); k++) {
            RowVectorXcd tmpOut = RowVectorXcd(meta->grid[0]);
            eigenFFT->inv(tmpOut, wavefunction->at(j).row(k));
            out.row(k) = tmpOut;
        }
        
        for (int k = 0; k < wavefunction->at(j).cols(); k++) {
            VectorXcd tmpOut = VectorXcd(meta->grid[1]);
            eigenFFT->inv(tmpOut, out.col(k));
            wavefunction->at(j).col(k) = tmpOut;
        }
        // wavefunction[j] = out;
    }
}

inline void MatrixData::FFT_internal::Backward_X(){

    for(int j = 0; j < wavefunction->size(); ++j){
        // FFT<double> fft;
        // fft.SetFlag(Eigen::FFT<double>::Unscaled);
        // MatrixXcd in = wavefunction[j];
        // MatrixXcd out = MatrixXcd(meta.grid[0],meta.grid[1]);
        
        #pragma omp parallel for
        for (int k = 0; k < wavefunction->at(j).rows(); k++) {
            RowVectorXcd tmpOut = RowVectorXcd(meta->grid[0]);
            eigenFFT->inv(tmpOut, wavefunction->at(j).row(k));
            wavefunction->at(j).row(k) = tmpOut;
        }

        // wavefunction[j] = out;
    }
}

inline void MatrixData::FFT_internal::Backward_Y(){

    for(int j = 0; j < wavefunction->size(); ++j){
        // FFT<double> fft;
        // fft.SetFlag(Eigen::FFT<double>::Unscaled);
        // MatrixXcd in = wavefunction[j];
        // MatrixXcd out = MatrixXcd(meta.grid[0],meta.grid[1]);
        
        #pragma omp parallel for
        for (int k = 0; k < wavefunction->at(j).cols(); k++) {
            VectorXcd tmpOut = VectorXcd(meta->grid[1]);
            eigenFFT->inv(tmpOut, wavefunction->at(j).col(k));
            wavefunction->at(j).col(k) = tmpOut;
        }
        // wavefunction[j] = out;
    }
}



inline void MatrixData::FFT_internal::Forward(MatrixXcd &DATA){

    // for(int j = 0; j < wavefunction.size(); ++j){
        // FFT<double> fft;
        // MatrixXcd in = wavefunction[j];
        MatrixXcd out_temp = MatrixXcd(DATA.rows(),DATA.cols());
        
        for (int k = 0; k < DATA.rows(); k++) {
            RowVectorXcd tmpOut = RowVectorXcd(DATA.cols());
            eigenFFT->fwd(tmpOut, DATA.row(k));
            out_temp.row(k) = tmpOut;
        }
        
        for (int k = 0; k < DATA.cols(); k++) {
            VectorXcd tmpOut = VectorXcd(DATA.rows());
            eigenFFT->fwd(tmpOut, out_temp.col(k));
            out_temp.col(k) = tmpOut;
        }
        // wavefunction[j] = out;
        DATA = out_temp;
    // }
}

inline void MatrixData::MetaData::print(){
    cout << endl
         << "MatrixData is dimensionless:" << isDimensionless << endl
         << "At step: " << steps << " with a time t = " << time << endl
         << "Gridsize: " << grid[0] << "x" << grid[1] << endl << endl
         << "Inital Values are: " << endl
         << "x = " << initCoord[0] << " y = " << initCoord[1] << endl
         << "h_x = " << initSpacing[0] << " h_y " << initSpacing[1] << endl << endl
         << "Actual Values are: " << endl
         << "x = " << coord[0] << " y = " << coord[1] << endl
         << "h_x = " << spacing[0] << " h_y " << spacing[1] << endl << endl;
}



#endif // EXP2D_MATRIXDATA_H__