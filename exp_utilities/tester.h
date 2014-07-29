class matrixData {
public:
    matrixData() : wavefunction(0), timeState(0), coordinateBoundaries(2) {}
    matrixData(const int &samplesize){
        wavefunction.resize(samplesize);
        timeState = 0;
        coordinateBoundaries.resize(2);
    };
    vector<MatrixXcd> * data() { return &wavefunction ;  }
    void setTime(const int &tmpTime){ 
        timeState = tmpTime;
    }
    void setCoord(const vector<double> &tmpCoord){
        if(coordinateBoundaries.size() == tmpCoord.size())
            coordinateBoundaries = tmpCoord;
        else
            cerr << "Wrong size in matrixData::setCoord()" << endl;
    }
    vector<MatrixXcd> getWavefct(){   return wavefunction;   }
    double getTime(){ return timeState;  }
    vector<double> getCoord(){    return coordinateBoundaries;    }

private:    
    vector<MatrixXcd> wavefunction;
    double timeState;
    vector<double> coordinateBoundaries;
};