#ifndef EXP2D_BINARTFILE_H__
#define EXP2D_BINARTFILE_H__

#include <list>
#include <string>
#include <fstream>
#include <stdint.h>
#include <vector>
#include <complex>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <unistd.h>
#include <dirent.h>
#include <string.h>

#include "tools.h"
#include "matrixdata.h"
#include "evaluation.h"

#include <hdf5.h>
#include <eigen3/Eigen/Dense>



using namespace std;

class binaryFile
{
public:
  enum mode {in, out, append};

protected:
  string filename;
  hid_t h5_file;
  hid_t h5_timegroup;
  // Options options;
  fstream stream;
  mode m;
  vector<int> time_list;


public:
  binaryFile(const string &file, mode nm);
  ~binaryFile();

  void close();

  bool appendSnapshot(const string &name, shared_ptr<MatrixData> pData, Options const &options);
  // bool appendSnapshot(const string &name, double time, const vector<RealGrid> &k);
  // bool appendSnapshot(const string &name, int snapShotTime, vector<ComplexGrid> &data, MatrixData::MetaData &meta, Options &options);

  bool appendEval( shared_ptr<Eval> results, Options const & opt );
  // bool appendDocString(const string &group, const string &docstring, double time);

  bool getSnapshot(const string &name, int time, shared_ptr<MatrixData> pData, Options &options);
  bool getLatestSnapshot(const string &name, shared_ptr<MatrixData> pData, Options &options);
  // bool getSnapshot(const string &name, double time, vector<RealGrid> &k);
  // bool getSnapshot(const string &name, int snapShotTime, vector<ComplexGrid> &data,MatrixData::MetaData &meta, Options &options);

  bool getEval(int snapShotTime, Eval &results, Options &options);

  const vector<int> getTimeList() const {return time_list;}

  // const Options & getOptions() const {return options;}

protected:

  binaryFile() {}
  bool checkTime(int snapShotTime);

  void writeMatrixData(const string &name, shared_ptr<MatrixData> pData, Options const &options );
  void readMatrixData(string const &name, shared_ptr<MatrixData> pData, Options &options);

  void writeMeta(hid_t &h5_group,MatrixData::MetaData &meta );
  void readMeta(hid_t &h5_group,MatrixData::MetaData &meta);

  void writeOptions(hid_t &h5_group,Options const & options);
  void readOptions(hid_t &h5_group,Options &options);



};

// bool initialize_binary_dir(string &dirname, const Options &opt);
// bool initialize_path_options(int argc, char** argv, Options &opt, vector<double> &snapshot_times,void *gen_times(vector<double> &snapshot_times)=NULL, bool verbose = false);
// list<string> get_file_list(const list<string> args);

template <typename T>
inline bool operator== (const vector<T> &v1, const vector<T> &v2)
{
  if(v1.size() == v2.size())
    {
      bool result = true;
      for(int i = 0; i < v1.size(); i++)
        {
          result = result && (v1[i] == v2[i]);
        }
      return result;
    }
  else
    return false;
}

inline bool operator== (const Options &p1, const Options &p2)
{
  return ((p1.N == p2.N) &&
          // (p1.grid[0] == p2.grid[0]) &&
          (p1.grid[1] == p2.grid[1]) &&
          (p1.grid[2] == p2.grid[2]) &&
          // (p1.grid[3] == p2.grid[3]) &&
          // (p1.klength[0] == p2.klength[0]) &&
          // (p1.klength[1] == p2.klength[1]) &&
          // (p1.klength[2] == p2.klength[2]) &&
          (p1.stateInformation[0] == p2.stateInformation[0]) &&
          (p1.stateInformation[1] == p2.stateInformation[1]) &&
          (p1.omega_x == p2.omega_x) &&
          (p1.omega_y == p2.omega_y) &&
          (p1.dispersion_x == p2.dispersion_x) &&
          (p1.dispersion_y == p2.dispersion_y) &&
          (p1.g == p2.g) &&
  		  (p1.min_x == p2.min_x) &&
  		  (p1.min_y == p2.min_y) &&
  		  // (p1.t_abs == p2.t_abs) &&
  		  (p1.exp_factor == p2.exp_factor) &&
  		  (p1.ITP_step == p2.ITP_step) &&
  		  (p1.RTE_step == p2.RTE_step) &&
  		  (p1.n_it_RTE == p2.n_it_RTE) &&
  		  (p1.samplesize == p2.samplesize) &&
  		  (p1.vortexnumber == p2.vortexnumber) &&
  		  (p1.runmode == p2.runmode) &&
  		  (p1.config == p2.config) &&
  		  (p1.workingdirectory == p2.workingdirectory));
}

inline bool operator!= (const Options &p1, const Options &p2)
{
  return !(p1 == p2);
}

#endif