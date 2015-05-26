#ifndef EXP2D_BINARTFILE_H__
#define EXP2D_BINARTFILE_H__

#include <list>
#include <string>
#include <fstream>
#include <stdint.h>
#include <vector>
#include <complex>
#include <hdf5.h>
#include <EXP2D_tools.h>
#include <EXP2D_MatrixData.h>
#include <EXP2D_evaluation.h>
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

  bool appendSnapshot(const string &name, MatrixData * const &pData, Options const &options);
  // bool appendSnapshot(const string &name, double time, const vector<RealGrid> &k);
  // bool appendSnapshot(const string &name, int snapShotTime, vector<ComplexGrid> &data, MatrixData::MetaData &meta, Options &options);

  bool appendEval( Eval &results, Options const & opt );
  // bool appendDocString(const string &group, const string &docstring, double time);

  bool getSnapshot(const string &name, int time, MatrixData* &pData, Options &options);
  bool getLatestSnapshot(const string &name, MatrixData* &pData, Options &options);
  // bool getSnapshot(const string &name, double time, vector<RealGrid> &k);
  // bool getSnapshot(const string &name, int snapShotTime, vector<ComplexGrid> &data,MatrixData::MetaData &meta, Options &options);

  bool getEval(int snapShotTime, Eval &results, Options &options);

  const vector<int> getTimeList() const {return time_list;}

  // const Options & getOptions() const {return options;}

protected:

  binaryFile() {}
  bool checkTime(int snapShotTime);

  void writeMatrixData(const string &name, MatrixData * const &pData, Options const &options );
  void readMatrixData(string const &name, MatrixData* &pData, Options &options);

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
          (p1.grid[0] == p2.grid[0]) &&
          (p1.grid[1] == p2.grid[1]) &&
          (p1.grid[2] == p2.grid[2]) &&
          (p1.grid[3] == p2.grid[3]) &&
          (p1.klength[0] == p2.klength[0]) &&
          (p1.klength[1] == p2.klength[1]) &&
          (p1.klength[2] == p2.klength[2]) &&
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

// // helper functions for reading and writing
// inline void write(ostream &o, int32_t n)
// {
//   o.write((char *) &n, sizeof(int32_t));
// }

// inline void write(ostream &o, uint32_t n)
// {
//   o.write((char *) &n, sizeof(uint32_t));
// }

// inline void write(ostream &o, float n)
// {
//   o.write((char *) &n, sizeof(float));
// }

// inline void write(ostream &o, double n)
// {
//   o.write((char *) &n, sizeof(double));
// }

// inline void read(istream &i, int32_t &n)
// {
//   i.read((char *) &n, sizeof(int32_t));
// }

// inline void read(istream &i, uint32_t &n)
// {
//   i.read((char *) &n, sizeof(uint32_t));
// }

// inline void read(istream &i, float &n)
// {
//   i.read((char *) &n, sizeof(float));
// }

// inline void read(istream &i, double &n)
// {
//   i.read((char *) &n, sizeof(double));
// }

// template<class T>
// inline void write(ostream &o, const complex<T> &c)
// {
//   o.write((const char *) &c, sizeof(complex<T>));
// }

// template<class T>
// inline void read(istream &i, complex<T> &c)
// {
//   i.read((char *) &c, sizeof(complex<T>));
// }

// template<class T>
// inline void write(ostream &o, const vector<T> &v)
// {
//   int32_t size = v.size();
//   write(o, size);
//   for(int i = 0; i < size; i++)
//     {
//       write(o, v[i]);
//     }
// }

// template<class T>
// inline void read(istream &i, vector<T> &v)
// {
//   int32_t size;
//   read(i, size);
//   v.resize(size);
//   for(int pos = 0; pos < size; pos++)
//     {
//       read(i, v[pos]);
//     }
// }

// template<class T>
// inline void write(ostream &o, const list<T> &l)
// {
//   list<T> temp_list = l;
//   int32_t size = temp_list.size();
//   write(o, size);
//   while(!temp_list.empty())
//     {
//       write(o, temp_list.front());
//       temp_list.pop_front();
//     }
// }

// template <class T>
// inline void read(istream &i, list<T> &l)
// {
//   int32_t size;
//   read(i, size);
//   l.clear();
//   for(int pos = 0; pos < size; pos++)
//     {
//       T temp;
//       read(i, temp);
//       l.push_back(temp);
//     }
// }

inline void write(ostream &stream, const Options &opt)
{
  write(stream, opt.N);
  write(stream, opt.stateInformation[0]);
  write(stream, opt.stateInformation[1]);
  write(stream, opt.omega_x);
  write(stream, opt.omega_y);
  write(stream, opt.dispersion_x);
  write(stream, opt.dispersion_y);
  write(stream, opt.min_x);
  write(stream, opt.min_y);
  // write(stream, opt.t_abs);
  write(stream, opt.exp_factor);
  write(stream, opt.ITP_step);
  write(stream, opt.RTE_step);
  write(stream, opt.n_it_RTE);
  write(stream, opt.samplesize);
  write(stream, opt.vortexnumber);
  write(stream, opt.g);
  write(stream, opt.grid[0]);
  write(stream, opt.grid[1]);
  write(stream, opt.grid[2]);
  write(stream, opt.grid[3]);
  write(stream, opt.klength[0]);
  write(stream, opt.klength[1]);
  write(stream, opt.klength[2]);
}

inline void read(istream &stream, Options &opt)
{
  read(stream, opt.N);
  read(stream, opt.stateInformation[0]);
  read(stream, opt.stateInformation[1]);
  read(stream, opt.omega_x);
  read(stream, opt.omega_y);
  read(stream, opt.dispersion_x);
  read(stream, opt.dispersion_y);
  read(stream, opt.min_x);
  read(stream, opt.min_y);
  // read(stream, opt.t_abs);
  read(stream, opt.exp_factor);
  read(stream, opt.ITP_step);
  read(stream, opt.RTE_step);
  read(stream, opt.n_it_RTE);
  read(stream, opt.samplesize);
  read(stream, opt.vortexnumber);
  read(stream, opt.g);
  read(stream, opt.grid[0]);
  read(stream, opt.grid[1]);
  read(stream, opt.grid[2]);
  read(stream, opt.grid[3]);
  read(stream, opt.klength[0]);
  read(stream, opt.klength[1]);
  read(stream, opt.klength[2]);
}

inline ostream & operator<< (ostream &o, const Options &opt)
{
  o << "Number of particles: " << opt.N << endl;
  o << "Interaction strength: " << opt.g << endl;
  o << endl;
  o << "Grid dimensions: " << opt.grid[1] << ", " << opt.grid[2] << endl;
  o << "K-Length: " << opt.klength[0] << ", " << opt.klength[1] <<  endl;
  o << endl;
  o << "Potential Frequencies: " << opt.omega_x << ", " << opt.omega_y << endl;
  o << "Dispersion Frequencies: " << opt.dispersion_x << ", " << opt.dispersion_y << " with an overall expansion factor: " << opt.exp_factor << endl;
  o << endl;
  o << "RTE Stepsize: " << opt.RTE_step << " for " << opt.n_it_RTE << " Steps." << endl;
  o << "Samplesize: " << opt.samplesize << endl;
  o << endl;
  return o;
}

#endif