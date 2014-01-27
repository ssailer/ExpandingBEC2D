#ifndef BH3BINARYFILE_H__
#define BH3BINARYFILE_H__

#include <list>
#include <string>
#include <fstream>
#include <stdint.h>
#include <vector>
#include <complex>
#include <hdf5.h>


using namespace std;

class ComplexGrid;
class RealGrid;

typedef struct {
	double timestepsize;
	double U;
	double N;
	int32_t grid[4];	
	double klength[3];
	vector<double> delta_t;
	vector<double> g;
} PathOptions;

class Bh3BinaryFile
{
	public:
		enum mode {in, out};
		
	protected:
		string filename;
		hid_t h5_file;
		hid_t h5_complexgriddata, h5_realgriddata;
		PathOptions options;
		fstream stream;
		mode m;
		int32_t num_snapshots_c;
		int32_t num_snapshots_r;
	public:
		Bh3BinaryFile(const string &file, const PathOptions &opt, mode nm);
		~Bh3BinaryFile();
		
		void close();
		
		bool append_snapshot(double time, const vector<ComplexGrid> &k, bool compression = false);
		bool append_snapshot(double time, const vector<RealGrid> &k, bool compression = false);		
		bool get_snapshot(double &time, vector<ComplexGrid> &k, int num);
		bool get_snapshot(double &time, vector<RealGrid> &k, int num);
		
		inline int get_num_snapshots_c() const {return num_snapshots_c;}
		inline int get_num_snapshots_r() const {return num_snapshots_r;}
		inline const PathOptions & get_options() const {return options;}
	protected:
		Bh3BinaryFile() {}

};

bool initialize_binary_dir(string &dirname, const PathOptions &opt);
bool initialize_path_options(int argc, char** argv, PathOptions &opt, vector<double> &snapshot_times,void *gen_times(vector<double> &snapshot_times)=NULL, bool verbose = false);
list<string> get_file_list(const list<string> args);
void convert_bin_file(const string &filename, bool forward, bool compression = false);

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

inline bool operator== (const PathOptions &p1, const PathOptions &p2)
{
	return ((p1.timestepsize == p2.timestepsize) &&
			(p1.grid[0] == p2.grid[0]) &&
			(p1.grid[1] == p2.grid[1]) &&
			(p1.grid[2] == p2.grid[2]) &&
			(p1.grid[3] == p2.grid[3]) &&
			(p1.klength[0] == p2.klength[0]) &&
			(p1.klength[1] == p2.klength[1]) &&
			(p1.klength[2] == p2.klength[2]) &&
			(p1.U == p2.U) &&
			(p1.N == p2.N) &&
			(p1.delta_t == p2.delta_t) &&
			(p1.g == p2.g));			
}

inline bool operator!= (const PathOptions &p1, const PathOptions &p2)
{
	return !(p1 == p2);
}

// helper functions for reading and writing
inline void write(ostream &o, int32_t n)
{
	o.write((char *) &n, sizeof(int32_t));
}

inline void write(ostream &o, uint32_t n)
{
	o.write((char *) &n, sizeof(uint32_t));
}

inline void write(ostream &o, float n)
{
	o.write((char *) &n, sizeof(float));
}

inline void write(ostream &o, double n)
{
	o.write((char *) &n, sizeof(double));
}

inline void read(istream &i, int32_t &n)
{
	i.read((char *) &n, sizeof(int32_t));
}

inline void read(istream &i, uint32_t &n)
{
	i.read((char *) &n, sizeof(uint32_t));
}

inline void read(istream &i, float &n)
{
	i.read((char *) &n, sizeof(float));
}

inline void read(istream &i, double &n)
{
	i.read((char *) &n, sizeof(double));
}

template<class T>
inline void write(ostream &o, const complex<T> &c)
{
	o.write((const char *) &c, sizeof(complex<T>));
}

template<class T>
inline void read(istream &i, complex<T> &c)
{
	i.read((char *) &c, sizeof(complex<T>));
}

template<class T>
inline void write(ostream &o, const vector<T> &v)
{
	int32_t size = v.size();
	write(o, size);
	for(int i = 0; i < size; i++)
	{
		write(o, v[i]);
	}
}

template<class T>
inline void read(istream &i, vector<T> &v)
{
	int32_t size;
	read(i, size);
	v.resize(size);
	for(int pos = 0; pos < size; pos++)
	{
		read(i, v[pos]);
	}
}

template<class T>
inline void write(ostream &o, const list<T> &l)
{
	list<T> temp_list = l;
	int32_t size = temp_list.size();
	write(o, size);
	while(!temp_list.empty())
	{
		write(o, temp_list.front());
		temp_list.pop_front();
	}
}

template <class T>
inline void read(istream &i, list<T> &l)
{
	int32_t size;
	read(i, size);
	l.clear();
	for(int pos = 0; pos < size; pos++)
	{
		T temp;
		read(i, temp);
		l.push_back(temp);
	}
}

inline void write(ostream &stream, const PathOptions &options)
{
	write(stream, options.N);
	write(stream, options.U);
	write(stream, options.timestepsize);
	write(stream, options.delta_t);
	write(stream, options.g);
	write(stream, options.grid[0]);
	write(stream, options.grid[1]);
	write(stream, options.grid[2]);
	write(stream, options.grid[3]);	
	write(stream, options.klength[0]);
	write(stream, options.klength[1]);
	write(stream, options.klength[2]);
}

inline void read(istream &stream, PathOptions &options)
{
	read(stream, options.N);
	read(stream, options.U);
	read(stream, options.timestepsize);
	read(stream, options.delta_t);
	read(stream, options.g);
	read(stream, options.grid[0]);
	read(stream, options.grid[1]);
	read(stream, options.grid[2]);
	read(stream, options.grid[3]);
	read(stream, options.klength[0]);
	read(stream, options.klength[1]);
	read(stream, options.klength[2]);
}

inline ostream & operator<< (ostream &o, const PathOptions &options)
{
	o << "Number of particles: " << options.N << endl;
	o << "Timestepsize: " << options.timestepsize << endl;
	o << "Interaction strength: " << options.U << endl;
	o << "Delta t: ";
	for(int i = 0; i < options.delta_t.size(); i++)
		o << options.delta_t[i] << ", ";
	o << endl;
	o << "g-Matrix: ";
	for(int i = 0; i < options.g.size(); i++)
		o << options.g[i] << ", ";
	o << endl;
	o << "Field components " << options.grid[0] << endl;
	o << "Grid dimensions: " << options.grid[1] << ", " << options.grid[2] << ", " << options.grid[3] << endl;
	o << "K-Length: " << options.klength[0] << ", " << options.klength[1] << ", " << options.klength[2] << endl;
	return o;
}

#endif
