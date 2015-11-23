#ifndef BH3_FILEMANAGEMENT_H__
#define BH3_FILEMANAGEMENT_H__

#include <string>
#include <bh3observables.h>


using namespace std;

bool initialize_evaluation_dir(string &dirname, const string &dirsuffix, const list<string> &files);

void save_averages(const string & dirname, const PathOptions &opt, int p, const vector<Bh3Evaluation::Averages> &av);
void save_averages(const string & dirname, const PathOptions &opt, int p, const vector<Bh3Evaluation::Averages> &av, const vector<Bh3Evaluation::Averages> &std_dev);
void append_fields(const string & dirname, const PathOptions &opt, int p, const Grid &phase, const Grid &zeros, const ComplexGrid &data);
void append_vortices(const string & dirname, const PathOptions &opt, int p, const Bh3Evaluation::PathResults &pres, const Bh3Evaluation::Averages &av, bool save0 = false);

void save_binary_averages(const string & dirname, const PathOptions &options, int p, const vector<Bh3Evaluation::Averages> &av, const vector<Bh3Evaluation::PathResults> &pres);
string start_binary_averages(const string & dirname, const PathOptions &options, int p, int32_t size);
void append_binary_averages(const string & filename, const Bh3Evaluation::Averages &av, const Bh3Evaluation::PathResults &pres);
void load_binary_averages(const string & filename, PathOptions &options, vector<Bh3Evaluation::Averages> &av, vector<Bh3Evaluation::PathResults> &pres);
list<string> get_average_file_list(const list<string> args);

#endif
