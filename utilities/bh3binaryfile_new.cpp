#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <unistd.h>
#include <dirent.h>
#include <string.h>

#include <bh3binaryfile.h>
#include <complexgrid.h>
#include <realgrid.h>
#include <hdf5.h>
#include <libconfig.h>

// Class implementation
Bh3BinaryFile::Bh3BinaryFile(const string &file, const PathOptions &opt, mode nm)
{
  filename = file;
  options = opt;
  m = nm;

  if(m == in || m == append)
    {
      struct stat buf;
      lstat(filename.c_str(), &buf);
      if(S_ISREG(buf.st_mode) && H5Fis_hdf5(filename.c_str())) // Nur normale Dateien ueberpruefen
        {
          if(m == in)
            h5_file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
          else if(m == append)
            h5_file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

          if(H5Aexists(h5_file, "PathOptions") > 0)
            {

              hid_t h5a_options;
              h5a_options = H5Aopen(h5_file, "PathOptions", H5P_DEFAULT);

              double temp_opt[10];

              H5Aread(h5a_options, H5T_IEEE_F64LE , temp_opt);

              //load PathOptions struct from file array. Don't forget to change appropriately when changing PathOptions struct
              options.timestepsize = temp_opt[0];
              options.N = temp_opt[1];
              options.U = temp_opt[2];
              for(int i = 0; i<4;i++)
                options.grid[i] = (uint32_t)temp_opt[i+3];
              for(int i= 0; i < 3; i++)
                options.klength[i] = temp_opt[i+7];

              H5Aclose(h5a_options);

              //load g_matrix from file
              if(H5Aexists(h5_file, "g_matrix") > 0)
                {
                  h5a_options = H5Aopen(h5_file, "g_matrix", H5P_DEFAULT);
                  int size = (int)(H5Aget_storage_size(h5a_options)/sizeof(double));
                  options.g.resize(size);

                  H5Aread(h5a_options, H5T_IEEE_F64LE , &(options.g.front()));

                  H5Aclose(h5a_options);
                }
              else
                options.g.resize(0);

              //load delta_t vector from file, if exists
              if(H5Aexists(h5_file, "delta_t") > 0)
                {
                  h5a_options = H5Aopen(h5_file, "delta_t", H5P_DEFAULT);
                  int size = (int)(H5Aget_storage_size(h5a_options)/sizeof(double));
                  options.delta_t.resize(size);

                  H5Aread(h5a_options, H5T_IEEE_F64LE , &(options.delta_t.front()));

                  H5Aclose(h5a_options);
                }
              else
                options.delta_t.resize(0);

              h5a_options = H5Aopen(h5_file, "timelist", H5P_DEFAULT);
              int size = (int)(H5Aget_storage_size(h5a_options)/sizeof(double));
              time_list.resize(size);
              H5Aread(h5a_options, H5T_IEEE_F64LE , &time_list.front());
              H5Aclose(h5a_options);

              if(m == append)
                {
                  if(options != opt)
                    cout << "WARNING for I/O operation on file: " << filename <<". Passed options do not fit file options in 'append' mode." << endl;
                }

            }
          else
            {
              cout << "WARNING for I/O operation on file: "<< filename << " occurred: No valid PathOptions attribute found in file." << endl;
            }
        }
      else
        {
          cout << "ERROR for I/O operation. File: "<< filename <<" is either not regular or not in hdf5 storage format" << endl;
        }
    }
  else
    {

      h5_file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      hid_t    h5a_options, dataspace;


      //load PathOptions struct to a reliable arrays. Don't forget to change appropriately when changing PathOptions struct
      double temp_opt[10];

      temp_opt[0] = options.timestepsize;
      temp_opt[1] = options.N;
      temp_opt[2] = options.U;
      for(int i = 0; i<4;i++)
        temp_opt[i+3] = (double)options.grid[i];
      for(int i= 0; i < 3; i++)
        temp_opt[i+7] = options.klength[i];

      //copy options to file
      hsize_t dimsf[1] = {10};
      dataspace = H5Screate_simple(1, dimsf, NULL);

      h5a_options = H5Acreate(h5_file, "PathOptions", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite (h5a_options, H5T_IEEE_F64LE, temp_opt);

      H5Aclose(h5a_options);
      H5Sclose(dataspace);

      //copy delta_t array to file if there is one
      if (options.delta_t.size()>0)
        {
          dimsf[0] = options.delta_t.size();
          dataspace = H5Screate_simple(1, dimsf, NULL);

          h5a_options = H5Acreate(h5_file, "delta_t", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
          H5Awrite (h5a_options, H5T_IEEE_F64LE, &(options.delta_t.front()));

          H5Aclose(h5a_options);
          H5Sclose(dataspace);
        }

      //copy g matrix to file

      if(options.g.size()>0)
        {
          dimsf[0] = options.g.size();
          dataspace = H5Screate_simple(1, dimsf, NULL);

          h5a_options = H5Acreate(h5_file, "g_matrix", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
          H5Awrite (h5a_options, H5T_IEEE_F64LE, &(options.g.front()));

          H5Aclose(h5a_options);
          H5Sclose(dataspace);
        }

    }
}

Bh3BinaryFile::~Bh3BinaryFile()
{
  close();
}

void Bh3BinaryFile::close()
{
  if(time_list.size()>0 && m==out)
    {
      hsize_t dimsf[] = {time_list.size()};
      hid_t dataspace = H5Screate_simple(1, dimsf, NULL);
      hid_t h5a = H5Acreate(h5_file, "timelist", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite (h5a, H5T_IEEE_F64LE, &time_list.front());
      H5Sclose(dataspace);
      H5Aclose(h5a);
    }
  else if(time_list.size()>0 && m == append)
    {
      if(H5Aexists(h5_file, "timelist"))
        {
          H5Adelete(h5_file, "timelist");
        }

      hsize_t dimsf[] = {time_list.size()};
      hid_t dataspace = H5Screate_simple(1, dimsf, NULL);
      hid_t h5a = H5Acreate(h5_file, "timelist", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite (h5a, H5T_IEEE_F64LE, &time_list.front());
      H5Sclose(dataspace);
      H5Aclose(h5a);
    }


  H5Fclose(h5_file);
  filename = "";
}

bool Bh3BinaryFile::check_time(double time)
{
  stringstream time_name;
  // time will always be formatted with two decimal digits
  time_name.setf(ios_base::fixed);
  time_name.precision(2);
  time_name << "/" << time;

  if(H5Lexists(h5_file, (time_name.str()).c_str(), H5P_DEFAULT))
    h5_timegroup = H5Gopen(h5_file, (time_name.str()).c_str(), H5P_DEFAULT);
  else
    {
      h5_timegroup = H5Gcreate(h5_file, (time_name.str()).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      time_list.push_back(time);
    }

  if(h5_timegroup >= 0)
    return true;
  else
    return false;
}


bool Bh3BinaryFile::append_snapshot(const string &name, double time, const vector<ComplexGrid> &k)
{
  if(m == in)
    {
      cout << "file "<< filename.c_str() << "is not in write mode" << endl;
      return false;
    }

  if(!check_time(time))
    {
      cout << "HDF5 Group error for time group: " << time << endl;
      H5Gclose(h5_timegroup);
      return false;
    }

  hid_t h5_CGgroup;

  if(H5Lexists(h5_timegroup, "ComplexGrid", H5P_DEFAULT))
    h5_CGgroup = H5Gopen(h5_timegroup, "ComplexGrid", H5P_DEFAULT);
  else
    h5_CGgroup = H5Gcreate(h5_timegroup, "ComplexGrid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


  if(H5Lexists(h5_CGgroup, name.c_str(), H5P_DEFAULT))
    {
      cout << "ComplexGrid " << name << " already exists for time " << time << ". I refuse to write." << endl;
      H5Gclose(h5_CGgroup);
      H5Gclose(h5_timegroup);
      return false;
    }

  if (k.size() == 1)
    {
      hid_t    dataset, dataspace, dset_create_props;

      hsize_t rank = (k[0].depth() > 1) ? 3 : ((k[0].height() > 1 && k[0].depth() == 1) ? 2 : ((k[0].width() > 1 && k[0].height() == 1 && k[0].depth() == 1) ? 1 : 0));
      rank += 1; //for internal dimension

      hsize_t *dimsf = (hsize_t *) malloc(rank*sizeof(hsize_t));

      if(rank == 4)
        {
          dimsf[0] = k[0].int_dim();
          dimsf[1] = k[0].width();
          dimsf[2] = k[0].height();
          dimsf[3] = 2*k[0].depth();
        }
      else if(rank == 3)
        {
          dimsf[0] = k[0].int_dim();
          dimsf[1] = k[0].width();
          dimsf[2] = 2*k[0].height();
        }
      else if(rank == 2)
        {
          dimsf[0] = k[0].int_dim();
          dimsf[1] = 2*k[0].width();
        }
      else
        dimsf[0] = 2*k[0].int_dim();

      dataspace = H5Screate_simple(rank, dimsf, NULL);

      dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list

      dataset = H5Dcreate(h5_CGgroup, name.c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT); //create data block in file

      H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)k[0].get_address()); //write grid data

      //clean up
      free(dimsf);
      H5Dclose(dataset);
      H5Pclose(dset_create_props);
      H5Sclose(dataspace);
    }
  else
    {

      hid_t h5_CGgroup_vecsub = H5Gcreate(h5_CGgroup, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); //create dataset with full space selection

      for (int i = 0; i < k.size(); i++)
        {
          hid_t    dataset, dataspace, dset_create_props, h5a_time;

          hsize_t rank = (k[i].depth() > 1) ? 3 : ((k[i].height() > 1 && k[i].depth() == 1) ? 2 : ((k[i].width() > 1 && k[i].height() == 1 && k[i].depth() == 1) ? 1 : 0));
          rank += 1; //for internal dimension

          hsize_t *dimsf = (hsize_t *) malloc(rank*sizeof(hsize_t));

          if(rank == 4)
            {
              dimsf[0] = k[i].int_dim();
              dimsf[1] = k[i].width();
              dimsf[2] = k[i].height();
              dimsf[3] = 2*k[i].depth();
            }
          else if(rank == 3)
            {
              dimsf[0] = k[i].int_dim();
              dimsf[1] = k[i].width();
              dimsf[2] = 2*k[i].height();
            }
          else if(rank == 2)
            {
              dimsf[0] = k[i].int_dim();
              dimsf[1] = 2*k[i].width();
            }
          else
            dimsf[0] = 2*k[0].int_dim();

          dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list
          dataspace = H5Screate_simple(rank, dimsf, NULL); //define space in file

          stringstream comp;
          comp << i;

          dataset = H5Dcreate(h5_CGgroup_vecsub, (comp.str()).c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT); //create data block in file

          H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)k[i].get_address()); //write grid data

          free(dimsf);
          H5Dclose(dataset);
          H5Pclose(dset_create_props);
          H5Sclose(dataspace);
        }

      if(options.delta_t.size())
        {
          //and a the following time intervalls
          hsize_t dimsf2[] = {options.delta_t.size()};
          hsize_t dataspace = H5Screate_simple(1, dimsf2, NULL);
          hsize_t h5a_time = H5Acreate(h5_CGgroup_vecsub, "delta_t", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
          H5Awrite(h5a_time, H5T_IEEE_F64LE, &options.delta_t.front());

          H5Aclose(h5a_time);
          H5Sclose(dataspace);
        }

      H5Gclose(h5_CGgroup_vecsub);
    }

  H5Gclose(h5_CGgroup);
  H5Gclose(h5_timegroup);

  return true;
}

bool Bh3BinaryFile::append_snapshot(const string &name, double time, const vector<RealGrid> &k)
{
  if(m == in)
    {
      cout << "file "<< filename.c_str() << "is not in write mode" << endl;
      return false;
    }

  if(!check_time(time))
    {
      cout << "HDF5 Group error for time group: " << time << endl;
      H5Gclose(h5_timegroup);

      return false;
    }

  hid_t h5_RGgroup;

  if(H5Lexists(h5_timegroup, "RealGrid", H5P_DEFAULT))
    h5_RGgroup = H5Gopen(h5_timegroup, "RealGrid", H5P_DEFAULT);
  else
    h5_RGgroup = H5Gcreate(h5_timegroup, "RealGrid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


  if(H5Lexists(h5_RGgroup, name.c_str(), H5P_DEFAULT))
    {
      cout << "RealGrid " << name << " already exists for time " << time << ". I refuse to write." << endl;
      H5Gclose(h5_RGgroup);
      H5Gclose(h5_timegroup);

      return false;
    }

  if (k.size() == 1)
    {
      hid_t    dataset, dataspace, dset_create_props;

      hsize_t rank = (k[0].fft_depth() > 1) ? 3 : ((k[0].fft_height() > 1 && k[0].fft_depth() == 1) ? 2 : ((k[0].fft_width() > 1 && k[0].fft_height() == 1 && k[0].fft_depth() == 1) ? 1 : 0));
      rank += 1; //for internal dimension

      hsize_t *dimsf = (hsize_t *) malloc(rank*sizeof(hsize_t));

      if(rank == 4)
        {
          dimsf[0] = k[0].int_dim();
          dimsf[1] = k[0].fft_width();
          dimsf[2] = k[0].fft_height();
          dimsf[3] = k[0].fft_depth();
        }
      else if(rank == 3)
        {
          dimsf[0] = k[0].int_dim();
          dimsf[1] = k[0].fft_width();
          dimsf[2] = k[0].fft_height();
        }
      else if(rank == 2)
        {
          dimsf[0] = k[0].int_dim();
          dimsf[1] = k[0].fft_width();
        }
      else
        dimsf[0] = 2*k[0].int_dim();

      dataspace = H5Screate_simple(rank, dimsf, NULL);

      dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list

      dataset = H5Dcreate(h5_RGgroup, name.c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT); //create data block in file

      H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, k[0].get_address()); //write grid data

      //clean up
      free(dimsf);
      H5Dclose(dataset);
      H5Pclose(dset_create_props);
      H5Sclose(dataspace);
    }
  else
    {
      hid_t h5_RGgroup_vecsub = H5Gcreate(h5_RGgroup, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); //create dataset with full space selection

      for (int i = 0; i < k.size(); i++)
        {
          hid_t    dataset, dataspace, dset_create_props, h5a_time;

          hsize_t rank = (k[i].fft_depth() > 1) ? 3 : ((k[i].fft_height() > 1 && k[i].fft_depth() == 1) ? 2 : ((k[i].fft_width() > 1 && k[i].fft_height() == 1 && k[i].fft_depth() == 1) ? 1 : 0));
          rank += 1; //for internal dimension

          hsize_t *dimsf = (hsize_t *) malloc(rank*sizeof(hsize_t));

          if(rank == 4)
            {
              dimsf[0] = k[i].int_dim();
              dimsf[1] = k[i].fft_width();
              dimsf[2] = k[i].fft_height();
              dimsf[3] = k[i].fft_depth();
            }
          else if(rank == 3)
            {
              dimsf[0] = k[i].int_dim();
              dimsf[1] = k[i].fft_width();
              dimsf[2] = k[i].fft_height();
            }
          else if(rank == 2)
            {
              dimsf[0] = k[i].int_dim();
              dimsf[1] = k[i].fft_width();
            }
          else
            dimsf[0] = k[0].int_dim();

          dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list
          dataspace = H5Screate_simple(rank, dimsf, NULL); //define space in file

          stringstream comp;
          comp << i;

          dataset = H5Dcreate(h5_RGgroup_vecsub, (comp.str()).c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT); //create data block in file

          H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)k[i].get_address()); //write grid data

          free(dimsf);
          H5Dclose(dataset);
          H5Pclose(dset_create_props);
          H5Sclose(dataspace);
        }

      //and a the following time intervalls
      if(options.delta_t.size() > 0)
        {
          hsize_t dimsf2[] = {options.delta_t.size()};
          hsize_t dataspace = H5Screate_simple(1, dimsf2, NULL);
          hsize_t h5a_time = H5Acreate(h5_RGgroup_vecsub, "delta_t", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
          H5Awrite(h5a_time, H5T_IEEE_F64LE, &options.delta_t.front());

          H5Aclose(h5a_time);
          H5Sclose(dataspace);
        }

      H5Gclose(h5_RGgroup_vecsub);
    }

  H5Gclose(h5_RGgroup);
  H5Gclose(h5_timegroup);

  return true;
}

bool Bh3BinaryFile::get_snapshot(const string &name, double time, vector<ComplexGrid> &k)
{
  if(m == out)
    {
      cout << "file "<< filename.c_str() << "is not in read mode" << endl;
      return false;
    }

  stringstream time_name;
  // time will always be formatted with two decimal digits
  time_name.setf(ios_base::fixed);
  time_name.precision(2);
  time_name << "/" << time;

  if(H5Lexists(h5_file, (time_name.str()).c_str(), H5P_DEFAULT))
    {
      h5_timegroup = H5Gopen(h5_file, (time_name.str()).c_str(), H5P_DEFAULT);

    }
  else
    {
      cout << "ERROR: HDF5 Group for time " << time << " does not exist" << endl;
      return false;
    }

  stringstream set_name;
  set_name << "ComplexGrid/" << name;

  if(!H5Lexists(h5_timegroup, (set_name.str()).c_str(), H5P_DEFAULT))
    {
      cout << "ERROR: HDF5 Dataset " << name << " does not exist for time " << time << " in ComplexGrid branch" << endl;

      H5Gclose(h5_timegroup);
      return false;
    }
  else
    {
      hid_t test_id = H5Oopen(h5_timegroup, (set_name.str()).c_str(), H5P_DEFAULT);

      if(H5Iget_type(test_id) == H5I_DATASET)
        {
          hid_t dataspace = H5Dget_space(test_id);

          hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
          hsize_t *dimf = (hsize_t *)malloc(rank*sizeof(hsize_t));
          H5Sget_simple_extent_dims(dataspace, dimf, NULL);

          if(rank == 4)
            {
              dimf[3] /= 2;
              k.assign(1, ComplexGrid(dimf[0], dimf[1], dimf[2], dimf[3]));
            }
          else if(rank == 3)
            {
              dimf[2] /= 2;
              k.assign(1, ComplexGrid(dimf[0], dimf[1], dimf[2], 1));
            }
          else if(rank == 2)
            {
              dimf[1] /= 2;
              k.assign(1, ComplexGrid(dimf[0], dimf[1], 1, 1));
            }
          else
            {
              dimf[0] /= 2;
              k.assign(1, ComplexGrid(dimf[0], 1 , 1, 1));
            }

          H5Dread(test_id,  H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)k[0].get_address());

          free(dimf);
          H5Dclose(test_id);
        }
      else
        {
          hsize_t vecsize;
          H5Gget_num_objs(test_id, &vecsize);
          k.resize(vecsize);

          for(int i = 0; i < vecsize; i++)
            {
              stringstream comp;
              comp << i;

              hid_t dataset = H5Dopen(test_id, (comp.str()).c_str(), H5P_DEFAULT);
              hid_t dataspace = H5Dget_space(dataset);

              hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
              hsize_t *dimf = (hsize_t *)malloc(rank*sizeof(hsize_t));
              H5Sget_simple_extent_dims(dataspace, dimf, NULL);

              if(rank == 4)
                {
                  dimf[3] /= 2;
                  ComplexGrid g(dimf[0], dimf[1], dimf[2], dimf[3]);
                  k[i] = g;
                }
              else if(rank == 3)
                {
                  dimf[2] /= 2;
                  ComplexGrid g(dimf[0], dimf[1], dimf[2], 1);
                  k[i] = g;
                }
              else if(rank == 2)
                {
                  dimf[1] /= 2;
                  ComplexGrid g(dimf[0], dimf[1], 1, 1);
                  k[i] = g;
                }
              else
                {
                  dimf[0] /= 2;
                  ComplexGrid g(dimf[0], 1 , 1, 1);
                  k[i] = g;
                }

              H5Dread(dataset,  H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)k[i].get_address());

              free(dimf);
              H5Dclose(dataset);
            }
          H5Gclose(test_id);
        }
    }

  H5Gclose(h5_timegroup);

  return true;
}


bool Bh3BinaryFile::get_snapshot(const string &name, double time, vector<RealGrid> &k)
{
  if(m == out)
    {
      cout << "file "<< filename.c_str() << "is not in read mode" << endl;
      return false;
    }

  stringstream time_name;
  // time will always be formatted with two decimal digits
  time_name.setf(ios_base::fixed);
  time_name.precision(2);
  time_name << "/" << time;

  if(H5Lexists(h5_file, (time_name.str()).c_str(), H5P_DEFAULT))
    {
      h5_timegroup = H5Gopen(h5_file, (time_name.str()).c_str(), H5P_DEFAULT);

    }
  else
    {
      cout << "ERROR: HDF5 Group for time " << time << " does not exist" << endl;
      return false;
    }

  stringstream set_name;
  set_name << "RealGrid/" << name;

  if(!H5Lexists(h5_timegroup, (set_name.str()).c_str(), H5P_DEFAULT))
    {
      cout << "ERROR: HDF5 Dataset " << name << " does not exist for time " << time << " in RealGrid branch" << endl;

      H5Gclose(h5_timegroup);
      return false;
    }
  else
    {
      hid_t test_id = H5Oopen(h5_timegroup, (set_name.str()).c_str(), H5P_DEFAULT);

      if(H5Iget_type(test_id) == H5I_DATASET)
        {
          hid_t dataspace = H5Dget_space(test_id);

          hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
          hsize_t *dimf = (hsize_t *)malloc(rank*sizeof(hsize_t));
          H5Sget_simple_extent_dims(dataspace, dimf, NULL);

          if(rank == 4)
            {
              dimf[3] = 2*((int)(dimf[3]/2) -1);
              k.assign(1, RealGrid(dimf[0], dimf[1], dimf[2], dimf[3]));
            }
          else if(rank == 3)
            {
              dimf[2] = 2*((int)(dimf[2]/2) -1);
              k.assign(1, RealGrid(dimf[0], dimf[1], dimf[2], 1));
            }
          else if(rank == 2)
            {
              dimf[1] = 2*((int)(dimf[1]/2) -1);
              k.assign(1, RealGrid(dimf[0], dimf[1], 1, 1));
            }
          else
            {
              dimf[0] = 2*((int)(dimf[0]/2) -1);
              k.assign(1, RealGrid(dimf[0], 1 , 1, 1));
            }

          H5Dread(test_id,  H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)k[0].get_address());

          free(dimf);
          H5Dclose(test_id);
        }
      else
        {
          hsize_t vecsize;
          H5Gget_num_objs(test_id, &vecsize);
          k.resize(vecsize);

          for(int i = 0; i < vecsize; i++)
            {
              stringstream comp;
              comp << i;

              hid_t dataset = H5Dopen(test_id, (comp.str()).c_str(), H5P_DEFAULT);
              hid_t dataspace = H5Dget_space(dataset);

              hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
              hsize_t *dimf = (hsize_t *)malloc(rank*sizeof(hsize_t));
              H5Sget_simple_extent_dims(dataspace, dimf, NULL);

              if(rank == 4)
                {
                  dimf[3] = 2*((int)(dimf[3]/2) -1);
                  RealGrid g(dimf[0], dimf[1], dimf[2], dimf[3]);
                  k[i] = g;
                }
              else if(rank == 3)
                {
                  dimf[2] = 2*((int)(dimf[2]/2) -1);
                  RealGrid g(dimf[0], dimf[1], dimf[2], 1);
                  k[i] = g;
                }
              else if(rank == 2)
                {
                  dimf[1] = 2*((int)(dimf[1]/2) -1);
                  RealGrid g(dimf[0], dimf[1], 1, 1);
                  k[i] = g;
                }
              else
                {
                  dimf[0] = 2*((int)(dimf[0]/2) -1);
                  RealGrid g(dimf[0], 1 , 1, 1);
                  k[i] = g;
                }

              H5Dread(dataset,  H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)k[i].get_address());

              free(dimf);
              H5Dclose(dataset);
            }
          H5Gclose(test_id);
        }
    }

  H5Gclose(h5_timegroup);

  return true;
}

bool Bh3BinaryFile::append_averages(const string &vec_name, double *vec, int vec_rank, int* vec_dim, double time)
{
  if(m == in)
    {
      cout << "file "<< filename.c_str() << "is not in write mode" << endl;
      return false;
    }

  if(!check_time(time))
    {
      cout << "HDF5 Group error for time group: " << time << endl;
      return false;
    }

  hid_t h5_avgroup;

  if(H5Lexists(h5_timegroup, "Averages", H5P_DEFAULT))
    h5_avgroup = H5Gopen(h5_timegroup, "Averages", H5P_DEFAULT);
  else
    h5_avgroup = H5Gcreate(h5_timegroup, "Averages", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


  if(H5Lexists(h5_avgroup, vec_name.c_str(), H5P_DEFAULT))
    {
      cout << "Averages " << vec_name << " already exists for time " << time << ". I refuse to write." << endl;
      return false;
    }



  hid_t dataset, dataspace, h5a_time, dset_create_props;

  dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list

  hsize_t *dimsf = (hsize_t*)malloc(vec_rank*sizeof(hsize_t));
  for(int l = 0; l < vec_rank; l++)
    {
      dimsf[l] = vec_dim[l];
    }

  dataspace = H5Screate_simple(1, dimsf,  NULL);

  dataset = H5Dcreate(h5_avgroup, vec_name.c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
  H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec); //write grid data

  //append grid time as attribute

  free(dimsf);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Pclose(dset_create_props);

  H5Gclose(h5_avgroup);
  H5Gclose(h5_timegroup);

  return true;
}

bool Bh3BinaryFile::append_docstring(const string &group, const string &docstring, double time)
{
  if(m == in)
    {
      cout << "file "<< filename.c_str() << "is not in write mode" << endl;
      return false;
    }

  if(!check_time(time))
    {
      cout << "HDF5 Group error for time group: " << time << endl;
      return false;
    }

  hid_t h5_avgroup;

  if(H5Lexists(h5_timegroup, group.c_str(), H5P_DEFAULT))
    h5_avgroup = H5Gopen(h5_timegroup, group.c_str(), H5P_DEFAULT);
  else
    h5_avgroup = H5Gcreate(h5_timegroup, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t dataset, dataspace, dset_create_props;

  dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list

  hsize_t dimsf[] = {1};

  dataspace = H5Screate_simple(1, dimsf,  NULL);
  hid_t atype = H5Tcopy(H5T_C_S1);
  hsize_t size = docstring.size();
  H5Tset_size(atype, size);

  dataset = H5Dcreate(h5_avgroup, "Docstring", atype, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
  H5Dwrite (dataset, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(docstring[0])); //write docstring

  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Pclose(dset_create_props);

  H5Gclose(h5_avgroup);
  H5Gclose(h5_timegroup);

  return true;
}



bool initialize_binary_dir(string &dname, const PathOptions &opt)
{
  stringstream dirname;
  dirname << dname << "_" << "TS" << opt.timestepsize
          << "_C" << opt.grid[0]
          << "_G" << opt.grid[1] << "_" << opt.grid[2] << "_" << opt.grid[3]
          << "_N" << opt.N
          << "_U" << opt.U
          << "_KL" << opt.klength[0] << "_" << opt.klength[1] << "_" << opt.klength[2];
  // Neuen Dateinamen garantieren

  DIR *dir = opendir(".");

  if(dir==NULL)
    {
      cout << "Couldn't open current directory!\n";
      return false;
    }

  dirent *entry;
  int run=0;
  char format[512];
  strcpy(format,dirname.str().c_str());
  strcat(format,"Run%d");
  while ( (entry=readdir(dir)) != NULL )
    {
      int temp;
      int result = sscanf(entry->d_name, format, &temp);
      if((result == 1) && (temp > run))
        run = temp;
    }
  run++;
  dirname << "Run" << run;
  closedir(dir);

  if(mkdir(dirname.str().c_str(),0755) != 0)
    {
      cout << "Couldn't create directory " << dirname.str() << " !\n";
      return false;
    }
  dirname << "/";
  ofstream description;
  try
    {
      description.open((dirname.str() + "Bh3Description.txt").c_str());
    }
  catch (exception &e)
    {
      cout << "Couldn't create one of the files!" << endl;
      return false;
    }

  stringstream ds;
  ds              << "Name			: Bh3-512" << endl
                  << "Time-StepSize		: " << opt.timestepsize << endl
                  << "Gridsize			: " << opt.grid[1] << "," << opt.grid[2] << "," << opt.grid[3] << endl
                  << "Number of particles		: " << opt.N << endl
                  << "Interaction strength (U)	: " << opt.U << endl
                  << "Comments			: Production run: vortex counting";
  description << ds.str();
  description.close();
  dname = dirname.str();
  return true;
}

bool initialize_path_options(int argc, char** argv, PathOptions &opt, vector<double> &snapshot_times, void *gen_times(vector<double> &snapshot_times), bool verbose)
{
  config_t conf;
  config_setting_t *setting;
  config_init(&conf);

  string config_file;

  if (argc > 1)
    config_file = argv[1];
  else
    config_file = "/home/karl/myPyModules/bh3default.conf";


  if (config_read_file(&conf, config_file.c_str()) == CONFIG_FALSE)
    {
      cout << "error " << config_error_text(&conf) << " occurred during read process" << endl;
      config_destroy(&conf);
      return false;
    }

  setting = config_lookup(&conf, "PathOptions");

  if (setting == NULL)
    {
      cout << "error: Setting \"PathOptions\" not found in " << config_setting_source_file(setting) << endl;
      config_destroy(&conf);
      return false;
    }

  if (config_setting_lookup_float(setting, "timestepsize", &opt.timestepsize) == CONFIG_FALSE)
    cout << "setting \"timestepsize\" not found" << endl;

  if (config_setting_lookup_float(setting, "U", &opt.U)== CONFIG_FALSE)
    cout << "setting \"U\" not found" << endl ;

  if (config_setting_lookup_float(setting, "N", &opt.N)== CONFIG_FALSE)
    cout << "setting \"N\" not found" << endl ;

  if (config_setting_lookup_int(setting, "internal_dimension", &opt.grid[0])== CONFIG_FALSE)
    cout << "setting \"internal_dimension\" not found" << endl ;

  setting = config_lookup(&conf, "PathOptions.grid_dims");
  if (setting == NULL)
    {
      cout << "error: Setting \"PathOptions.grid_dims\" not found in " << config_setting_source_file(setting)<< endl;
      config_destroy(&conf);
      return false;
    }
  for (int i = 0; i < 3; i++)
    {
      opt.grid[i+1] = config_setting_get_int_elem (setting, i);
    }

  setting = config_lookup(&conf, "PathOptions.klengths");
  if (setting == NULL)
    {
      cout << "error: Setting \"PathOptions.klength\" not found in " << config_setting_source_file(setting)<< endl;
      config_destroy(&conf);
      return false;
    }
  for (int i = 0; i < 3; i++)
    opt.klength[i] = config_setting_get_float_elem (setting, i);


  setting = config_lookup(&conf, "PathOptions.gmatrix");
  if (setting == NULL)
    {
      cout << "error: Setting \"PathOptions.gmatrix\" not found in " << config_setting_source_file(setting)<< endl;
      config_destroy(&conf);
      return false;
    }
  int index = config_setting_length(setting);
  if(index != opt.grid[0]*opt.grid[0])
    {
      cout << "###############################################################################"<<endl;
      cout << "# WARNING: dimension of gmatrix in file does not match the internal dimension #"<<endl;
      cout << "###############################################################################"<<endl;
    }

  opt.g.resize(index);
  for (int i = 0; i<index; i++)
    opt.g[i] = config_setting_get_float_elem (setting, i);

  setting = config_lookup(&conf, "PathOptions.delta_t");
  if (setting != NULL)
    {
      index = config_setting_length(setting);
      opt.delta_t.resize(index);
      for (int i = 0; i<index; i++)
        opt.delta_t[i] = config_setting_get_float_elem (setting, i);
    }
  else
    opt.delta_t.resize(0);

  //print PathOptions if desired
  if(verbose)
    cout << opt << endl;

  //generate Snapshots either from given function gen_times
  if (gen_times != NULL)
    gen_times(snapshot_times);
  //or from config file
  else
    {
      setting = config_lookup(&conf, "SnapshotTimes");
      if (setting != NULL)
        {
          index = config_setting_length(setting);
          snapshot_times.resize(index);
          for (int i = 0; i<index; i++)
            {
              snapshot_times[i] = config_setting_get_float_elem (setting, i);
            }
        }
      //if no snapshots are provided, neither from config file nor from function, make snapshot from initial field
      else
        {
          snapshot_times.resize(1);
          snapshot_times[0] = 0.;
        }
    }

  config_destroy(&conf);
  return true;

}


list<string> get_file_list(const list<string> args)
{
  list<string> files;
  PathOptions options;

  // iterate the arguments and check them for valid Bh3-files
  for(list<string>::const_iterator bin = args.begin(); bin != args.end(); ++bin)
    {
      struct stat buffer;
      lstat((*bin).c_str(), &buffer);
      if(S_ISREG(buffer.st_mode))
        {
          Bh3BinaryFile b(*bin, options, Bh3BinaryFile::in);
          if(b.get_timelist().size() > 0)
            files.push_back(*bin);
        }
      else
        {
          string bindir = *bin;
          if(bindir[bindir.length() - 1] != '/')
            bindir = bindir + "/";

          // Dateiname der Binaerdateien fuer Analyse bestimmen
          DIR *dir = opendir(bindir.c_str());
          if (dir == NULL)
            {
              cout << "Couldn't open directory with snapshots!" << endl;
              return list<string>();
            }
          dirent *entry;
          while( (entry = readdir(dir)) != NULL)
            {
              string filename = bindir + entry->d_name;
              Bh3BinaryFile b(filename, options, Bh3BinaryFile::in);
              if(b.get_timelist().size() > 0)
                files.push_back(filename);
            }
          closedir(dir);
        }
    }
  return files;
}
