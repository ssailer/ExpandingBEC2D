#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <unistd.h>
#include <dirent.h>
#include <string.h>

#include <EXP2D_binaryfile.h>
// #include <complexgrid.h>
// #include <realgrid.h>
#include <hdf5.h>
#include <eigen3/Eigen/Dense>
// #include <libconfig.h>

// Class implementation
binaryFile::binaryFile(const string &file, mode nm)
{
  filename = file;
  // options = opt;
  m = nm;

  if(m == in || m == append)
    {
      struct stat buf;
      lstat(filename.c_str(), &buf);
      if(S_ISREG(buf.st_mode) && H5Fis_hdf5(filename.c_str())) // Nur normale Dateien ueberpruefen
        { 
          // cerr << "Reached ERROR location #1" << endl;
          if(m == in)
            h5_file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
          else if(m == append)
            h5_file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

          if(H5Aexists(h5_file, "SnapshotTimes") > 0)
            {
              // cerr << "Reached ERROR location #2" << endl;

              hid_t h5a_options;
              // h5a_options = H5Aopen(h5_file, "Options", H5P_DEFAULT);

              // double tmpOpt1[24];

              // H5Aread(h5a_options, H5T_IEEE_F64LE , tmpOpt1);

              // //load Options struct from file array. Don't forget to change appropriately when changing Options struct
              // options.N = tmpOpt1[0];
              // for(int i= 0; i < 3; i++){
              // 	options.klength[i] = tmpOpt1[i+1];
              // }
              // options.stateInformation[0] = tmpOpt1[4];
              // options.stateInformation[1] = tmpOpt1[5];
              // options.omega_x = complex<double>(tmpOpt1[6],0.0);
              // options.omega_y = complex<double>(tmpOpt1[7],0.0);
              // options.dispersion_x = complex<double>(tmpOpt1[8],0.0);
              // options.dispersion_x = complex<double>(tmpOpt1[9],0.0);
              // options.min_x = tmpOpt1[10];
              // options.min_y = tmpOpt1[11];
              // options.t_abs = complex<double>(tmpOpt1[12],0.0);
              // options.exp_factor = complex<double>(tmpOpt1[13],0.0);
              // options.g = tmpOpt1[14];
              // options.ITP_step = tmpOpt1[15];
              // options.RTE_step = tmpOpt1[16];

              // for(int i = 0; i<4;i++)
              //   options.grid[i] = (uint32_t)tmpOpt1[i+17];

              // options.n_it_RTE = (int)tmpOpt1[21];
              // options.samplesize = (int)tmpOpt1[22];
              // options.vortexnumber = (int)tmpOpt1[23];              

              // H5Aclose(h5a_options);

              // h5a_options = H5Aopen(h5_file, "strings", H5P_DEFAULT);




              // string tmpOpt2[5];

              // H5Aread(h5a_options,H5T_C_S1, tmpOpt2);
              // 	options.runmode = tmpOpt2[0];
              // 	options.name = tmpOpt2[1];
              // 	options.config = tmpOpt2[2];
              // 	options.workingdirectory = tmpOpt2[3];
              // 	options.workingfile = tmpOpt2[4];

              // H5Aclose(h5a_options);

              // //load g_matrix from file
              // if(H5Aexists(h5_file, "g_matrix") > 0)
              //   {
              //     h5a_options = H5Aopen(h5_file, "g_matrix", H5P_DEFAULT);
              //     int size = (int)(H5Aget_storage_size(h5a_options)/sizeof(double));
              //     options.g.resize(size);

              //     H5Aread(h5a_options, H5T_IEEE_F64LE , &(options.g.front()));

              //     H5Aclose(h5a_options);
              //   }
              // else
              //   options.g.resize(0);

              // //load delta_t vector from file, if exists
              // if(H5Aexists(h5_file, "delta_t") > 0)
              //   {
              //     h5a_options = H5Aopen(h5_file, "delta_t", H5P_DEFAULT);
              //     int size = (int)(H5Aget_storage_size(h5a_options)/sizeof(double));
              //     options.delta_t.resize(size);

              //     H5Aread(h5a_options, H5T_IEEE_F64LE , &(options.delta_t.front()));

              //     H5Aclose(h5a_options);
              //   }
              // else
              //   options.delta_t.resize(0);

              h5a_options = H5Aopen(h5_file, "SnapshotTimes", H5P_DEFAULT);
              int size = (int)(H5Aget_storage_size(h5a_options)/sizeof(double));
              time_list.resize(size);
              H5Aread(h5a_options, H5T_IEEE_F64LE , &time_list.front());
              H5Aclose(h5a_options);

              // if(m == append)
              //   {
              //     if(options != opt)
              //       cout << "WARNING for I/O operation on file: " << filename <<". Passed options do not fit file options in 'append' mode." << endl;
              //   }
              // cerr << "Reached ERROR location #3" << endl;

            }
          else
            {
              cout << "WARNING for I/O operation on file: "<< filename << " occurred: No valid SnapshotTimes attribute found in file." << endl;
            }
        }
      else
        {
          cout << "ERROR for I/O operation. File: "<< filename <<" is either not regular or not in hdf5 storage format" << endl;
        }
    }
  else
    {
        // cerr << "Reached ERROR location #4" << endl;

      h5_file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      // hid_t    h5a_options, dataspace;


      // //load Options struct to a reliable arrays. Don't forget to change appropriately when changing Options struct
      // double tmpOpt1[24];

      // tmpOpt1[0] = options.N;
      // for(int i= 0; i < 3; i++)
      // 	tmpOpt1[i+1] = options.klength[i];

      // tmpOpt1[4] = options.stateInformation[0];
      // tmpOpt1[5] = options.stateInformation[1];
      // tmpOpt1[6] = options.omega_x.real();
      // tmpOpt1[7] = options.omega_y.real();
      // tmpOpt1[8] = options.dispersion_x.real();
      // tmpOpt1[9] = options.dispersion_y.real();
      // tmpOpt1[10] = options.min_x;
      // tmpOpt1[11] = options.min_y;
      // tmpOpt1[12] = options.t_abs.real();
      // tmpOpt1[13] = options.exp_factor.real();
      // tmpOpt1[14] = options.g;
      // tmpOpt1[15] = options.ITP_step;
      // tmpOpt1[16] = options.RTE_step;

      // for(int i = 0; i<4;i++)
      //   tmpOpt1[i+17] = options.grid[i];

      // tmpOpt1[21] = options.n_it_RTE;
      // tmpOpt1[22] = options.samplesize;
      // tmpOpt1[23] = options.vortexnumber;  

      // //copy options to file
      // hsize_t dimsf[1] = {24};
      // dataspace = H5Screate_simple(1, dimsf, NULL);

      // h5a_options = H5Acreate(h5_file, "Options", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      // H5Awrite (h5a_options, H5T_IEEE_F64LE, tmpOpt1);

      // H5Aclose(h5a_options);
      // H5Sclose(dataspace);

      // string tmpOpt2[5];

      // tmpOpt2[0] = options.runmode;
      // tmpOpt2[1] = options.name;
      // tmpOpt2[2] = options.config;
      // tmpOpt2[3] = options.workingdirectory;
      // tmpOpt2[4] = options.workingfile;

      // dimsf[0] = {5};
      // dataspace = H5Acreate(h5_file, "strings", H5T_C_S1, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      // H5Awrite(h5a_options, H5T_C_S1, tmpOpt2);
      // H5Aclose(h5a_options);
      // H5Sclose(dataspace);

      // //copy delta_t array to file if there is one
      // if (options.delta_t.size()>0)
      //   {
      //     dimsf[0] = options.delta_t.size();
      //     dataspace = H5Screate_simple(1, dimsf, NULL);

      //     h5a_options = H5Acreate(h5_file, "delta_t", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      //     H5Awrite (h5a_options, H5T_IEEE_F64LE, &(options.delta_t.front()));

      //     H5Aclose(h5a_options);
      //     H5Sclose(dataspace);
      //   }

      // //copy g matrix to file

      // if(options.g.size()>0)
      //   {
      //     dimsf[0] = options.g.size();
      //     dataspace = H5Screate_simple(1, dimsf, NULL);

      //     h5a_options = H5Acreate(h5_file, "g_matrix", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      //     H5Awrite (h5a_options, H5T_IEEE_F64LE, &(options.g.front()));

      //     H5Aclose(h5a_options);
      //     H5Sclose(dataspace);
      //   }

    }
}

binaryFile::~binaryFile()
{
  close();
}

void binaryFile::close()
{
  if(time_list.size()>0 && m==out)
    {      
      hsize_t dimsf[] = {time_list.size()};
      hid_t dataspace = H5Screate_simple(1, dimsf, NULL);
      hid_t h5a = H5Acreate(h5_file, "SnapshotTimes", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite (h5a, H5T_IEEE_F64LE, &time_list.front());
      H5Sclose(dataspace);
      H5Aclose(h5a);
      // cerr << "Reached ERROR location #5" << endl;
    }
  else if(time_list.size()>0 && m == append)
    {
      
      if(H5Aexists(h5_file, "SnapshotTimes"))
        {
          H5Adelete(h5_file, "SnapshotTimes");
        }

      hsize_t dimsf[] = {time_list.size()};
      hid_t dataspace = H5Screate_simple(1, dimsf, NULL);
      hid_t h5a = H5Acreate(h5_file, "SnapshotTimes", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite (h5a, H5T_IEEE_F64LE, &time_list.front());
      H5Sclose(dataspace);
      H5Aclose(h5a);
      // cerr << "Reached ERROR location #6" << endl;
    }


  H5Fclose(h5_file);
  filename = "";
}

bool binaryFile::checkTime(double snapShotTime)
{
  stringstream time_name;
  // snapShotTime will always be formatted with two decimal digits
  time_name.setf(ios_base::fixed);
  time_name.precision(2);
  time_name << snapShotTime;

  if(H5Lexists(h5_file, (time_name.str()).c_str(), H5P_DEFAULT)){
    h5_timegroup = H5Gopen(h5_file, (time_name.str()).c_str(), H5P_DEFAULT);
    // cerr << "Reached ERROR location #7" << endl;
  }
  else
    {
      h5_timegroup = H5Gcreate(h5_file, (time_name.str()).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      time_list.push_back(snapShotTime);
      // cerr << "Reached ERROR location #8" << endl;
    }

  if(h5_timegroup >= 0)
    return true;
  else
    return false;
}


bool binaryFile::appendSnapshot(const string &name, double snapShotTime, const vector<MatrixXcd> &k, Options &options)
{
  if(m == in)
    {
      cout << "file "<< filename.c_str() << "is not in write mode" << endl;
      return false;
    }

  if(!checkTime(snapShotTime))
    {
      cout << "HDF5 Group error for snapShotTime group: " << snapShotTime << endl;
      H5Gclose(h5_timegroup);
      return false;
    }

    // cerr << "Reached ERROR location #9" << endl;

  // hid_t h5_EMGroup;

  // if(H5Lexists(h5_timegroup, "EigenMatrix", H5P_DEFAULT))
  //   h5_EMGroup = H5Gopen(h5_timegroup, "EigenMatrix", H5P_DEFAULT);
  // else
  //   h5_EMGroup = H5Gcreate(h5_timegroup, "EigenMatrix", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // // cerr << "Reached ERROR location #10" << endl;
  // if(H5Lexists(h5_EMGroup, name.c_str(), H5P_DEFAULT))
  //   {
  //     cout << "EigenMatrix " << name << " already exists for snapShotTime " << snapShotTime << ". I refuse to write." << endl;
  //     H5Gclose(h5_EMGroup);
  //     H5Gclose(h5_timegroup);
  //     return false;
  //   }

  // if (k.size() == 1)
  //   {
  //     hid_t    dataset, dataspace, dset_create_props;

  //     hsize_t rank = 2;
  //     // rank += 1; //for internal dimension

  //     hsize_t *dimsf = (hsize_t *) malloc(rank*sizeof(hsize_t));

  //     dimsf[0] = 2*k[0].rows();
  //     dimsf[1] = 2*k[0].cols();


  //     dataspace = H5Screate_simple(rank, dimsf, NULL);

  //     dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list

  //     dataset = H5Dcreate(h5_EMGroup, name.c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT); //create data block in file

  //     H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)k[0].data()); //write grid data

  //     //clean up
  //     free(dimsf);
  //     H5Dclose(dataset);
  //     H5Pclose(dset_create_props);
  //     H5Sclose(dataspace);
  //   }
  // else
  //   {

      hid_t h5_EMGroup_vecsub = H5Gcreate(h5_timegroup, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); //create dataset with full space selection
      // cerr << "Reached ERROR location #11" << endl;
      for (int i = 0; i < k.size(); i++)
        {
          hid_t    dataset, dataspace, dset_create_props;

          hsize_t rank = 2;
          // rank += 1; //for internal dimension

          hsize_t *dimsf = (hsize_t *) malloc(rank*sizeof(hsize_t));

          dimsf[0] = 2*k[i].rows();
          dimsf[1] = k[i].cols();


          dset_create_props = H5Pcreate(H5P_DATASET_CREATE); //create a default creation property list
          dataspace = H5Screate_simple(rank, dimsf, NULL); //define space in file

          stringstream comp;
          comp << i;

          dataset = H5Dcreate(h5_EMGroup_vecsub, (comp.str()).c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT); //create data block in file

          H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)k[i].data()); //write grid data

          free(dimsf);
          H5Dclose(dataset);
          H5Pclose(dset_create_props);
          H5Sclose(dataspace);
        }
        // cerr << "Reached ERROR location #12" << endl;




      H5Gclose(h5_EMGroup_vecsub);
    // }
  // H5Gclose(h5_EMGroup);

      hid_t    h5a_options, dataspace;
      double tmpOpt1[24];

      tmpOpt1[0] = options.N;
      for(int i= 0; i < 3; i++)
        tmpOpt1[i+1] = options.klength[i];

      tmpOpt1[4] = options.stateInformation[0];
      tmpOpt1[5] = options.stateInformation[1];
      tmpOpt1[6] = options.omega_x.real();
      tmpOpt1[7] = options.omega_y.real();
      tmpOpt1[8] = options.dispersion_x.real();
      tmpOpt1[9] = options.dispersion_y.real();
      tmpOpt1[10] = options.min_x;
      tmpOpt1[11] = options.min_y;
      tmpOpt1[12] = options.t_abs.real();
      tmpOpt1[13] = options.exp_factor.real();
      tmpOpt1[14] = options.g;
      tmpOpt1[15] = options.ITP_step;
      tmpOpt1[16] = options.RTE_step;

      for(int i = 0; i<4;i++)
        tmpOpt1[i+17] = options.grid[i];

      tmpOpt1[21] = options.n_it_RTE;
      tmpOpt1[22] = options.samplesize;
      tmpOpt1[23] = options.vortexnumber;  

      //copy options to file
      hsize_t dimsf[1] = {24};
      dataspace = H5Screate_simple(1, dimsf, NULL);

      h5a_options = H5Acreate(h5_timegroup, "Options", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite (h5a_options, H5T_IEEE_F64LE, tmpOpt1);

      H5Aclose(h5a_options);
      H5Sclose(dataspace);
      // cerr << "Reached ERROR location #13" << endl;

  H5Gclose(h5_timegroup);

  return true;
}


bool binaryFile::getSnapshot(const string &name, double snapShotTime, vector<MatrixXcd> &k, Options &options)
{
  if(m == out)
    {
      cout << "file "<< filename.c_str() << "is not in read mode" << endl;
      return false;
    }

  stringstream time_name;
  // snapShotTime will always be formatted with two decimal digits
  time_name.setf(ios_base::fixed);
  time_name.precision(2);
  time_name << snapShotTime;

  if(H5Lexists(h5_file, (time_name.str()).c_str(), H5P_DEFAULT))
    {
      h5_timegroup = H5Gopen(h5_file, (time_name.str()).c_str(), H5P_DEFAULT);

    }
  else
    {
      cout << "ERROR: HDF5 Group for time " << snapShotTime << " does not exist" << endl;
      return false;
    }

  stringstream set_name;
  set_name << name;

  if(!H5Lexists(h5_timegroup, (set_name.str()).c_str(), H5P_DEFAULT))
    {
      cout << "ERROR: HDF5 Dataset " << name << " does not exist for time " << snapShotTime << " in EigenMatrix branch" << endl;

      H5Gclose(h5_timegroup);
      return false;
    }
  else
    {
      hid_t test_id = H5Oopen(h5_timegroup, (set_name.str()).c_str(), H5P_DEFAULT);

              hid_t h5a_options;
              h5a_options = H5Aopen(h5_timegroup, "Options", H5P_DEFAULT);

              double tmpOpt1[24];

              H5Aread(h5a_options, H5T_IEEE_F64LE , tmpOpt1);

              //load Options struct from file array. Don't forget to change appropriately when changing Options struct
              options.N = tmpOpt1[0];
              for(int i= 0; i < 3; i++){
               options.klength[i] = tmpOpt1[i+1];
              }
              options.stateInformation[0] = tmpOpt1[4];
              options.stateInformation[1] = tmpOpt1[5];
              options.omega_x = complex<double>(tmpOpt1[6],0.0);
              options.omega_y = complex<double>(tmpOpt1[7],0.0);
              options.dispersion_x = complex<double>(tmpOpt1[8],0.0);
              options.dispersion_x = complex<double>(tmpOpt1[9],0.0);
              options.min_x = tmpOpt1[10];
              options.min_y = tmpOpt1[11];
              options.t_abs = complex<double>(tmpOpt1[12],0.0);
              options.exp_factor = complex<double>(tmpOpt1[13],0.0);
              options.g = tmpOpt1[14];
              options.ITP_step = tmpOpt1[15];
              options.RTE_step = tmpOpt1[16];

              for(int i = 0; i<4;i++)
                options.grid[i] = (uint32_t)tmpOpt1[i+17];

              options.n_it_RTE = (int)tmpOpt1[21];
              options.samplesize = (int)tmpOpt1[22];
              options.vortexnumber = (int)tmpOpt1[23];              

              H5Aclose(h5a_options);

      // if(H5Iget_type(test_id) == H5I_DATASET)
      //   {
      //     hid_t dataspace = H5Dget_space(test_id);

      //     hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
      //     hsize_t *dimf = (hsize_t *)malloc(rank*sizeof(hsize_t));
      //     H5Sget_simple_extent_dims(dataspace, dimf, NULL);

      //     dimf[0] /= 2;
      //     dimf[1] /= 2;
      //     k[0].resize(dimf[0],dimf[1]);
          

      //     H5Dread(test_id,  H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)k[0].data());

      //     free(dimf);
      //     H5Dclose(test_id);
      //   }
      // else
      //   {
          hsize_t vecsize;
          H5Gget_num_objs(test_id, &vecsize);
          k.resize(vecsize);

          for(int i = 0; i < vecsize; i++)
            {
              k[i] = MatrixXcd(options.grid[1],options.grid[2]);
              stringstream comp;
              comp << i;

              hid_t dataset = H5Dopen(test_id, (comp.str()).c_str(), H5P_DEFAULT);
              hid_t dataspace = H5Dget_space(dataset);

              hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
              hsize_t *dimf = (hsize_t *)malloc(rank*sizeof(hsize_t));
              H5Sget_simple_extent_dims(dataspace, dimf, NULL);

          	  dimf[0] /= 2;
          	  // dimf[1] /= 2;
              k[i].resize(dimf[0],dimf[1]);

              H5Dread(dataset,  H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)k[i].data());

              free(dimf);
              H5Dclose(dataset);
            }
          H5Gclose(test_id);
        // }
    }



  H5Gclose(h5_timegroup);

  return true;
}


bool binaryFile::appendEval(const string &vec_name, double *vec, int vec_rank, int* vec_dim, double snapShotTime)
{
  if(m == in)
    {
      cout << "file "<< filename.c_str() << "is not in write mode" << endl;
      return false;
    }

  if(!checkTime(snapShotTime))
    {
      cout << "HDF5 Group error for time group: " << snapShotTime << endl;
      return false;
    }

  hid_t h5_avgroup;

  if(H5Lexists(h5_timegroup, "Averages", H5P_DEFAULT))
    h5_avgroup = H5Gopen(h5_timegroup, "Averages", H5P_DEFAULT);
  else
    h5_avgroup = H5Gcreate(h5_timegroup, "Averages", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


  if(H5Lexists(h5_avgroup, vec_name.c_str(), H5P_DEFAULT))
    {
      cout << "Averages " << vec_name << " already exists for time " << snapShotTime << ". I refuse to write." << endl;
      return false;
    }



  hid_t dataset, dataspace, dset_create_props;

  dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list

  hsize_t *dimsf = (hsize_t*)malloc(vec_rank*sizeof(hsize_t));
  for(int l = 0; l < vec_rank; l++)
    {
      dimsf[l] = vec_dim[l];
    }

  dataspace = H5Screate_simple(1, dimsf,  NULL);

  dataset = H5Dcreate(h5_avgroup, vec_name.c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
  H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec); //write grid data

  //append grid snapShotTime as attribute

  free(dimsf);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Pclose(dset_create_props);

  H5Gclose(h5_avgroup);
  H5Gclose(h5_timegroup);

  return true;
}