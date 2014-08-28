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
  filename = "runData/" + file;
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
              int size = (int)(H5Aget_storage_size(h5a_options)/sizeof(int));
              time_list.resize(size);
              H5Aread(h5a_options, H5T_STD_I32LE , &time_list.front());
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
      string dirname = "runData";
      struct stat st;
      if(stat(dirname.c_str(),&st) != 0){
        mkdir(dirname.c_str(),0755);
      }

    //   string fileNameListName = "runData/fileNameList.dat";
    //   // struct stat buffer;   
    //   // if(stat (filename.c_str(), &buffer) != 0){
    //     ofstream fileNameListFile;
    //     fileNameListFile.open(fileNameListName.c_str(), ios::out | ios::app);
    //     fileNameListFile << std::left << file << endl;
    //     fileNameListFile.close();
    // // } 



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
      hid_t h5a = H5Acreate(h5_file, "SnapshotTimes", H5T_STD_I32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite (h5a, H5T_STD_I32LE, &time_list.front());
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
      hid_t h5a = H5Acreate(h5_file, "SnapshotTimes", H5T_STD_I32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite (h5a, H5T_STD_I32LE, &time_list.front());
      H5Sclose(dataspace);
      H5Aclose(h5a);
      // cerr << "Reached ERROR location #6" << endl;
    }


  H5Fclose(h5_file);
  filename = "";
}

bool binaryFile::checkTime(int snapShotTime)
{
  // stringstream time_name;
  // // snapShotTime will always be formatted with two decimal digits
  // time_name.setf(ios_base::fixed);
  // time_name.precision(2);
  // time_name << snapShotTime;

  string time_name = to_string(snapShotTime);

  if(H5Lexists(h5_file, time_name.c_str(), H5P_DEFAULT)){
    h5_timegroup = H5Gopen(h5_file, time_name.c_str(), H5P_DEFAULT);
    // cerr << "Reached ERROR location #7" << endl;
  }
  else
    {
      h5_timegroup = H5Gcreate(h5_file, time_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      time_list.push_back(snapShotTime);
      // cerr << "Reached ERROR location #8" << endl;
    }

  if(h5_timegroup >= 0)
    return true;
  else
    return false;
}


bool binaryFile::appendSnapshot(const string &name, int snapShotTime, MatrixData* const &pData, Options &options)
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
      for (int i = 0; i < pData->wavefunction.size(); i++)
        {
          hid_t    dataset, dataspace, dset_create_props;

          hsize_t rank = 2;
          // rank += 1; //for internal dimension

          hsize_t *dimsf = (hsize_t *) malloc(rank*sizeof(hsize_t));

          dimsf[0] = 2 * pData->wavefunction[i].rows();
          dimsf[1] = pData->wavefunction[i].cols();


          dset_create_props = H5Pcreate(H5P_DATASET_CREATE); //create a default creation property list
          dataspace = H5Screate_simple(rank, dimsf, NULL); //define space in file

          stringstream comp;
          comp << i;

          dataset = H5Dcreate(h5_EMGroup_vecsub, (comp.str()).c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT); //create data block in file

          H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)pData->wavefunction[i].data()); //write grid data

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

      tmpOpt1[21] = options.potFactor;
      tmpOpt1[22] = options.samplesize;
      tmpOpt1[23] = options.vortexnumber;  

      //copy options to file
      hsize_t dimsf[1] = {24};
      dataspace = H5Screate_simple(1, dimsf, NULL);

      h5a_options = H5Acreate(h5_timegroup, "Options", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite (h5a_options, H5T_IEEE_F64LE, tmpOpt1);

      H5Aclose(h5a_options);
      H5Sclose(dataspace);

      hid_t h5a_meta;
      dimsf[0] = 9;
      dataspace = H5Screate_simple(1, dimsf, NULL);

      h5a_meta = H5Acreate(h5_timegroup, "Meta", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(h5a_meta, H5T_IEEE_F64LE, pData->meta.data());
      H5Aclose(h5a_meta);
      H5Sclose(dataspace);
      // cerr << "Reached ERROR location #13" << endl;

  H5Gclose(h5_timegroup);

  return true;
}


bool binaryFile::getSnapshot(const string &name, int snapShotTime, MatrixData* &pData, Options &options)
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

              options.potFactor = tmpOpt1[21];
              options.samplesize = (int)tmpOpt1[22];
              options.vortexnumber = (int)tmpOpt1[23];              

              H5Aclose(h5a_options);

              hid_t h5a_meta;
              h5a_meta = H5Aopen(h5_timegroup, "Meta", H5P_DEFAULT);
              H5Aread(h5a_meta, H5T_IEEE_F64LE, pData->meta.data());
              H5Aclose(h5a_meta);
              pData->meta.arrayToData();

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
          pData->wavefunction.resize(vecsize);

          for(int i = 0; i < vecsize; i++)
            {
              pData->wavefunction[i] = MatrixXcd(pData->meta.grid[0], pData->meta.grid[1]);
              stringstream comp;
              comp << i;

              hid_t dataset = H5Dopen(test_id, (comp.str()).c_str(), H5P_DEFAULT);
              hid_t dataspace = H5Dget_space(dataset);

              hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
              hsize_t *dimf = (hsize_t *)malloc(rank*sizeof(hsize_t));
              H5Sget_simple_extent_dims(dataspace, dimf, NULL);

          	  // dimf[0] /= 2;
          	  // dimf[1] /= 2;
              // k[i].resize(dimf[0],dimf[1]);

              H5Dread(dataset,  H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)pData->wavefunction[i].data());

              free(dimf);
              H5Dclose(dataset);
            }
          H5Gclose(test_id);
        // }
    }



  H5Gclose(h5_timegroup);

  return true;
}

bool binaryFile::appendEval(int snapShotTime, Options options, MatrixData::MetaData meta, Eval results)
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

    if(!H5Lexists(h5_timegroup, "Options", H5P_DEFAULT)){
          hid_t    h5a_options, dataspace_opt;
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

      tmpOpt1[21] = options.potFactor;
      tmpOpt1[22] = options.samplesize;
      tmpOpt1[23] = options.vortexnumber;  

      //copy options to file
      hsize_t dimsf[1] = {24};
      dataspace_opt = H5Screate_simple(1, dimsf, NULL);

      h5a_options = H5Acreate(h5_timegroup, "Options", H5T_IEEE_F64LE, dataspace_opt, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite (h5a_options, H5T_IEEE_F64LE, tmpOpt1);

      H5Aclose(h5a_options);
      H5Sclose(dataspace_opt);

      if(!H5Lexists(h5_timegroup, "Meta", H5P_DEFAULT)){

      hid_t h5a_meta;
      dimsf[0] = 9;
      dataspace_opt = H5Screate_simple(1, dimsf, NULL);

      h5a_meta = H5Acreate(h5_timegroup, "Meta", H5T_IEEE_F64LE, dataspace_opt, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(h5a_meta, H5T_IEEE_F64LE, meta.data());
      H5Aclose(h5a_meta);
    }

      H5Sclose(dataspace_opt);
    }

  

  hid_t h5_observables;

  if(H5Lexists(h5_timegroup, "Observables", H5P_DEFAULT))
    h5_observables = H5Gopen(h5_timegroup, "Observables", H5P_DEFAULT);
  else
    h5_observables = H5Gcreate(h5_timegroup, "Observables", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /// PREPARE THE DOUBLE ARRAYS
      string vec1Name = "Averages";
      int vec1Rank = 11;
      double vec1[11];
      vec1[0] = results.totalResult.Ekin;
      vec1[1] = results.totalResult.particle_count;
      vec1[2] = results.totalResult.healing_length;
      vec1[3] = results.totalResult.volume;
      vec1[4] = results.totalResult.density;
      vec1[5] = results.totalResult.aspectRatio;
      vec1[6] = results.totalResult.aspectRatioAngle;
      vec1[7] = results.totalResult.r_max;
      vec1[8] = results.totalResult.r_min;
      vec1[9] = results.totalResult.r_max_phi;
      vec1[10] = results.totalResult.r_min_phi;

  if(H5Lexists(h5_observables, vec1Name.c_str(), H5P_DEFAULT))
    {
      cout << "Observables " << vec1Name << " already exists for time " << snapShotTime << ". I refuse to write." << endl;
      return false;
    }

  hid_t dataset, dataspace, dset_create_props;

  dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list

  hsize_t *dimsf = (hsize_t*)malloc(sizeof(hsize_t));
    dimsf[0] = vec1Rank;


  dataspace = H5Screate_simple(1, dimsf,  NULL);

  dataset = H5Dcreate(h5_observables, vec1Name.c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
  H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec1); //write grid data

  //append grid snapShotTime as attribute


  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Pclose(dset_create_props);



      
  string vec4Name = "KVector";
  string vec5Name = "OccupationNumber";
  int vec4Rank = results.totalResult.k.size();
  int vec5Rank = results.totalResult.number.size();
  vector<double> vec4;
  vector<double> vec5;


  int sharedRank;
  if(vec4Rank == vec5Rank){
    sharedRank = vec4Rank;
  } else {
    cout << "Ranks of Kvector and OccupationNumber are not the same!" << endl;
    sharedRank = (vec4Rank >= vec5Rank) ? vec4Rank : vec5Rank;
  }


  for(int i = 0; i < sharedRank; i++){
    if((results.totalResult.k(i) != 0) || (results.totalResult.number(i) != 0))
      vec4.push_back(results.totalResult.k(i));
      vec5.push_back(results.totalResult.number(i));
  }
  vec4Rank = vec4.size();
  vec5Rank = vec5.size();

  cout << "Vector Ranks: " << vec4Rank << " " << vec5Rank << endl;

  if(H5Lexists(h5_observables, vec4Name.c_str(), H5P_DEFAULT))
    {
      cout << "Observables " << vec4Name << " already exists for time " << snapShotTime << ". I refuse to write." << endl;
      return false;
    }


  dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list

  
    dimsf[0] = vec4Rank;


  dataspace = H5Screate_simple(1, dimsf,  NULL);

  dataset = H5Dcreate(h5_observables, vec4Name.c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
  H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec4.data()); //write grid data

  //append grid snapShotTime as attribute


  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Pclose(dset_create_props);

  if(H5Lexists(h5_observables, vec5Name.c_str(), H5P_DEFAULT))
    {
      cout << "Observables " << vec5Name << " already exists for time " << snapShotTime << ". I refuse to write." << endl;
      return false;
    }


  dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list

    dimsf[0] = vec5Rank;


  dataspace = H5Screate_simple(1, dimsf,  NULL);

  dataset = H5Dcreate(h5_observables, vec5Name.c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
  H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec5.data()); //write grid data

  //append grid snapShotTime as attribute


  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Pclose(dset_create_props);

  
  





  string vec2Name = "Contour";
  int samples = results.contour.size();
  int vec2Ranks[samples];
  int arraysize = 1 + samples; // first elements are the number of samples, and the size of each sample
  double *vec2;

  for(int k = 0; k < samples; k++){
    vec2Ranks[k] = results.contour[k].size();
    arraysize += 2 * vec2Ranks[k]; // x and y coordinate of each contourpoint summed over all samples
  }

  vec2 = (double*)malloc(sizeof(double)*arraysize);
  vec2[0] = samples;

  int l = 1 + samples;
  for(int k = 0; k < samples; k++){
    vec2[k+1] = vec2Ranks[k];    
    for(std::unordered_set<Coordinate<int32_t>,Hash>::const_iterator it = results.contour[k].begin(); it != results.contour[k].end(); ++it){
      vec2[l] = it->x();
      vec2[l+1] = it->y();
      l+=2;
    }
  }

    if(H5Lexists(h5_observables, vec2Name.c_str(), H5P_DEFAULT))
    {
      cout << "Observables " << vec2Name << " already exists for time " << snapShotTime << ". I refuse to write." << endl;
      return false;
    }

    dset_create_props = H5Pcreate(H5P_DATASET_CREATE); // default creation property list;
    dimsf[0] = arraysize;

    dataspace = H5Screate_simple(1,dimsf,NULL);
    dataset = H5Dcreate(h5_observables,vec2Name.c_str(),H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
    H5Dwrite(dataset,H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec2); 

    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Pclose(dset_create_props);
    free(vec2);

  string vec3Name = "Vortices";
  samples = results.pres.size();
  int vec3Ranks[samples];
  arraysize = 1 + samples; // first elements are the number of samples, and the size of each sample
  double *vec3;

  for(int k = 0; k < samples; k++){
    vec3Ranks[k] = results.pres[k].vlist.size();
    arraysize += 2 * vec3Ranks[k]; // x and y coordinate of each vortex summed over all samples
  }

  vec3 = (double*)malloc(sizeof(double)*arraysize);
  vec3[0] = samples;

  l = 1 + samples;
  for(int k = 0; k < samples; k++){
    vec3[k+1] = vec3Ranks[k];    
    for(std::list<VortexData>::const_iterator it = results.pres[k].vlist.begin(); it != results.pres[k].vlist.end(); ++it){
      vec3[l] = it->x.x();
      vec3[l+1] = it->x.y();
      l+=2;
    }
  }

    if(H5Lexists(h5_observables, vec3Name.c_str(), H5P_DEFAULT))
    {
      cout << "Observables " << vec3Name << " already exists for time " << snapShotTime << ". I refuse to write." << endl;
      return false;
    }

    dset_create_props = H5Pcreate(H5P_DATASET_CREATE); // default creation property list;
    dimsf[0] = arraysize;

    dataspace = H5Screate_simple(1,dimsf,NULL);
    dataset = H5Dcreate(h5_observables,vec3Name.c_str(),H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
    H5Dwrite(dataset,H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec3); 

    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Pclose(dset_create_props);
    free(vec3);


  free(dimsf);
  H5Gclose(h5_observables);
  H5Gclose(h5_timegroup);

  return true;
}