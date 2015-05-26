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
  m = nm;
 
  if(m == in || m == append)
	{
	  struct stat buf;
	  lstat(filename.c_str(), &buf);
	  if(S_ISREG(buf.st_mode)) // Nur normale Dateien ueberpruefen
		{ 
		  if(H5Fis_hdf5(filename.c_str())){
		  	if(m == in)
				h5_file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			  	else if(m == append)
				h5_file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	
			  	if(H5Aexists(h5_file, "SnapshotTimes") > 0){
	
					hid_t h5a_times;
	
					h5a_times = H5Aopen(h5_file, "SnapshotTimes", H5P_DEFAULT);
					int size = (int)(H5Aget_storage_size(h5a_times)/sizeof(int));
					time_list.resize(size);
					H5Aread(h5a_times, H5T_STD_I32LE , &time_list.front());
					H5Aclose(h5a_times);
					}
				  	else
					{
					  cout << "WARNING for I/O operation on file: "<< filename << " occurred: No valid SnapshotTimes attribute found in file." << endl;
				}
			}
			else {
				cout << "ERROR for I/O operation. File: "<< filename <<" is not in hdf5 storage format" << endl;
			}
		}
	  else
		{
		  cout << "ERROR for I/O operation. File: "<< filename <<" is not regular" << endl;
		}
	}
  else
	{
	  string dirname = "runData";
	  struct stat st;
	  if(stat(dirname.c_str(),&st) != 0){
		mkdir(dirname.c_str(),0755);
	  }
	  h5_file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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

bool binaryFile::checkTime(int snapShotTime){

  string time_name = to_string(snapShotTime);

  if(H5Lexists(h5_file, time_name.c_str(), H5P_DEFAULT)){
	h5_timegroup = H5Gopen(h5_file, time_name.c_str(), H5P_DEFAULT);
  }
  else{
	h5_timegroup = H5Gcreate(h5_file, time_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	time_list.push_back(snapShotTime);
  }

  if(h5_timegroup >= 0)
	return true;
  else
	return false;
}




void binaryFile::writeMatrixData(const string &name, MatrixData const * const &pData ){

	hid_t h5_EMGroup_vecsub = H5Gcreate(h5_timegroup, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); //create dataset with full space selection

	for (int i = 0; i < pData->wavefunction.size(); i++){
		hid_t dataset, dataspace, dset_create_props;
	
		hsize_t rank = 2;
	
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
	H5Gclose(h5_EMGroup_vecsub);


}

void binaryFile::writeMeta(MatrixData::MetaData &meta ){
	if(!H5Aexists(h5_timegroup, "Meta")){
		hsize_t dimsf[1] = {(hsize_t)meta.size()};
		hid_t dataspace = H5Screate_simple(1, dimsf, NULL);
	
		hid_t h5a_meta = H5Acreate(h5_timegroup, "Meta", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(h5a_meta, H5T_IEEE_F64LE, meta.data());
		H5Aclose(h5a_meta);
		H5Sclose(dataspace);
	}
}

void binaryFile::writeOptions(Options const & options){
  if(!H5Aexists(h5_timegroup, "Options")){
	hid_t    h5a_options, dataspace;

	double tmpOpt1[30];

	tmpOpt1[0] = options.N;
	for(int i= 0; i < 3; i++)
	  tmpOpt1[i+1] = options.klength[i];

	 tmpOpt1[4] = options.stateInformation[0];
	 tmpOpt1[5] = options.stateInformation[1];
	 tmpOpt1[6] = options.stateInformation[2];
	 tmpOpt1[7] = options.omega_x.real();
	 tmpOpt1[8] = options.omega_y.real();
	 tmpOpt1[9] = options.omega_w.real();
	 tmpOpt1[10] = options.dispersion_x.real();
	 tmpOpt1[11] = options.dispersion_y.real();
	 tmpOpt1[12] = options.min_x;
	 tmpOpt1[13] = options.min_y;
	 tmpOpt1[14] = options.min_z;
	 // tmpOpt1[15] = options.t_abs.real();
	 tmpOpt1[16] = options.exp_factor.real();
	 tmpOpt1[17] = options.g;
	 tmpOpt1[18] = options.ITP_step;
	 tmpOpt1[19] = options.RTE_step;


	for(int i = 0; i<4;i++)
	  tmpOpt1[i+20] = options.grid[i];

	tmpOpt1[24] = options.potFactor;
	tmpOpt1[25] = options.samplesize;
	tmpOpt1[26] = options.vortexnumber;
	tmpOpt1[27] = options.vortexspacing;
	tmpOpt1[28] = options.Ag;
	tmpOpt1[29] = options.OmegaG;    

	//copy options to file
	hsize_t dimsf[1] = {30};
	dataspace = H5Screate_simple(1, dimsf, NULL);

	h5a_options = H5Acreate(h5_timegroup, "Options", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite (h5a_options, H5T_IEEE_F64LE, tmpOpt1);

	H5Aclose(h5a_options);
	H5Sclose(dataspace);
  }

  
}


bool binaryFile::appendSnapshot(const string &name, MatrixData * const &pData, Options const &options){
	int snapShotTime = pData->meta.steps;
  if(m == in){
	cout << "file "<< filename.c_str() << "is not in write mode" << endl;
	return false;
  }

  if(!checkTime(snapShotTime)){
	cout << "HDF5 Group error for snapShotTime group: " << snapShotTime << endl;
	H5Gclose(h5_timegroup);
	return false;
  }

  writeMatrixData(name, pData);
  writeMeta(pData->meta);
  writeOptions(options);

  H5Gclose(h5_timegroup);

  return true;
}

void binaryFile::readOptions(Options &options){
	
	hid_t h5a_options;
	h5a_options = H5Aopen(h5_timegroup, "Options", H5P_DEFAULT);

	// int sizeOfOptions = (int)(H5Aget_storage_size(h5a_options)/sizeof(double));

	double tmpOpt1[30];

	H5Aread(h5a_options, H5T_IEEE_F64LE , tmpOpt1);

	//load Options struct from file array. Don't forget to change appropriately when changing Options struct
	options.N = tmpOpt1[0];
	for(int i= 0; i < 3; i++){
	 options.klength[i] = tmpOpt1[i+1];
	}
	options.stateInformation[0] = tmpOpt1[4];
	options.stateInformation[1] = tmpOpt1[5];
	options.stateInformation[2] = tmpOpt1[6];
	options.omega_x = complex<double>(tmpOpt1[7],0.0);
	options.omega_y = complex<double>(tmpOpt1[8],0.0);
	options.omega_w = complex<double>(tmpOpt1[9],0.0);
	options.dispersion_x = complex<double>(tmpOpt1[10],0.0);
	options.dispersion_y = complex<double>(tmpOpt1[11],0.0);
	options.min_x = tmpOpt1[12];
	options.min_y = tmpOpt1[13];
	options.min_z = tmpOpt1[14];
	// options.t_abs = complex<double>(tmpOpt1[15],0.0);
	options.exp_factor = complex<double>(tmpOpt1[16],0.0);
	options.g = tmpOpt1[17];
	options.ITP_step = tmpOpt1[18];
	options.RTE_step = tmpOpt1[19];

	for(int i = 0; i<4;i++)
	  options.grid[i] = (uint32_t)tmpOpt1[i+20];

	options.potFactor = tmpOpt1[24];
	options.samplesize = (int)tmpOpt1[25];
	options.vortexnumber = (int)tmpOpt1[26];
	options.vortexspacing = (int)tmpOpt1[27];
	options.Ag = tmpOpt1[28];
	options.OmegaG = tmpOpt1[29];              

	H5Aclose(h5a_options);
}

void binaryFile::readMeta(MatrixData::MetaData &meta){

	hid_t h5a_meta;
	h5a_meta = H5Aopen(h5_timegroup, "Meta", H5P_DEFAULT);
	H5Aread(h5a_meta, H5T_IEEE_F64LE, meta.data());
	H5Aclose(h5a_meta);
	meta.arrayToData(); 
}

void binaryFile::readMatrixData(string const &name, MatrixData* &pData){

  stringstream set_name;
  set_name << name;

  if(!H5Lexists(h5_timegroup, (set_name.str()).c_str(), H5P_DEFAULT)){
	cout << "ERROR: HDF5 Dataset " << name << " does not exist for name " << name << " in EigenMatrix branch" << endl;
  }else{
	hid_t h5_EMGroup_vecsub = H5Oopen(h5_timegroup, (set_name.str()).c_str(), H5P_DEFAULT);
		  
	hsize_t vecsize;
	H5Gget_num_objs(h5_EMGroup_vecsub, &vecsize);
	pData->wavefunction.resize(vecsize);

	for(int i = 0; i < vecsize; i++){
	  pData->wavefunction[i] = MatrixXcd(pData->meta.grid[0], pData->meta.grid[1]);
	  stringstream comp;
	  comp << i;

	  hid_t dataset = H5Dopen(h5_EMGroup_vecsub, (comp.str()).c_str(), H5P_DEFAULT);
	  hid_t dataspace = H5Dget_space(dataset);

	  hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
	  hsize_t *dimf = (hsize_t *)malloc(rank*sizeof(hsize_t));
	  H5Sget_simple_extent_dims(dataspace, dimf, NULL);

	  H5Dread(dataset,  H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)pData->wavefunction[i].data());

	  free(dimf);
	  H5Dclose(dataset);
	}
	H5Gclose(h5_EMGroup_vecsub);
  }
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

  if(H5Lexists(h5_file, (time_name.str()).c_str(), H5P_DEFAULT)){
	h5_timegroup = H5Gopen(h5_file, (time_name.str()).c_str(), H5P_DEFAULT);
  }else{
	cout << "ERROR: HDF5 Group for time " << snapShotTime << " does not exist" << endl;
	return false;
  }

  readOptions(options);
  readMeta(pData->meta);
  readMatrixData(name, pData);

  H5Gclose(h5_timegroup);

  return true;
}

bool binaryFile::getLatestSnapshot(const string &name, MatrixData* &pData, Options &options){
	if(m == out){
		cout << "file " << filename.c_str() << " is not in read mode" << endl;
		return false;
	}

	stringstream time_name;
	time_name.setf(ios_base::fixed);
	time_name.precision(2);
	time_name << time_list.back();

	if(H5Lexists(h5_file, (time_name.str()).c_str(), H5P_DEFAULT)){
		h5_timegroup = H5Gopen(h5_file, (time_name.str()).c_str(), H5P_DEFAULT);
  	}else{
		cout << "ERROR: HDF5 Group for latest time " << time_list.back() << " does not exist" << endl;
		return false;
  	}

  	readOptions(options);
  	readMeta(pData->meta);
  	readMatrixData(name, pData);

  	H5Gclose(h5_timegroup);

}









bool binaryFile::appendEval(Eval &results, Options const & options){
	int snapShotTime = results.data.meta.steps;

  if(m == in){
	cout << "file "<< filename.c_str() << "is not in write mode" << endl;
	return false;
  }

  if(!checkTime(snapShotTime)){
	cout << "HDF5 Group error for time group: " << snapShotTime << endl;
	return false;
  }

  writeOptions(options);
  writeMeta(results.data.meta);  

  hid_t h5_observables;

  if(H5Lexists(h5_timegroup, "Observables", H5P_DEFAULT))
	h5_observables = H5Gopen(h5_timegroup, "Observables", H5P_DEFAULT);
  else
	h5_observables = H5Gcreate(h5_timegroup, "Observables", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	/// PREPARE THE DOUBLE ARRAYS
	string vec1Name = "Averages";
	int vec1Rank = 13;
	double vec1[13];
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
	vec1[11] = results.totalResult.Rx;
	vec1[12] = results.totalResult.Ry;

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
	string vec6Name = "AspectRatios";
	int vec4Rank = results.totalResult.k.size();
	int vec5Rank = results.totalResult.number.size();
	int vec6Rank = results.totalResult.fixedAspectRatio.size();
	vector<double> vec4;
	vector<double> vec5;
	vector<double> vec6;

	int sharedRank;
	if(vec4Rank == vec5Rank){
		sharedRank = vec4Rank;
	} else {    
		sharedRank = (vec4Rank >= vec5Rank) ? vec4Rank : vec5Rank;
		cout << "Ranks of Kvector and OccupationNumber Input are not the same!" << vec4Rank <<" "<< vec5Rank << endl;
	}

	for(int i = 0; i < sharedRank; i++){
		if(results.totalResult.number(i) != 0.0){
			vec4.push_back(results.totalResult.k(i));
			vec5.push_back(results.totalResult.number(i));
		}
	}
	for(int i = 0; i < vec6Rank; i++){
		vec6.push_back(results.totalResult.fixedAspectRatio(i));
	}
	vec4Rank = vec4.size();
	vec5Rank = vec5.size();
	vec6Rank = vec6.size();


  	if(H5Lexists(h5_observables, vec4Name.c_str(), H5P_DEFAULT)){
		cout << "Observables " << vec4Name << " already exists for time " << snapShotTime << ". I refuse to write." << endl;
  	}else {	
		dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list
		
		dimsf[0] = vec4Rank;
	
		dataspace = H5Screate_simple(1, dimsf,  NULL);
	
		dataset = H5Dcreate(h5_observables, vec4Name.c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
		H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec4.data()); //write grid data
	
		//append grid snapShotTime as attribute
	
		H5Dclose(dataset);
		H5Sclose(dataspace);
		H5Pclose(dset_create_props);
	}

  	if(H5Lexists(h5_observables, vec5Name.c_str(), H5P_DEFAULT)){
		cout << "Observables " << vec5Name << " already exists for time " << snapShotTime << ". I refuse to write." << endl;
	} else {	
		dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list
	
		dimsf[0] = vec5Rank;
	
		dataspace = H5Screate_simple(1, dimsf,  NULL);
	
		dataset = H5Dcreate(h5_observables, vec5Name.c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
		H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec5.data()); //write grid data
	
		//append grid snapShotTime as attribute
	
		H5Dclose(dataset);
		H5Sclose(dataspace);
		H5Pclose(dset_create_props);
	}

 	if(H5Lexists(h5_observables, vec6Name.c_str(), H5P_DEFAULT)){
		cout << "Observables " << vec6Name << " already exists for time " << snapShotTime << ". I refuse to write." << endl;
  	} else {	
		dset_create_props = H5Pcreate (H5P_DATASET_CREATE); //create a default creation property list
		
		dimsf[0] = vec6Rank;
		
		dataspace = H5Screate_simple(1, dimsf,  NULL);
		
		dataset = H5Dcreate(h5_observables, vec6Name.c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
		H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec6.data()); //write grid data
		
		//append grid snapShotTime as attribute
		
		H5Dclose(dataset);
		H5Sclose(dataspace);
		H5Pclose(dset_create_props);
	}



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

  if(H5Lexists(h5_observables, vec2Name.c_str(), H5P_DEFAULT)){
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

  if(H5Lexists(h5_observables, vec3Name.c_str(), H5P_DEFAULT)){
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

  string vec7Name = "VortexWinding";
  samples = results.pres.size();
  int vec7Ranks[samples];
  arraysize = 1 + samples;
  double *vec7;

  for(int k = 0; k < samples; k++){
  	vec7Ranks[k] = results.pres[k].vlist.size();
  	arraysize += vec7Ranks[k];
  }

  vec7 = (double*)malloc(sizeof(double)*arraysize);
  vec7[0] = samples;

  l = 1 + samples;
  for(int k = 0; k < samples; k++){
	vec7[k+1] = vec7Ranks[k];    
	for(std::list<VortexData>::const_iterator it = results.pres[k].vlist.begin(); it != results.pres[k].vlist.end(); ++it){
	  vec7[l] = it->n;
	  l++;
	}
  }


  if(H5Lexists(h5_observables, vec7Name.c_str(), H5P_DEFAULT)){
	cout << "Observables " << vec7Name << " already exists for time " << snapShotTime << ". I refuse to write." << endl;
	return false;
  }

  dset_create_props = H5Pcreate(H5P_DATASET_CREATE); // default creation property list;
  dimsf[0] = arraysize;

  dataspace = H5Screate_simple(1,dimsf,NULL);
  dataset = H5Dcreate(h5_observables,vec7Name.c_str(),H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT);
  H5Dwrite(dataset,H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec7); 

  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Pclose(dset_create_props);
  free(vec7);

  free(dimsf);
  H5Gclose(h5_observables);
  H5Gclose(h5_timegroup);

  return true;
}


bool binaryFile::getEval(int snapShotTime, Eval &results, Options &options){

	if(m == out){
		cout << "file "<< filename.c_str() << "is not in read mode" << endl;
		return false;
	}

	string time_name = to_string(snapShotTime);

	if(H5Lexists(h5_file, (time_name).c_str(), H5P_DEFAULT)){
		h5_timegroup = H5Gopen(h5_file, (time_name).c_str(), H5P_DEFAULT);
	}else{
		cout << "ERROR: HDF5 Group for time " << snapShotTime << " does not exist" << endl;
		return false;
	}

	readOptions(options);
	readMeta(results.data.meta);

	hid_t h5_observables;
	
	if(H5Lexists(h5_timegroup, "Observables", H5P_DEFAULT)){
		h5_observables = H5Gopen(h5_timegroup, "Observables", H5P_DEFAULT);
	}else{
		cout << "Cannot open Group Observables for time " << snapShotTime << ", it doesn't exist." << endl;
		return false;
	}

	hid_t dataset, dataspace ;
	hsize_t *dimf;
	hsize_t rank;
	
	string vec1Name = "Averages";
	if(!H5Lexists(h5_observables, vec1Name.c_str(), H5P_DEFAULT)){
		cout << "Observables " << vec1Name << " doesn't exists for time " << snapShotTime << ". Cannot read." << endl;
	} else {		
		int vec1Rank = 13;
		vector<double> vec1(vec1Rank);
		
		dataset = H5Dopen(h5_observables, vec1Name.c_str(), H5P_DEFAULT);
		dataspace = H5Dget_space(dataset);
		
		rank = H5Sget_simple_extent_ndims(dataspace);
		
		dimf = (hsize_t*)malloc(rank*sizeof(hsize_t));
		H5Sget_simple_extent_dims(dataspace, dimf, NULL);
		
		H5Dread(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec1.data());		
		H5Dclose(dataset);
		H5Sclose(dataspace);
	
	  	results.totalResult.Ekin = vec1[0];
		results.totalResult.particle_count = vec1[1];
		results.totalResult.healing_length = vec1[2];
		results.totalResult.volume = vec1[3];
		results.totalResult.density = vec1[4];
		results.totalResult.aspectRatio = vec1[5];
		results.totalResult.aspectRatioAngle = vec1[6];
		results.totalResult.r_max = vec1[7];
		results.totalResult.r_min = vec1[8];
		results.totalResult.r_max_phi = vec1[9];
		results.totalResult.r_min_phi = vec1[10];
		results.totalResult.Rx = vec1[11];
		results.totalResult.Ry = vec1[12];
	}

	 bool readState[2];

	string vec4Name = "KVector";
	int vec4Rank;	
	vector<double> vec4;
	if(!H5Lexists(h5_observables, vec4Name.c_str(), H5P_DEFAULT)){
		cout << "Observables " << vec4Name << " doesn't exist for time " << snapShotTime << ". Cannot read." << endl;
		readState[0] = false;
	} else {		
		dataset = H5Dopen(h5_observables, vec4Name.c_str(), H5P_DEFAULT);
		dataspace = H5Dget_space(dataset);
		
		rank = H5Sget_simple_extent_ndims(dataspace);
		
		if(rank != 1){
			cout << "reading KVector: wrong dimensions.";
		}
		
		dimf = (hsize_t*)malloc(rank*sizeof(hsize_t));
		H5Sget_simple_extent_dims(dataspace, dimf, NULL);
		vec4.resize(dimf[0]);
		vec4Rank = dimf[0];
		
		H5Dread(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec4.data());
		
		H5Dclose(dataset);
		H5Sclose(dataspace);
		readState[0] = true;
	}

	string vec5Name = "OccupationNumber";
	int vec5Rank;
	vector<double> vec5;
	if(!H5Lexists(h5_observables, vec5Name.c_str(), H5P_DEFAULT)){
		cout << "Observables " << vec5Name << " doesn't exists for time " << snapShotTime << ". Cannot read." << endl;
		readState[1] = false;
	} else {
		dataset = H5Dopen(h5_observables, vec5Name.c_str(), H5P_DEFAULT);
		dataspace = H5Dget_space(dataset);
		
		rank = H5Sget_simple_extent_ndims(dataspace);
		
		if(rank != 1){
			cout << "reading Occupation Numbers: wrong dimensions.";
		}
		
		dimf = (hsize_t*)malloc(rank*sizeof(hsize_t));
		H5Sget_simple_extent_dims(dataspace, dimf, NULL);
		vec5.resize(dimf[0]);
		vec5Rank = dimf[0];
		
		H5Dread(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec5.data());
		
		H5Dclose(dataset);
		H5Sclose(dataspace);
		readState[1] = true;
	}
	if(readState[0] && readState[1]){
		int sharedRank;
		if(vec4Rank == vec5Rank){
			sharedRank = vec4Rank;
	
			results.totalResult.k.resize(vec4Rank);
			results.totalResult.number.resize(vec5Rank);
	
			for(int i = 0; i < sharedRank; i++){
				results.totalResult.k(i) = vec4[i];
				results.totalResult.number(i) = vec5[i];	
			}
	
		} else {    
			sharedRank = (vec4Rank >= vec5Rank) ? vec4Rank : vec5Rank;
			cout << "Size of Kvector and OccupationNumber Input are not the same!" << vec4Rank <<" "<< vec5Rank << endl;
		}
	}

	string vec6Name = "AspectRatios";
	
	vector<double> vec6;
	if(!H5Lexists(h5_observables, vec6Name.c_str(), H5P_DEFAULT)){
		cout << "Observables " << vec6Name << " doesn't exists for time " << snapShotTime << ". Cannot read." << endl;
	} else {
		int vec6Rank;
		dataset = H5Dopen(h5_observables, vec6Name.c_str(), H5P_DEFAULT);
		dataspace = H5Dget_space(dataset);
		
		rank = H5Sget_simple_extent_ndims(dataspace);
		
		if(rank != 1){
			cout << "reading Aspect Ratios: wrong dimensions.";
		}
		
		dimf = (hsize_t*)malloc(rank*sizeof(hsize_t));
		H5Sget_simple_extent_dims(dataspace, dimf, NULL);
		vec6.resize(dimf[0]);
		vec6Rank = dimf[0];
		
		H5Dread(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec6.data());
		
		H5Dclose(dataset);
		H5Sclose(dataspace);

		results.totalResult.fixedAspectRatio.resize(vec6Rank);
		for(int i = 0; i < vec6Rank; i++){
			results.totalResult.fixedAspectRatio(i) = vec6[i];
		}
	}

	string vec2Name = "Contour";
	vector<double> vec2;

	if(!H5Lexists(h5_observables, vec2Name.c_str(), H5P_DEFAULT)){
		cout << "Observables " << vec2Name << " doesn't exists for time " << snapShotTime << ". Cannot read." << endl;
	} else {
		dataset = H5Dopen(h5_observables, vec2Name.c_str(), H5P_DEFAULT);
		dataspace = H5Dget_space(dataset);
		rank = H5Sget_simple_extent_ndims(dataspace);
		if(rank != 1){
			cout << "reading Contour: wrong dimensions.";
		}
		dimf = (hsize_t*)malloc(rank*sizeof(hsize_t));
		H5Sget_simple_extent_dims(dataspace, dimf, NULL);
		vec2.resize(dimf[0]);
		H5Dread(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec2.data());
		H5Dclose(dataset);
		H5Sclose(dataspace);
		int samples = vec2[0];
		int vec2Ranks[samples];
		results.contour.resize(samples);
		for(int i = 0; i < samples; i++){
			vec2Ranks[i] = vec2[i+1];
		}
		
		int l = 1 + samples;
		for(int k = 0; k < samples; k++){
			results.contour[k].clear();
			for(int j = 0; j < vec2Ranks[k]; j++){
				Coordinate<int32_t> c(vec2[l],vec2[l+1],0,options.grid[1],options.grid[2],options.grid[3]);
				l+=2;
				results.contour[k].insert(c);
			}
		}
	}

	string vec3Name = "Vortices";
	vector<double> vec3;

	if(!H5Lexists(h5_observables, vec3Name.c_str(), H5P_DEFAULT)){
	cout << "Observables " << vec3Name << " doesn't exists for time " << snapShotTime << ". Cannot read." << endl;
	} else {
		dataset = H5Dopen(h5_observables, vec3Name.c_str(), H5P_DEFAULT);
		dataspace = H5Dget_space(dataset);
		rank = H5Sget_simple_extent_ndims(dataspace);
		if(rank != 1){
			cout << "reading Vortices: wrong dimensions";
		}
		dimf = (hsize_t*)malloc(rank*sizeof(hsize_t));
		H5Sget_simple_extent_dims(dataspace, dimf, NULL);
		vec3.resize(dimf[0]);
		H5Dread(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec3.data());
		H5Dclose(dataset);
		H5Sclose(dataspace);
		int samples = vec3[0];
		int vec3Ranks[samples];
		results.pres.resize(samples);
		for(int i = 0; i < samples; i++){
			vec3Ranks[i] = vec3[i+1];
		}

		string vec7Name = "VortexWinding";
		vector<double> vec7;
	
		if(!H5Lexists(h5_observables, vec7Name.c_str(), H5P_DEFAULT)){
		cout << "Observables " << vec7Name << " doesn't exists for time " << snapShotTime << ". Cannot read." << endl;
		} else {
			dataset = H5Dopen(h5_observables, vec7Name.c_str(), H5P_DEFAULT);
			dataspace = H5Dget_space(dataset);
			rank = H5Sget_simple_extent_ndims(dataspace);
			if(rank != 1){
				cout << "reading VortexWinding: wrong dimensions";
			}
			dimf = (hsize_t*)malloc(rank*sizeof(hsize_t));
			H5Sget_simple_extent_dims(dataspace, dimf, NULL);
			vec7.resize(dimf[0]);
			H5Dread(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec7.data());
			H5Dclose(dataset);
			H5Sclose(dataspace);
			int vec7Ranks[samples];
			for(int i = 0; i < samples; i++){
				vec7Ranks[i] = vec7[i+1];
				if(vec7Ranks[i] != vec3Ranks[i]){
					cout << "Different number of elements in Vortices and VortexWinding, ignoring VortexWinding";
				}
			}
		}


		
		int l = 1 + samples;
		for(int k = 0; k < samples; k++){
			results.pres[k].vlist.clear();
			for(int j = 0; j < vec3Ranks[k]; j++){
				VortexData c;
				Coordinate<int32_t> x(vec3[l],vec3[l+1],0,options.grid[1],options.grid[2],options.grid[3]);
				c.x = x;
				c.n = (int) vec7[j];
				// c.x.x() = vec3[l];
				// c.x.y() = vec3[l+1];
				l+=2;
				results.pres[k].vlist.push_back(c);
			}
		}
	}

	free(dimf);
	H5Gclose(h5_timegroup);

	return true;
}



// bool binaryFile::appendSnapshot(const string &name, int snapShotTime, vector<ComplexGrid> &data, MatrixData::MetaData &meta, Options &options){

//   if(m == in){
// 	cout << "file "<< filename.c_str() << "is not in write mode" << endl;
// 	return false;
//   }

//   if(!checkTime(snapShotTime)){
// 	cout << "HDF5 Group error for snapShotTime group: " << snapShotTime << endl;
// 	H5Gclose(h5_timegroup);
// 	return false;
//   }

//   hid_t h5_EMGroup_vecsub = H5Gcreate(h5_timegroup, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); //create dataset with full space selection

//   for (int i = 0; i < data.size(); i++){
// 	hid_t    dataset, dataspace, dset_create_props;

// 	hsize_t rank = 3;

// 	hsize_t *dimsf = (hsize_t *) malloc(rank*sizeof(hsize_t));

// 	dimsf[0] = 2 * data[i].width();
// 	dimsf[1] = data[i].height();
// 	dimsf[2] = data[i].depth();

// 	dset_create_props = H5Pcreate(H5P_DATASET_CREATE); //create a default creation property list
// 	dataspace = H5Screate_simple(rank, dimsf, NULL); //define space in file

// 	stringstream comp;
// 	comp << i;

// 	dataset = H5Dcreate(h5_EMGroup_vecsub, (comp.str()).c_str(), H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, dset_create_props, H5P_DEFAULT); //create data block in file

// 	H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)data[i].get_address()); //write grid data

// 	free(dimsf);
// 	H5Dclose(dataset);
// 	H5Pclose(dset_create_props);
// 	H5Sclose(dataspace);
//   }
//   H5Gclose(h5_EMGroup_vecsub);

//   if(!H5Lexists(h5_timegroup, "Options", H5P_DEFAULT)){
// 	hid_t    h5a_options, dataspace;
// 	double tmpOpt1[30];

// 	tmpOpt1[0] = options.N;
// 	for(int i= 0; i < 3; i++)
// 	  tmpOpt1[i+1] = options.klength[i];

// 	 tmpOpt1[4] = options.stateInformation[0];
// 	 tmpOpt1[5] = options.stateInformation[1];
// 	 tmpOpt1[6] = options.stateInformation[2];
// 	 tmpOpt1[7] = options.omega_x.real();
// 	 tmpOpt1[8] = options.omega_y.real();
// 	 tmpOpt1[9] = options.omega_w.real();
// 	 tmpOpt1[10] = options.dispersion_x.real();
// 	 tmpOpt1[11] = options.dispersion_y.real();
// 	 tmpOpt1[12] = options.min_x;
// 	 tmpOpt1[13] = options.min_y;
// 	 tmpOpt1[14] = options.min_z;
// 	 tmpOpt1[15] = options.t_abs.real();
// 	 tmpOpt1[16] = options.exp_factor.real();
// 	 tmpOpt1[17] = options.g;
// 	 tmpOpt1[18] = options.ITP_step;
// 	 tmpOpt1[19] = options.RTE_step;


// 	for(int i = 0; i<4;i++)
// 	  tmpOpt1[i+20] = options.grid[i];

// 	tmpOpt1[24] = options.potFactor;
// 	tmpOpt1[25] = options.samplesize;
// 	tmpOpt1[26] = options.vortexnumber;
// 	tmpOpt1[27] = options.vortexspacing;
// 	tmpOpt1[28] = options.Ag;
// 	tmpOpt1[29] = options.OmegaG;    

// 	//copy options to file
// 	hsize_t dimsf[1] = {30};
// 	dataspace = H5Screate_simple(1, dimsf, NULL);

// 	h5a_options = H5Acreate(h5_timegroup, "Options", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
// 	H5Awrite (h5a_options, H5T_IEEE_F64LE, tmpOpt1);

// 	H5Aclose(h5a_options);
// 	H5Sclose(dataspace);

// 	hid_t h5a_meta;
// 	dimsf[0] = 9;
// 	dataspace = H5Screate_simple(1, dimsf, NULL);

// 	h5a_meta = H5Acreate(h5_timegroup, "Meta", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
// 	H5Awrite(h5a_meta, H5T_IEEE_F64LE, meta.data());
// 	H5Aclose(h5a_meta);
// 	H5Sclose(dataspace);
//   }

//   H5Gclose(h5_timegroup);

//   return true;
// }

// bool binaryFile::getSnapshot(const string &name, int snapShotTime, vector<ComplexGrid> &data,MatrixData::MetaData &meta, Options &options)
// {
//   if(m == out)
// 	{
// 		cout << "file "<< filename.c_str() << "is not in read mode" << endl;
// 		return false;
// 	}

//   stringstream time_name;
//   // snapShotTime will always be formatted with two decimal digits
//   time_name.setf(ios_base::fixed);
//   time_name.precision(2);
//   time_name << snapShotTime;

//   if(H5Lexists(h5_file, (time_name.str()).c_str(), H5P_DEFAULT)){
// 	h5_timegroup = H5Gopen(h5_file, (time_name.str()).c_str(), H5P_DEFAULT);
//   }else{
// 	cout << "ERROR: HDF5 Group for time " << snapShotTime << " does not exist" << endl;
// 	return false;
//   }

//   stringstream set_name;
//   set_name << name;

//   if(!H5Lexists(h5_timegroup, (set_name.str()).c_str(), H5P_DEFAULT)){
// 	cout << "ERROR: HDF5 Dataset " << name << " does not exist for time " << snapShotTime << " in ComplexGrid branch" << endl;

// 	H5Gclose(h5_timegroup);
// 	return false;
//   }else{
// 	hid_t test_id = H5Oopen(h5_timegroup, (set_name.str()).c_str(), H5P_DEFAULT);


	
// 	hid_t h5a_options;
// 	h5a_options = H5Aopen(h5_timegroup, "Options", H5P_DEFAULT);

// 	// int sizeOfOptions = (int)(H5Aget_storage_size(h5a_options)/sizeof(double));

// 	double tmpOpt1[30];

// 	H5Aread(h5a_options, H5T_IEEE_F64LE , tmpOpt1);

// 	//load Options struct from file array. Don't forget to change appropriately when changing Options struct
// 	options.N = tmpOpt1[0];
// 	for(int i= 0; i < 3; i++){
// 	 options.klength[i] = tmpOpt1[i+1];
// 	}
// 	options.stateInformation[0] = tmpOpt1[4];
// 	options.stateInformation[1] = tmpOpt1[5];
// 	options.stateInformation[2] = tmpOpt1[6];
// 	options.omega_x = complex<double>(tmpOpt1[7],0.0);
// 	options.omega_y = complex<double>(tmpOpt1[8],0.0);
// 	options.omega_w = complex<double>(tmpOpt1[9],0.0);
// 	options.dispersion_x = complex<double>(tmpOpt1[10],0.0);
// 	options.dispersion_y = complex<double>(tmpOpt1[11],0.0);
// 	options.min_x = tmpOpt1[12];
// 	options.min_y = tmpOpt1[13];
// 	options.min_z = tmpOpt1[14];
// 	options.t_abs = complex<double>(tmpOpt1[15],0.0);
// 	options.exp_factor = complex<double>(tmpOpt1[16],0.0);
// 	options.g = tmpOpt1[17];
// 	options.ITP_step = tmpOpt1[18];
// 	options.RTE_step = tmpOpt1[19];

// 	for(int i = 0; i<4;i++)
// 	  options.grid[i] = (uint32_t)tmpOpt1[i+20];

// 	options.potFactor = tmpOpt1[24];
// 	options.samplesize = (int)tmpOpt1[25];
// 	options.vortexnumber = (int)tmpOpt1[26];
// 	options.vortexspacing = (int)tmpOpt1[27];
// 	options.Ag = tmpOpt1[28];
// 	options.OmegaG = tmpOpt1[29];              

// 	H5Aclose(h5a_options);

// 	hid_t h5a_meta;
// 	h5a_meta = H5Aopen(h5_timegroup, "Meta", H5P_DEFAULT);
// 	H5Aread(h5a_meta, H5T_IEEE_F64LE, meta.data());
// 	H5Aclose(h5a_meta);
// 	meta.arrayToData();  
			  
// 	hsize_t vecsize;
// 	H5Gget_num_objs(test_id, &vecsize);
// 	data.resize(vecsize);

// 	for(int i = 0; i < vecsize; i++){
// 	  data[i] = ComplexGrid(options.grid[0],options.grid[1],options.grid[2],options.grid[3]);
// 	  stringstream comp;
// 	  comp << i;

// 	  hid_t dataset = H5Dopen(test_id, (comp.str()).c_str(), H5P_DEFAULT);
// 	  hid_t dataspace = H5Dget_space(dataset);

// 	  hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
// 	  hsize_t *dimf = (hsize_t *)malloc(rank*sizeof(hsize_t));
// 	  H5Sget_simple_extent_dims(dataspace, dimf, NULL);

// 	  H5Dread(dataset,  H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double *)data[i].get_address());

// 	  free(dimf);
// 	  H5Dclose(dataset);
// 	}
// 	H5Gclose(test_id);
//   }
//   H5Gclose(h5_timegroup);

//   return true;
// }