#ifndef H5OBJ_HPP
#define H5OBJ_HPP
#include "armadillo"
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdlib.h>
#include <openssl/md5.h>
#include <iostream>
#include <vector>
#include <string>
#include "hdf5.h"
#include "hdf5_hl.h"

using namespace std;
using namespace arma;

hid_t H5finst(const char* h5file,bool isTrunc){
  hid_t tempid;
  if(isTrunc){
    tempid = H5Fcreate(h5file,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    return(tempid);
  } else{ 
    tempid = H5Fopen(h5file,H5F_ACC_RDONLY,H5P_DEFAULT);
    return(tempid);
  }
}

int H5fclose(hid_t file_id){
  
  return(H5Fclose(file_id));
}

void readh5mat(hid_t file_id,const char*fieldname,size_t startcol,size_t colnum,size_t startrow,size_t rownum, Mat<double> &outmat){


  double* tempmat;
  herr_t status;
  hsize_t dimsout[2]={colnum,rownum};
  hsize_t offset[2]={startcol,startrow};
  hsize_t offset_out[2]={0,0};
  hsize_t dimsm[2]={colnum,rownum};
  hid_t matspace,matdata,memspace;
  tempmat = outmat.memptr();
  matdata = H5Dopen2(file_id,fieldname,H5P_DEFAULT);
  matspace = H5Dget_space(matdata);
  status = H5Sselect_hyperslab(matspace,H5S_SELECT_SET,offset,NULL,dimsm,NULL);
  memspace = H5Screate_simple(2,dimsm,NULL);
  status = status| H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,dimsm,NULL);
  status = status | H5Dread(matdata,H5T_NATIVE_DOUBLE,memspace,matspace,H5P_DEFAULT,tempmat);
  status = status | H5Dclose(matdata);
  status = status | H5Sclose(matspace);
  status = status | H5Sclose(memspace);

  if(status!=0){
    cerr<<"Error reading matrix!"<<endl;
    exit(255);
  }
}
vector<string> readcharhdf5(hid_t file_id,const char* fieldname,size_t startrow, size_t rownum){
  hsize_t rownamedims[1]={rownum};
  char ** temprownames;
  temprownames = (char**) malloc(rownum*sizeof(char*));
  hsize_t dimsout[1]={rownum};
  hsize_t offset[1]={startrow};
  hsize_t offset_out[1]={0};
  hsize_t dimsm[1]={rownum};
  hid_t rowset,rowspace,memspace;
  hid_t fstrtype;
  herr_t status;

  fstrtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(fstrtype,H5T_VARIABLE);
  rowset = H5Dopen2(file_id,fieldname,H5P_DEFAULT);
  rowspace = H5Dget_space(rowset);

  status = H5Sselect_hyperslab(rowspace,H5S_SELECT_SET,offset,NULL,dimsm,NULL);
  memspace = H5Screate_simple(1,dimsm,NULL);
  status = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,dimsm,NULL);
  status = H5Dread(rowset,fstrtype,memspace,rowspace,H5P_DEFAULT,temprownames);
  status = H5Dclose(rowset);
  status = H5Sclose(rowspace);
  status = H5Sclose(memspace);
  vector<string> outputstring(temprownames,temprownames+rownum);
  free(temprownames);
  return(outputstring);
}
herr_t writecharparamhd5(hid_t file_id, const char *fieldname, const char* value){

  hid_t root,dataspace,att,fstrtype;
  hsize_t dim[1]={1};
  herr_t ret;
  root = H5Gopen(file_id,"/",H5P_DEFAULT);
  dataspace = H5Screate_simple(1,dim,NULL);
  fstrtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(fstrtype,H5T_VARIABLE);
  att = H5Acreate(root,fieldname,fstrtype,dataspace,H5P_DEFAULT,H5P_DEFAULT);
  ret = H5Awrite(att,fstrtype,&value);
  ret = ret| H5Sclose(dataspace);
  ret = ret|H5Sclose(att);
  
  return(ret);
}
  

herr_t writemathdf5(hid_t file_id, const char* fieldname, mat matrix){

  inplace_trans(matrix);
  hsize_t matdims[2];
  hid_t matdataspace;
  hid_t matdataset;
  herr_t status;
  double *matptr = matrix.memptr();
  

  matdims[0]=matrix.n_cols;
  matdims[1]=matrix.n_rows;

  matdataspace = H5Screate_simple(2,matdims,NULL);
  matdataset = H5Dcreate(file_id,fieldname,H5T_NATIVE_DOUBLE,matdataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status = H5Dwrite(matdataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,matptr);
  status = status | H5Sclose(matdataspace);
  status = status | H5Dclose(matdataset);
  return(status);
}  
  

herr_t writecharhdf5(hid_t file_id,const char* fieldname,const vector<string> rowdata){
  hsize_t rownamedims[1];
  hid_t fstrtype;
  hid_t rownamespace,rownamedata;
  herr_t status;
  char** rawarray;

  rownamedims[0]=rowdata.size();
  rawarray = (char**) malloc(rowdata.size()*sizeof(char*));

  for(int i=0; i<rowdata.size(); i++)
    {
      rawarray[i] = (char*)malloc(rowdata[i].size()+1*sizeof(char));
      strcpy(rawarray[i],rowdata[i].c_str());
    }

  fstrtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(fstrtype,H5T_VARIABLE);
  rownamespace = H5Screate_simple(1,rownamedims,NULL);
  rownamedata = H5Dcreate(file_id,fieldname,fstrtype,rownamespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status = H5Dwrite(rownamedata,fstrtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,rawarray);
  status = status|H5Dclose(rownamedata);
  status = status|H5Sclose(rownamespace);
  for(int i=0; i<rowdata.size(); i++){
    free(rawarray[i]);
  }
  free(rawarray);
  return(status);
}




unsigned long get_size_by_fd(int fd){
  struct stat statbuf;
  if(fstat(fd,&statbuf)<0) exit(-1);
  return(statbuf.st_size);
}

string hashfile(const char* filename){
  int file_descript;
  unsigned char* result;
  char tstr[3];
  result = (unsigned char*) malloc(MD5_DIGEST_LENGTH*sizeof(unsigned char));
  unsigned long file_size;
  char* file_buffer;
  file_descript = open(filename,O_RDONLY);
  if(file_descript<0){
    cerr<<"Bad file descriptor!"<<endl;
    exit(-1);
  }
  file_size = get_size_by_fd(file_descript);

  file_buffer = (char*)mmap(0,file_size,PROT_READ,MAP_SHARED,file_descript,0);

  MD5((unsigned char*) file_buffer, file_size,result);
  munmap(file_buffer,file_size);
  string strresult;
  for( int i=0; i<MD5_DIGEST_LENGTH; i++){
    sprintf(tstr,"%02x",result[i]);
    strresult = strresult+string(tstr);
  }
  return(strresult);
}

#endif
