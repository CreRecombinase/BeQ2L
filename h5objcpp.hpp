#ifndef H5OBJ_HPP
#define H5OBJ_HPP
#include "armadillo"
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <string>
#include "H5Cpp.h"

using namespace std;
using namespace arma;
using namespace H5;

//Functions to facilitate the reading and writing of data to and from HDF5

  
int readh5mat(H5File file,const H5std_string fieldname,size_t colchunk, size_t colchunksize, size_t rowchunk, size_t rowchunksize, Mat<double> &outmat){
  //This function takes an OPEN hdf5 file, a field name, a start column and start row, and the address of an already allocated armadillo matrix and reads to the matrix

  //Data is stored in row major order in HDF5, and column major order in Armadillo, so we have to store the column data as rows

  hsize_t startcol = colchunk*colchunksize;
  hsize_t startrow = rowchunk*rowchunksize;
  
  hsize_t colnum = colchunksize;
  hsize_t rownum = rowchunksize;
  
  double* tempmat;
  hsize_t total_dims[2];                     //We find the size of the HDF5 dataset to check that what we're asking for isn't larger than the whole dataset
  hsize_t dims_out[2]={rownum,colnum};      //We want a data structure that is rowsizexcolumnsize in memory
  hsize_t offset_h5[2]={startcol,startrow}; //Our starting rows and columns in the HDF5 matrix
  hsize_t offset_out[2]={0,0};              //No offset in memory
  hsize_t dims_h5[2]={colnum,rownum};       //The HDF5 data is stored as columns and then rows
  tempmat = outmat.memptr();                //The address of our matrix

  try{
    DataSet dataset = file.openDataSet(fieldname);
    DataSpace dataspace = dataset.getSpace();
    dataspace.getSimpleExtentDims(total_dims,NULL);
    int changsize=0;

    if(colchunksize*colchunk+colchunksize>total_dims[0]){

      dims_out[1]=total_dims[0]-(colchunksize*colchunk);
      dims_h5[0]=total_dims[0]-(colchunksize*colchunk);
      outmat.set_size(rowchunksize,dims_out[1]);
    }
    if(rowchunksize*rowchunk+rowchunksize>total_dims[1]){

      dims_out[0] = total_dims[1]-(rowchunksize*rowchunk);
      dims_h5[1] = total_dims[1]-(rowchunksize*rowchunk);
      outmat.set_size(dims_out[0],colchunksize);
    }
    dataspace.selectHyperslab(H5S_SELECT_SET,dims_h5,offset_h5);  //Ask for HDF5 slab of size and offset we specified above
    
    DataSpace memspace(2,dims_out);                               
    memspace.selectHyperslab(H5S_SELECT_SET,dims_out,offset_out); //Specify the memory we're going to read in to 
    
    dataset.read(tempmat,PredType::NATIVE_DOUBLE,memspace,dataspace);  //Read the data

    dataspace.close();
    memspace.close();
    dataset.close();   //Close everything up
    
  }
  catch(FileIException error)
    {
      error.printError();
      cerr<<"FileIException error"<<endl;
      int tendp = offset_h5[0]+dims_h5[0];
      int endp[2];
      endp[0] = tendp;
      tendp = offset_h5[1]+dims_h5[1];
      endp[1] = tendp;
      
      cerr<<"Wanted to read from "<<offset_h5[0]<<"to "<<endp[0]<<"(Dim1) And from "<<offset_h5[1]<<" to "<<endp[1]<<" (Dim2)"<<endl;
      cerr<<"dims_out:"<<dims_out[0]<<","<<dims_out[1]<<endl;
      cerr<<"offset_h5:"<<offset_h5[0]<<","<<offset_h5[1]<<endl;
      cerr<<"dims_h5:"<<dims_h5[0]<<","<<dims_h5[1]<<endl;
      throw(1);
    }
  catch (DataSetIException error)
    {
      error.printError();
      cerr<<"DataSetIException error"<<endl;
      int tendp = offset_h5[0]+dims_h5[0];
      int endp[2];
      endp[0] = tendp;
      tendp = offset_h5[1]+dims_h5[1];
      endp[1] = tendp;
      cerr<<"Wanted to read from "<<offset_h5[0]<<"to "<<endp[0]<<"(Dim1) And from "<<offset_h5[1]<<" to "<<endp[1]<<" (Dim2)"<<endl;
      cerr<<"Fieldname:"<<fieldname<<endl;
      cerr<<"dims_out:"<<dims_out[0]<<","<<dims_out[1]<<endl;
      cerr<<"offset_h5:"<<offset_h5[0]<<","<<offset_h5[1]<<endl;
      cerr<<"dims_h5:"<<dims_h5[0]<<","<<dims_h5[1]<<endl;
      throw(1);
    }
  catch (DataSpaceIException error)
    {
      cerr<<"DataSpaceIException error"<<endl;
      int tendp = offset_h5[0]+dims_h5[0];
      int endp[2];
      endp[0] = tendp;
      tendp = offset_h5[1]+dims_h5[1];
      endp[1] = tendp;
      cerr<<"Wanted to read from "<<offset_h5[0]<<"to "<<endp[0]<<"(Dim1) And from "<<offset_h5[1]<<" to "<<endp[1]<<" (Dim2)"<<endl;
      cerr<<"Fieldname:"<<fieldname<<endl;
      cerr<<"dims_out:"<<dims_out[0]<<","<<dims_out[1]<<endl;
      cerr<<"offset_h5:"<<offset_h5[0]<<","<<offset_h5[1]<<endl;
      cerr<<"dims_h5:"<<dims_h5[0]<<","<<dims_h5[1]<<endl;
      error.printError();
      throw(1);
    }
  return(0);
}

vector<hsize_t> getdims(H5File file, H5std_string datasetname){

  vector<hsize_t> retvec(2);
  DataSet dataset = file.openDataSet(datasetname);
  DataSpace dataspace = dataset.getSpace();
  dataspace.getSimpleExtentDims(&retvec[0],NULL);
  dataset.close();
  dataspace.close();
  return(retvec);
}



int writemathdf5(H5File file, const H5std_string fieldname, mat &matrix){
  //This function takes an OPEN HDF5 file, a field, and the address of a prepopulated matrix and writes it to the opened file

  //It is important to transpose the data before writing it, so that it will read in the correct dimensions.  This is due to the array ordering differences between HDF5 and Armadillo

  hsize_t matdims[2];
  double *matptr = matrix.memptr();
  

  matdims[0]=matrix.n_cols;
  matdims[1]=matrix.n_rows;
  
  try{
    
    DataSpace filespace(2,matdims);
    
    DataSet dataset = file.createDataSet(fieldname,PredType::NATIVE_DOUBLE,filespace);
    dataset.write(matptr,PredType::NATIVE_DOUBLE);
    filespace.close();
    dataset.close();
  }

  catch(FileIException error)
    {
      cerr<<"Error writing file:"<<endl;
      error.printError();
      throw(2);
    }
  catch(DataSetIException error)
    {
      cerr<<"Error writing file:"<<endl;
      error.printError();
      throw(2);
    }
  catch(DataSpaceIException error)
    {
      cerr<<"Error writing file:"<<endl;
      error.printError();
      throw(2);
    }
  return(0);
}  


  
vector<string> readcharhdf5(H5File file,const H5std_string fieldname,size_t chunk, size_t chunksize){


  hsize_t chunkstart = chunksize*chunk;
  hsize_t rownamedims[1]={chunksize};
  char ** temprownames;
  temprownames = (char**) malloc(chunksize*sizeof(char*));
  hsize_t dims_out[1]={chunksize};
  hsize_t offset_h5[1]={chunkstart};
  hsize_t offset_out[1]={0};
  hsize_t dims_h5[1]={chunksize};
  hid_t rowset,rowspace,memspace;
  hid_t fstrtype;
  herr_t status;

  try{
    
    Exception::dontPrint();
    StrType stringt(0,H5T_VARIABLE);
    
    DataSet dataset = file.openDataSet(fieldname);
    DataSpace dataspace = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET,dims_h5,offset_h5);
    DataSpace memspace(1,dims_h5,NULL);
    memspace.selectHyperslab(H5S_SELECT_SET,dims_h5,offset_out);
    dataset.read(temprownames,stringt,memspace,dataspace);
  }
  catch(FileIException error)
    {
      cerr<<"Error reading rownames: "<<fieldname<<endl;
      error.printError();
      throw(3);
    }
  catch (DataSetIException error)
    {
      cerr<<"Error reading rownames: "<<fieldname<<endl;
      error.printError();
      throw(3);

    }
  catch (DataSpaceIException error)
    {
      cerr<<"Error reading rownames: "<<fieldname<<endl;
      error.printError();
      throw(3);
    }
  vector<string> outputstring(temprownames,temprownames+chunksize); 
  free(temprownames);
  return(outputstring);
}
int writecharparamhd5(H5File file, const char* datasetnamei,const char* attributenamei, const string valuei){
  hsize_t dim[1]={1};

  try{
    Exception::dontPrint();
    H5std_string datasetname(datasetnamei);
    H5std_string attributename(attributenamei);
    H5std_string value(valuei);
    
 
    StrType stringt(0,H5T_VARIABLE);
    DataSet dataset = file.openDataSet(datasetname);
    DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
    Attribute myatt_in = dataset.createAttribute(attributename,stringt,attr_dataspace);
    myatt_in.write(stringt,value);
    myatt_in.close();
    attr_dataspace.close();
    dataset.close();
    stringt.close();

  }

  catch(FileIException error)
    {
      cerr<<"Error writing parameter: "<<attributenamei<<endl;
      error.printError();
      throw(4);
    }
  catch (DataSetIException error)
    {
      cerr<<"Error writing parameter: "<<attributenamei<<endl;
      error.printError();
      throw(4);

    }
  catch (DataSpaceIException error)
    {
      cerr<<"Error writing parameter: "<<attributenamei<<endl;
      error.printError();
      throw(4);
    }
  return(0);
}


int writecharhdf5(H5File file,const H5std_string fieldname,const vector<string> matnames){
  
  hsize_t matnamedims[1];
  hid_t fstrtype;
  char** rawarray;

  matnamedims[0]=matnames.size()-1;
  rawarray = (char**) malloc(matnames.size()*sizeof(char*));

  for(int i=0; i<matnames.size(); i++)
    {
      rawarray[i] = (char*)malloc(matnames[i].size()+1*sizeof(char));
      strcpy(rawarray[i],matnames[i].c_str());
    }
  
  try{
     Exception::dontPrint();
    DataSpace dataspace(1,matnamedims);
    StrType stringt(PredType::C_S1,H5T_VARIABLE);
    
    DataSet dataset = file.createDataSet(fieldname,stringt,dataspace);
    
    dataset.write(rawarray,stringt);
    dataset.close();
    dataspace.close();
    stringt.close();
  }
  catch(FileIException error)
    {
      cerr<<"Error writing rownames(File error) "<<fieldname<<endl;
      error.printError();
      throw(5);
    }

  catch(DataSpaceIException error)
    {
      cerr<<"Error writing rownames(Dataspace error) "<<fieldname<<endl;
      error.printError();
      throw(5);
    }
  catch(DataSetIException error)
    {
      cerr<<"Error writing rownames(DataSetError) "<<fieldname<<endl;
      error.printError();
      throw(5);
    }

  for(int i=0; i<matnames.size(); i++){
    free(rawarray[i]);
  }
  free(rawarray);
  return(0);
}

  



#endif
