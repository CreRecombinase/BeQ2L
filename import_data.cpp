#include "armadillo"
#include "importdata.hpp"
#include <iostream>
using namespace std;    
using namespace arma;

int main(int argc, char* argv[])
{
  //Data comes in with one column per case, but is then transposed
  if(argc<6){
    cout<<"usage "<<argv[0]<<" matrixA.txt Arows.txt matrixB.txt Brows.txt output.h5"<<endl;
    return(-1);
  }
  hid_t file_id;
  char* matAfile,*matBfile,*rowAfile,*rowBfile,*outputfile;
  matAfile=argv[1];
  rowAfile = argv[2];
  matBfile = argv[3];
  rowBfile = argv[4];
  outputfile = argv[5];
  
  vector<string> rowAnames;
  vector<string> rowBnames;
  mat A;
  mat B;
  int Arows,Acols;
  int Brows,Bcols;

  string hashrowA = hashfile(rowAfile);
  string hashrowB = hashfile(rowBfile);
  string hashmatA = hashfile(matAfile);
  string hashmatB = hashfile(matBfile);
  cout<<"Hash of file "<<rowAfile<<" is "<<hashrowA<<endl;
  
  A.load(matAfile,raw_ascii);
  rowAnames = readnames(rowAfile);
  rowBnames = readnames(rowBfile);
  inplace_trans(A);
  Arows=A.n_rows;
  Acols=A.n_cols;
  B.load(matBfile,raw_ascii);
  inplace_trans(B);
  Brows=B.n_rows;
  Bcols=B.n_cols;

  file_id = H5finst(outputfile,true);
  A.print();
  cout<<"Matrix A has "<<A.n_rows<<" rows and "<<A.n_cols<<" columns"<<endl;
  writemathdf5(file_id,"matA",A);
  writemathdf5(file_id,"matB",B);
  writecharhdf5(file_id,"Acolnames",rowAnames);
  writecharhdf5(file_id,"Bcolnames",rowBnames);
  writecharparamhd5(file_id,"Anames_md5",hashrowA.c_str());
  writecharparamhd5(file_id,"Bnames_md5",hashrowB.c_str());
  writecharparamhd5(file_id,"Amat_md5",hashmatA.c_str());
  writecharparamhd5(file_id,"Bmat_md5",hashmatB.c_str());
  
  H5fclose(file_id);
  
  return(0);
}
    
