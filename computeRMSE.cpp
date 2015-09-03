#include <iostream>
#include "matmethods.hpp"
#include "h5objcpp.hpp"

using namespace std;

int main(int argc, char* argv[]){
  
  char* h5file;
  char* outputfile;
  H5File file;
  int chunkstart, chunknum,Achunksize,Bchunksize;
  int totalAchunks,totalBchunks,totalchunks;
  int kfold,bsi;
  mkl_set_num_threads(1);
  int rownum;
  if(argc!=9){
    cerr<<"Usage:RMSE h5file chunkstart chunknum Achunksize Bchunksize bsi kfold outputfile"<<endl;
    cerr<<"argc: "<<argc<<endl;
    return(2);
  }
  h5file = argv[1];
  chunkstart = atoi(argv[2]);
  chunknum = atoi(argv[3]);
  Achunksize = atoi(argv[4]);
  Bchunksize = atoi(argv[5]);
  bsi = atoi(argv[6]);
  kfold = atoi(argv[7]);
  outputfile = argv[8];
  
  

  H5std_string Amatfield = "matA";
  H5std_string Bmatfield= "matB";

  try{
    file = H5File(h5file,H5F_ACC_RDONLY);
  } 
  catch(...)
    {
      cerr<<"Can't open file: "<<h5file<<endl;
      throw 3;
    }

    
  
  vector <hsize_t> sizevec(2);
  sizevec = getdims(file,Amatfield);
  rownum = sizevec[1];

  totalAchunks = ceil((double) sizevec[0]/(double) Achunksize);

  sizevec = getdims(file,Bmatfield);
  totalBchunks = ceil((double) sizevec[0]/(double) Bchunksize);
  totalchunks = totalAchunks*totalBchunks;
  cout<<"Total Chunks="<<totalchunks<<endl;
  
  int chunkstop = min(chunkstart+chunknum,totalchunks)+1;
  
  for( int i=chunkstart; i<chunkstop; i++){
    if(i>=totalchunks){
      return(0);
    }
    int Achunk = i/totalBchunks;
    int Bchunk = i%totalBchunks;

    mat A=readh5mat(file,Amatfield,Achunk,Achunksize,0,rownum);
    mat B=readh5mat(file,Bmatfield,Bchunk,Bchunksize,0,rownum);
    cout<<"Starting chunk: "<<i<<endl;
    KfoldCV(A,B,kfold,i,bsi,outputfile,10);
  }
  return(0);
}

    



  
