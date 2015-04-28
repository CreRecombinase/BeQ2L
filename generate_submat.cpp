#include "armadillo"
#include "importdata.hpp"
#include "matmethods.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <utility>
#include <iterator>
#include "H5Cpp.h"
#include "h5objcpp.hpp"
using namespace std;
using namespace arma;
using namespace H5;

int main(int argc, char* argv[]){

  char* h5infile;
  char* snpfilename;
  char* genefilename;
  char* h5outfile;
  int chunksize;


  H5File infile;
  H5File outfile;


  h5infile = argv[1];
  snpfilename = argv[2];
  genefilename = argv[3];
  chunksize = atoi(argv[4]);
  h5outfile = argv[5];
  H5std_string Genefield = "matA";
  H5std_string Genelistfield = "Acolnames";
  H5std_string Snpfield = "matB";
  H5std_string Snplistfield = "Bcolnames";

  infile = H5File(h5infile,H5F_ACC_RDONLY);
  outfile = H5File(h5outfile,H5F_ACC_TRUNC);


  cout<<"Reading in genelist"<<endl;
  ifstream genefile(genefilename);
  vector<string> genelist;
  copy(istream_iterator<string>(genefile),
	    istream_iterator<string>(),
	    back_inserter(genelist));
  cout<<"Sorting genelist"<<endl;
  sort(genelist.begin(),genelist.end());

  vector <hsize_t> sizevec(2);
  sizevec = getdims(infile,Genefield);
  int rownum = sizevec[1];

  int Achunks = ceil((double) sizevec[0]/(double) chunksize);
    
  vector<string> intergenes;
  mat intergenemat;
  cout<<"Starting reading and subsetting of data"<<endl;
  for(int i=0; i<Achunks; i++){
    cout<<"Reading genelist from h5file"<<endl;
    vector<string> filegenes = readcharhdf5(infile,Genelistfield,i,chunksize);
    unordered_map<string,int> geneindex;
    pair<string,int> tp;
    cout<<"Making index of genes"<<endl;
    for(int j = 0; j!=filegenes.size(); j++){
      tp = make_pair(filegenes[j],j);
      geneindex.insert(tp);
    }
    cout<<"Sorting genelist"<<endl;
    sort(filegenes.begin(),filegenes.end());
    cout<<"Computing Intersect"<<endl;
    vector<string> intersectv(filegenes.end()-filegenes.begin());
    vector<string>::iterator intsct = set_intersection(filegenes.begin(),filegenes.end(),genelist.begin(),genelist.end(),intersectv.begin());
    intersectv.resize(intsct-intersectv.begin());
    uvec intersectIndices(intersectv.size());
    cout<<"Generating index for intersect"<<endl;
    for(vector<string>::iterator it=intersectv.begin(); it !=intersectv.end(); ++it){
      intersectIndices[it-intersectv.begin()]=geneindex[*it];
    }
    cout<<"Sorting index"<<endl;
    sort(intersectIndices.begin(),intersectIndices.end());
    
    if(intersectv.size()>0){
      intergenes.insert( intergenes.end(), intersectv.begin(), intersectv.end() );
      mat genemat = readh5mat(infile,Genefield,i,chunksize,0,rownum);
      if(intergenemat.n_rows>0){
	mat sgenemat = genemat.rows(intersectIndices);
	intergenemat.insert_rows(intergenemat.n_rows,sgenemat);
      }
      else{
	intergenemat= genemat.rows(intersectIndices);
      }
    }
  }
  cout<<intergenemat.n_rows<<"x"<<intergenemat.n_cols<<endl;
  cout<<intergenes.size()<<endl;
  
  
  return(0);
}








//  sizevec = getdims(h5infile,"matB");
// int Bchunks = ceil((double) sizevec[0]/(double) chunksize);
