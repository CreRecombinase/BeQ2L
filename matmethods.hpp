#ifndef MATMETHODS_HPP
#define MATMETHODS_HPP
#include "mkl.h"
#include "mkl_vsl.h"
#include <algorithm>
#include "armadillo"

using namespace arma;
using namespace std;

mat CorMats(const mat &A, const mat &B){
  mat C = cor(A,B);
  C.elem(find_nonfinite(C)).zeros();
  return(C);
}

mat CVCorMAD(const mat &Est, const mat &TestA, const mat &TestB){
  mat testC = CorMats(TestA,TestB);
  mat Errmat = abs(Est-testC);
  return(Errmat);
}

umat GenBoot (size_t colsize, size_t bootstrapnumber){

  umat samplemat = randi<umat>(bootstrapnumber,colsize,distr_param(0,colsize-1));
 
  return(samplemat);
}

mat BootstrapSample(mat &A,umat &Bootmat,int i){

  return(A.rows(Bootmat.row(i)));
}

void Mdim(const mat &cmat,const string matname){
  cout<<matname<<" has "<<cmat.n_rows<<"rows and "<<cmat.n_cols<<" cols"<<endl;
}

mat BootCorMedianNS(mat &A, mat &B,mat &C,umat &BootMat,int bsi){
  
  int status;
  mat medians(A.n_cols,B.n_cols);
  if(A.n_rows!=BootMat.n_cols){
    cerr<<"Boot(incorrect indices):  A.n_rows= "<<A.n_rows<<" BM.n_cols= "<<BootMat.n_cols<<endl;
    throw 15;
  }
  int n= bsi;


  double params;

  
  cerr<<"Starting Bootstrap:"<<bsi<<endl;
  //First iteration of bootstrap
  mat tA(A.n_rows,A.n_cols);
  mat tB(B.n_rows,B.n_cols);
  mat tC(A.n_cols,B.n_cols);
  for(int i=0; i<bsi; ++i){
    tA = BootstrapSample(A,BootMat,i);
    tB = BootstrapSample(B,BootMat,i);
    
    tC = CorMats(tA,tB);
    C.col(i) = vectorise(tC,0);
  }
  medians = reshape(median(C,1),medians.n_rows,medians.n_cols);
  
  cerr<<"Bootstrap finished"<<endl;
  return(medians);
}






mat BootCorMedian(mat &A, mat &B,cube &C,cube &quants,umat &BootMat,int bsi){
  
  int status;
  mat medians(A.n_cols,B.n_cols);
  if(A.n_rows!=BootMat.n_cols){
    cerr<<"Boot(incorrect indices):  A.n_rows= "<<A.n_rows<<" BM.n_cols= "<<BootMat.n_cols<<endl;
    throw 15;
  }
  int n= bsi;
  int samplenum = C.n_slices;
  VSLSSTaskPtr task;
  MKL_INT q_order_n=1;
  double q_order[q_order_n];
  double params;
  MKL_INT p,nparams,xstorage;
  p = C.n_rows*C.n_cols;
  xstorage = VSL_SS_MATRIX_STORAGE_COLS;
  params = 0.001;
  q_order[0] = 0.5;
  
  status = vsldSSNewTask(&task,&p,&samplenum,&xstorage,C.memptr(),0,NULL);
  status = vsldSSEditStreamQuantiles(task,&q_order_n,q_order,quants.memptr(),&nparams,&params);
  
  cerr<<"Starting Bootstrap:"<<bsi<<endl;
  //First iteration of bootstrap
  mat tA(A.n_rows,A.n_cols);
  mat tB(B.n_rows,B.n_cols);
  int m=0;
  for(int i=0; i<bsi; ++i){
    tA = BootstrapSample(A,BootMat,i);
    tB = BootstrapSample(B,BootMat,i);
    
    m = i%C.n_slices;
    C.slice(m) = CorMats(tA,tB);
    
    if(m==(C.n_slices-1)){
      try{
	status = vsldSSCompute(task,VSL_SS_QUANTS,VSL_SS_METHOD_FAST);
      }
      catch(...){
	cerr<<"Something happened with computing the quantiles"<<endl;
      }
    }
  }
  samplenum=0;
  status = vsldSSCompute(task,VSL_SS_STREAM_QUANTS,VSL_SS_METHOD_SQUANTS_ZW);
  status = vslSSDeleteTask(&task);
  if(n%2==0){
    medians = 0.5*quants(span(0,0),span(),span())+0.5*quants(span(1,1),span(),span());
  }else{
    medians = quants(span(0,0),span(),span());
  }
  cerr<<"Bootstrap finished"<<endl;
  return(medians);
}



uvec trainindex(const size_t totalsize,const int ksize,const int k){
  uvec returnvec(totalsize-ksize);
  int j=-1;
  for(int i=0; i<totalsize; i++){
    if(i <(k*ksize)||i>((k+1)*ksize)-1){
      j++;
      returnvec(j)=i;
    }
  }
  if(j+1!=totalsize-ksize){
    cerr<<"j="<<j<<" returnsize= "<<totalsize-ksize-1<<endl;
  }
  return(returnvec);
}
  
uvec testindex(const size_t totalsize,const  int ksize,const  int k){
  uvec returnvec(ksize);
  int j=-1;
  for( int i=0; i<totalsize; i++){
    if(i>=(k*ksize)&&i<=((k+1)*ksize-1)){
      j++;
      returnvec(j)=i;

    }
  }
  if(j+1!=ksize){
    cerr<<"j="<<j<<" ksize= "<<ksize<<endl;
  }
  return(returnvec);
}


void KfoldCV (const mat &A,const mat &B, const int kfold, const int chunknum, const int bsi,const string outfilename,const int bsichunksize){


  int kfoldIterations = kfold;
  int iter_k = (int) floor((double)A.n_rows/(double)kfold);
  if(kfoldIterations<1){
    cerr<<"Must have at least one iteration of cross-validation, not "<<kfoldIterations<<endl;
    throw 12;
  }
  
  cout<<"Begin "<<iter_k<<"-fold cross validation ("<<kfoldIterations<<" iterations). There are "<<A.n_rows<<" samples"<<endl;
  uvec testi = testindex(A.n_rows,iter_k,0);
  uvec traini = trainindex(A.n_rows,iter_k,0);

  
  mat testA,testB,trainA,trainB;

  cout<<"Initiating variables."<<endl;
  mat Point(A.n_cols,B.n_cols);
  mat C(A.n_cols*B.n_cols,bsi);
  cout<<"Generating Bootmat (Dimensions should be "<<bsi<<"x"<<traini.n_elem<<")"<<endl;
  umat BootMat = GenBoot(traini.n_elem,bsi);
  cout<<"Starting kfold: 0"<<endl;
  mat MedianMat(A.n_cols,B.n_cols,fill::zeros);
  mat BootMAD(A.n_cols,B.n_cols,fill::zeros);
  mat pointMAD(A.n_cols,B.n_cols,fill::zeros);


  
  for(int i=0; i<kfoldIterations; i++){
    cout<<"Starting kfold: "<<i<<endl;
    testi = testindex(A.n_rows,iter_k,i);
    traini = trainindex(A.n_rows,iter_k,i);
    testA = A.rows(testi);
    testB = B.rows(testi);
    trainA = A.rows(traini);
    trainB = B.rows(traini);

    Point = CorMats(trainA,trainB);
 
    MedianMat = BootCorMedianNS(trainA,trainB,C,BootMat,bsi);
   
    BootMAD=BootMAD+CVCorMAD(MedianMat,testA,testB);
   
    pointMAD += CVCorMAD(Point,testA,testB);
  
  }
  pointMAD = pointMAD/kfoldIterations;
  BootMAD/=kfoldIterations;

  double pointSumMAD = accu(pointMAD)/pointMAD.n_elem;
  double SumMAD = accu(BootMAD)/BootMAD.n_elem;
  ofstream outputfile;

  outputfile.open(outfilename.c_str(),std::ofstream::out|std::ofstream::app);
  outputfile.setf(ios::fixed,ios::floatfield);

  outputfile.precision(10);  
  if(chunknum==0){
    outputfile<<"Chunk\tSize\tbsi\tPointMAD\tBootMAD"<<endl;
  }
  outputfile<<chunknum<<"\t"<<Point.n_elem<<"\t"<<bsi<<"\t"<<pointSumMAD<<"\t"<<SumMAD<<"\t"<<endl;
  
  outputfile.close();
}
    

#endif


