#ifndef MATMETHODS_HPP
#define MATMETHODS_HPP
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/timer/timer.hpp>
#include "armadillo"

using namespace arma;
using namespace std;
using boost::timer::auto_cpu_timer;


umat GenBoot (size_t colsize, size_t bootstrapnumber){
  boost::random::mt19937 rng;

  boost::random::uniform_int_distribution<> samples(0,colsize-1);

  umat samplemat(bootstrapnumber,colsize);
  
  for(umat::iterator it=samplemat.begin(); it!=samplemat.end();++it){
    (*it) = samples(rng);
  }
  return(samplemat);
}

void doBoot(mat &A, mat &B,cube &C,umat &BootMat, mat &Summaries){
  auto_cpu_timer timerInd;
  timerInd.stop();
  auto_cpu_timer timerCor;

  timerCor.stop();
  auto_cpu_timer timerMedian;
  timerMedian.stop();
  
  if(A.n_rows!=BootMat.n_cols){

    for(int i=1; i<C.n_slices; ++i){
      cout<<"Boot(incorrect indices): "<<i<<endl;

      C.slice(i) = cor(A.rows(BootMat(i,span(0,A.n_rows))),B.rows(BootMat(i,span(0,B.n_rows))));
    }
  }
  else{
    
    timerInd.resume();
    mat tA = A.rows(BootMat.row(0));
    mat tB = B.rows(BootMat.row(0));
    timerInd.stop();

    timerCor.resume();
    C.slice(0) = cor(tA,tB);
    timerCor.stop();
    for(int i=0; i<C.n_slices; ++i){
      cout<<"Boot: "<<i<<endl;
      timerInd.resume();
      tA = A.rows(BootMat.row(i));
      tB = B.rows(BootMat.row(i));
      timerInd.stop();
      timerCor.resume();
      C.slice(i) = cor(tA,tB);
      timerCor.stop();
    }
  }
  timerMedian.resume();
  mat S = mat(C.n_rows,C.n_slices);
  for(int i=0; i< C.n_cols; i++){
    S = C(span::all,span(i,i),span::all);
    Summaries.col(i)=median(S,1);
  }
  timerMedian.stop();
  timerInd.report();
  timerCor.report();
  timerMedian.report();
  
}

uvec trainindex(const size_t totalsize,const int ksize,const int k){
  uvec returnvec(totalsize-ksize);
  int j=0;
  for(int i=0; i<totalsize; i++){
    if(i <(k*ksize)||i>((k+1)*ksize)-1){
      returnvec(j)=i;
      j++;
    }
  }
  return(returnvec);
}
  

uvec testindex(const size_t totalsize,const  int ksize,const  int k){
  uvec returnvec(ksize);
  int j=0;
  for( int i=0; i<totalsize; i++){
    if(i>=(k*ksize)&&i<=((k+1)*ksize-1)){
      returnvec(j)=i;
      j++;
    }
  }
  return(returnvec);
}

double RMSE(mat &Est, const mat &True){
  Est = sqrt(pow((Est-True),2));
  Est.elem(find_nonfinite(Est)).zeros();
  return(mean(mean(Est)));
}
    
  



void KfoldCV (const mat &A,const mat &B, const int kfold, const int chunknum, const int bsi,const string file){

  mat testA,testB,trainA,trainB;
  int kfoldIterations = (int) ceil((double)A.n_rows/(double)kfold);
  if(kfoldIterations<1){
    cerr<<"Must have at least one iteration of cross-validation, not "<<kfoldIterations<<endl;
    throw 12;
  }
  
  int iter_k = kfold;


  cout<<"Begin "<<iter_k<<"-fold cross validation. There are "<<A.n_rows<<" samples"<<endl;
  uvec testi = testindex(A.n_rows,iter_k,0);
  uvec traini = trainindex(A.n_rows,iter_k,0);

  try{
    testA = A.rows(testi);
    testB = B.rows(testi);

    trainA = A.rows(traini);
    trainB = B.rows(traini);
  } catch(...){
    cerr<<"Error indexing training and testing sets:"<<endl;
    cerr<<"Training: "<<traini<<endl;
    cerr<<"Testing: "<<testi<<endl;
    throw 4;
  }
  cout<<"Initiating variables and starting kfold: 0"<<endl;
  mat Point(A.n_cols,B.n_cols);
  mat Summaries(A.n_cols,B.n_cols);
  umat BootMat = GenBoot(trainA.n_rows,bsi);
  cube C(A.n_cols,B.n_cols,bsi);
  doBoot(trainA,trainB,C,BootMat,Summaries);
  Point = cor(trainA,trainB);
  mat TrueCor = cor(testA,testB);
  double pointSumRMSE = RMSE(Point,TrueCor);
  double bootSumRMSE = RMSE(Summaries,TrueCor);
  for(int i=1; i<kfoldIterations; i++){
    cout<<"Starting kfold: "<<i<<endl;
    testi = testindex(A.n_rows,iter_k,0);
    traini = trainindex(A.n_rows,iter_k,0);
    testA = A.rows(testi);
    testB = B.rows(testi);
    trainA = A.rows(traini);
    trainB = B.rows(traini);
    doBoot(trainA,trainB,C,BootMat,Summaries);
    Point = cor(trainA,trainB);
    TrueCor = cor(testA,testB);
    pointSumRMSE += RMSE(Point,TrueCor);
    bootSumRMSE += RMSE(Summaries,TrueCor);
  }
  pointSumRMSE = pointSumRMSE/kfoldIterations;
  bootSumRMSE = bootSumRMSE/kfoldIterations;

  ofstream outputfile;
  outputfile.open(file,std::ofstream::out|std::ofstream::app);
  
  if(chunknum==0){
    outputfile<<"Chunk\tSize\tbsi\tPointSumRMSE\tBootSumRMSE"<<endl;
  }
  outputfile<<chunknum<<"\t"<<Point.n_elem<<"\t"<<bsi<<"\t"<<pointSumRMSE<<"\t"<<bootSumRMSE<<endl;
  outputfile.close();
}
    
    
    


    
    
    
    

    
  




  
  

  


//void ComputeQuantiles(const running_stat_vec<rowvec> &stats, const mat &Co, mat &Quant, const vector<double> quantiles){
  

  
  


#endif


