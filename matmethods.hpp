#ifndef MATMETHODS_HPP
#define MATMETHODS_HPP
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include "armadillo"

using namespace arma;
using namespace std;



umat GenBoot (size_t colsize, size_t bootstrapnumber){
  boost::random::mt19937 rng;

  boost::random::uniform_int_distribution<> samples(0,colsize-1);

  umat samplemat(bootstrapnumber,colsize);
  
  for(umat::iterator it=samplemat.begin(); it!=samplemat.end();++it){
    (*it) = samples(rng);
  }
  return(samplemat);
}

mat BootstrapSample(mat &A,umat &Bootmat,int i){
  return(A.rows(Bootmat.row(i)));
}

void CumMeanBoot(mat &A, mat &B,mat &C,umat &BootMat, mat &means, mat &vars,int bsi){


  running_stat_vec<vec> stats;


  
  if(A.n_rows!=BootMat.n_cols){
    cerr<<"Boot(incorrect indices):  A.n_rows= "<<A.n_rows<<" BM.n_cols= "<<BootMat.n_cols<<endl;
    throw 15;
  }
  else{

    //First iteration of bootstrap
    mat tA = BootstrapSample(A,BootMat,0);
    mat tB = BootstrapSample(B,BootMat,0);
    
    C = cor(tA,tB);
    stats(vectorise(C));
    
    for(int i=1; i<bsi; ++i){
      //      cout<<"Boot: "<<i<<endl;
      tA = BootstrapSample(A,BootMat,i);
      tB = BootstrapSample(B,BootMat,i);
      C = cor(tA,tB);
      stats(vectorise(C));
    }
  }

  means = stats.mean();
  means.reshape(C.n_rows,C.n_cols);
  vars = stats.var();
  vars.reshape(C.n_rows,C.n_cols);
 
}
  

void doBoot(mat &A, mat &B,cube &C,umat &BootMat, cube &Summaries){
  
  if(A.n_rows!=BootMat.n_cols){
    for(int i=1; i<C.n_slices; ++i){
      cerr<<"Boot(incorrect indices): "<<i<<endl;
      throw 15;
      //      C.slice(i) = cor(A.rows(BootMat(i,span(0,A.n_rows))),B.rows(BootMat(i,span(0,B.n_rows))));
    }
  }
  else{

    //First iteration of bootstrap
    mat tA = BootstrapSample(A,BootMat,0);
    mat tB = BootstrapSample(B,BootMat,0);

    C.slice(0) = cor(tA,tB);

    for(int i=1; i<C.n_slices; ++i){
      cout<<"Boot: "<<i<<endl;
      tA = BootstrapSample(A,BootMat,i);
      tB = BootstrapSample(B,BootMat,i);
      C.slice(i) = cor(tA,tB);
    }
  }

  for(int i=0; i<C.n_slices; i++){
    mat S = mat(C.n_rows,C.n_slices);
    for(int j=0; j< C.n_cols; j++){
      S = C(span::all,span(i,i),span(0,i));
      Summaries.slice(i).col(j)=median(S,1);
    }
  }
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

double RMSE(mat &Est, const mat &True){
  Est = sqrt(pow((Est-True),2));
  Est.elem(find_nonfinite(Est)).zeros();
  return(mean(mean(Est)));
}
    
  



void KfoldCV (const mat &A,const mat &B, const int kfold, const int chunknum, const int bsi,const string file){

  mat testA,testB,trainA,trainB;
  int kfoldIterations = (int) floor((double)A.n_rows/(double)kfold);
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
  mat means(A.n_cols,B.n_cols);
  mat vars(A.n_cols,B.n_cols);
  mat C(A.n_cols,B.n_cols);
  
  umat BootMat = GenBoot(trainA.n_rows,bsi);
  CumMeanBoot(trainA,trainB,C,BootMat,means,vars,bsi);
  Point = cor(trainA,trainB);
  mat TrueCor = cor(testA,testB);
  double pointSumRMSE = RMSE(Point,TrueCor);
  double bootSumRMSE = RMSE(means,TrueCor);
  
  for(int i=1; i<kfoldIterations; i++){
    cout<<"Starting kfold: "<<i<<endl;
    testi = testindex(A.n_rows,iter_k,i);
    traini = trainindex(A.n_rows,iter_k,i);
    testA = A.rows(testi);
    testB = B.rows(testi);
    trainA = A.rows(traini);
    trainB = B.rows(traini);
    Point = cor(trainA,trainB);
    CumMeanBoot(trainA,trainB,C,BootMat,means,vars,bsi);
    TrueCor = cor(testA,testB);
    pointSumRMSE += RMSE(Point,TrueCor);
    bootSumRMSE += RMSE(means,TrueCor);
  }
  pointSumRMSE = pointSumRMSE/kfoldIterations;
  bootSumRMSE = bootSumRMSE/kfoldIterations;
  
  ofstream outputfile;
  outputfile.open(file,std::ofstream::out|std::ofstream::app);
  
  if(chunknum==0){
    outputfile<<"Chunk\tSize\tbsi\tPointSumRMSE\tBootSumRMSE"<<endl;
  }
  outputfile<<chunknum<<"\t"<<Point.n_elem<<"\t"<<"\t"<<pointSumRMSE<<"\t"<<bootSumRMSE<<endl;
  
  outputfile.close();
}
    
    
    


    
    
    
    

    
  




  
  

  


//void ComputeQuantiles(const running_stat_vec<rowvec> &stats, const mat &Co, mat &Quant, const vector<double> quantiles){
  

  
  


#endif


