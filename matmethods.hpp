#ifndef MATMETHODS_HPP
#define MATMETHODS_HPP
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <algorithm>
#include "armadillo"

using namespace arma;
using namespace std;



double RMSE(const mat &Est, const mat &TrueA, const mat &TrueB){
  uvec zeros;
  double res;
  mat Aimp = Est*TrueB.t();
  vec  rsq(TrueA.n_rows);
  for(int i=0; i<TrueA.n_rows; i++){
    rsq(i)=as_scalar(cor(Aimp.row(i),TrueA.col(i)));
  }
  zeros =  find_nonfinite(rsq);
  rsq.elem(find_nonfinite(rsq)).zeros();
  rsq = square(rsq);
  if(zeros.n_elem!=rsq.n_elem){
    res = sum(rsq)/(rsq.n_elem-zeros.n_elem);
  }
  else{
    res = 0;
  }
  return(res);
}

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



double CumMeanBoot(mat &A, mat &B,mat &C,umat &BootMat,int bsi,mat &TrueA,mat &TrueB){

  double res;
  running_stat_vec<vec> stats;
  mat means(A.n_cols,B.n_cols);
  if(A.n_rows!=BootMat.n_cols){
    cerr<<"Boot(incorrect indices):  A.n_rows= "<<A.n_rows<<" BM.n_cols= "<<BootMat.n_cols<<endl;
    throw 15;
  }
  
  cerr<<"Starting Bootstrap:"<<endl;
  //First iteration of bootstrap
  mat tA = BootstrapSample(A,BootMat,0);
  mat tB = BootstrapSample(B,BootMat,0);
  
  vec Bv = var(tB,1,0);
  C = cov(tA,tB);
  C.each_row() /= Bv;
  C.cols(find(Bv==0)).zeros();
  uvec stillbad = find_nonfinite(C);
  if(stillbad.n_elem>0){
    cerr<<"Error with Logic(Boot0):"<<stillbad.n_elem<<endl;
    uvec badb = find(Bv==0);
    cerr<<Bv.n_elem<<endl;
    cerr<<C.n_cols<<"x"<<C.n_rows<<endl;
    cerr<<badb.n_elem<<endl;
    throw 12;
    
  }
  stats(vectorise(C));
  
  for(int i=1; i<bsi; ++i){
    //      cout<<"Boot: "<<i<<endl;
    try{
      tA = BootstrapSample(A,BootMat,i);
      tB = BootstrapSample(B,BootMat,i);
    }
    catch(...){
      cerr<<"Can't create bootstraps!"<<endl;
      throw 2;
    }
    try{
      Bv = var(tB,1,0);
      C = cov(tA,tB);
      C.each_row()/=Bv;
      C.cols(find(Bv==0)).zeros();
      stillbad = find_nonfinite(C);
      if(stillbad.n_elem>0){
	cerr<<"Error with Logic(Boot1):"<<stillbad.n_elem<<endl;
	uvec badb = find(Bv==0);
	cerr<<Bv<<endl;
	cerr<<badb.n_elem<<endl;
	throw 11;
      }
    }
    catch(...){
      cerr<<"Error computing correlation"<<endl;
      throw 3;
    }
    
    
  }
  
  cerr<<"Bootstrap finished"<<endl;
  stats(vectorise(C));
  means = stats.mean();	
  means.reshape(C.n_rows,C.n_cols);
  res = RMSE(means,TrueA,TrueB);
  
  return(res);
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


void KfoldCV (const mat &A,const mat &B, const int kfold, const int chunknum, const int bsi,const string file){

  mat testA,testB,trainA,trainB;
  int kfoldIterations = kfold;
  int iter_k = (int) floor((double)A.n_rows/(double)kfold);
  if(kfoldIterations<1){
    cerr<<"Must have at least one iteration of cross-validation, not "<<kfoldIterations<<endl;
    throw 12;
  }
  
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
  cout<<"Initiating variables."<<endl;
  mat Point(A.n_cols,B.n_cols);
  mat C(A.n_cols,B.n_cols);

  umat BootMat = GenBoot(trainA.n_rows,bsi);
  cout<<"Starting kfold: 0"<<endl;
  double Rsq = CumMeanBoot(trainA,trainB,C,BootMat,bsi,testA,testB);

 
  Point = cov(trainA,trainB);
  vec Bv = var(trainB,1,0);
  Point.each_row()/=Bv;
  Point.cols(find(Bv==0)).zeros();
  uvec stillbad = find_nonfinite(Point);
  if(stillbad.n_elem>0){
    cerr<<"Logic error(Point0):"<<stillbad.n_elem<<endl;
    uvec badb = find(Bv==0);
    cerr<<Bv<<endl;
    cerr<<badb.n_elem<<endl;
    throw 15;
  }
  
  double  pointSumRMSE = RMSE(Point,testA,testB);
  
  for(int i=1; i<kfoldIterations; i++){
    cout<<"Starting kfold: "<<i<<endl;
    testi = testindex(A.n_rows,iter_k,i);
    traini = trainindex(A.n_rows,iter_k,i);
    testA = A.rows(testi);
    testB = B.rows(testi);
    trainA = A.rows(traini);
    trainB = B.rows(traini);
    Bv = var(trainB,1,0);
    Point = cov(trainA,trainB);
    Point.each_row()/=Bv;
    Point.cols(find(Bv==0)).zeros();
    stillbad = find_nonfinite(C);
    if(stillbad.n_elem>0){
      cerr<<"Logic error(Point1):"<<stillbad.n_elem<<endl;
      uvec badb = find(Bv==0);
      cerr<<Bv<<endl;
      cerr<<badb.n_elem<<endl;
      throw 14;
      
    }
    Rsq =Rsq+CumMeanBoot(trainA,trainB,C,BootMat,bsi,testA,testB);
    pointSumRMSE += RMSE(Point,testA,testB);
  }
  pointSumRMSE = pointSumRMSE/kfoldIterations;
  Rsq /=kfoldIterations;
  
  ofstream outputfile;
  outputfile.open(file,std::ofstream::out|std::ofstream::app);
  outputfile.setf(ios::fixed,ios::floatfield);
  outputfile.precision(10);  
  if(chunknum==0){
    outputfile<<"Chunk\tSize\tbsi\tPointMeanR2\tBootMeanR2"<<endl;
  }
  outputfile<<chunknum<<"\t"<<Point.n_elem<<"\t"<<bsi<<"\t"<<pointSumRMSE<<"\t"<<Rsq<<"\t"<<endl;
  
  outputfile.close();
}
    
    
  

  
  


#endif


