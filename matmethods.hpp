#ifndef MATMETHODS_HPP
#define MATMETHOIDS_HPP
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include "armadillo"

Mat<int> GenBoot (size_t colsize, size_t bootstrapnumber){
  boost::random::mt19937 rng;

  boost::random::uniform_int_distribution<> samples(0,colsize-1);

  Mat<int> samplemat(bootstrapnumber,colsize);
  
  for(Mat<int>::iterator it=samplemat.begin(); it!=samplemat.end();++it){
    (*it) = samples(rng);
  }
  return(samplemat);
}


void ComputeQuantiles(const running_stat_vec<rowvec> &stats, const mat &Co, mat &Quant){
  





