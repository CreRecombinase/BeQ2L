

int main(int argc, char* argv[]){
  
  char* paramfile;
  int chunkstart,chunknum;
  
  fparam param(paramfile,chunkstart,chunknum);
  
  Mat<int> BootMat = GenBoot(param.rownum,param.bootnum);
  mat A(param.rownum,param.Acols(1));
  mat BA(param.rownum,param.Acols(1));
  mat B(param.rownum,param.Bcols(1));
  mat BB(param.rownum,param.Bcols(1));
  mat Co(param.Acols(1),param.Bcols(1));
  if(param.streaming){
    running_stat_vec<rowvec> stats;
  } else{ 
  cube C(param.Acols(1),param.Bcols(1),param.bootnum);
  }
  cube Quant(param.Acols(1).param.Bcols(1),param.quantiles.size());


  
  for( int i=param.chunkstart; i<param.chunknum; i++){
    if((param.Acols(i)!=A.n_cols)||(param.Bcols(i)!=B.n_cols)){
      A.set_size(param.rownum,param.Acols(i));
      B.set_size(param.rownum,param.Bcols(i));
      if(param.streaming){
	stats.set_size(param.Acols(i),param.Bcols(i));
      } else {
	C.set_size(param.Acols(i),param.Bcols(i),param.bootnum);

      }
      Quant.set_size(param.Acols(i),param.Bcols(i),param.quantiles.size());
      Co.set_size(param.Acols(i),param.Bcols(i));
      BA.set_size(param.rownum,param.Acols(i));
      BB.set_size(param.rownum,param.Bcols(i));
    }

    
    readh5mat(param.h5id,param.Amatfield,param.startA(i),param.Acols(i),0,param.rownum,A);
    readh5mat(param.h5id,param.Bmatfield,param.startB(i),param.Bcols(i),0,param.rownum,B);
    Co=cor(A,B);
    for(in j=0; j<param.bootnum; j++){
      BA=A.rows(BootMat.col(j));
      BB=B.rows(BootMat.col(j));
      if( param.streaming){
	Co=cor(BA,BB);
	stats(vectorise(cor(BA,BB),1));
      }else{
	C.slice(j)=cor(BA,BB);
      }
    }

    if(param.streaming){
      ComputeQuantiles(stats,Co,Quant);
    }else{
      ComputeQuantiles(C,Co,&Quant);
    }
    filterQuantiles(param,Quant);

    writeQuantiles(param,Quant);
    
    return(0);
    
  }

  
    
    
    
