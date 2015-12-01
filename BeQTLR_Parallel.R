#Script to Run BeQTLR on test data on hpc cluster
library(BeQTLR)
library(rhdf5)
args <- commandArgs(trailingOnly=T)
if(length(args)!=9){
  print(length(args))
  stop("Usage: Rscript BeQTLR_Parallel.R hdf5file chunkstart chunkstop Achunksize Bchunksize BootstrapNumber TestSetSize KfoldCount Outfile")
}
h5file <- args[1]
chunkstart <- as.integer(args[2])
chunkstop <- as.integer(args[3])
Achunksize <- as.integer(args[4])
Bchunksize <- as.integer(args[5])
bsi <- as.integer(args[6])
testsize <- as.integer(args[7])
kfoldcount <- as.integer(args[8])
outfile <- args[9]


#The first step is to read in the hdf5 SNP and expresion data
h5getdim <- function(h5file,spacename){
  fid <- H5Fopen(h5file)
  did <- H5Dopen(fid,spacename)
  hid <- H5Dget_space(did)
  sid <- H5Sget_simple_extent_dims(hid)$size
  H5Dclose(did)
  H5Fclose(fid)
  return(sid)
}
genchunks <- function(chunk,Achunksize,Afeatures,Bchunksize,Bfeatures){
  Achunktotal <- ceiling(Afeatures/Achunksize)
  Bchunktotal <- ceiling(Bfeatures/Bchunksize)
  chunktotal <- Achunktotal*Bchunktotal
  chunkmat <- matrix(1:chunktotal,Achunktotal,Bchunktotal)
  Achunk <- row(chunkmat)[chunk]
  Bchunk <- col(chunkmat)[chunk]
  seq(from=0,to=Afeatures,by=Achunksize)
  diff(seq(from=0,to=Bfeatures,by=Bchunksize))
  AllA <- 1:Afeatures
  AllB <- 1:Bfeatures
  Aind <- ((Achunk-1)*Achunksize):(min(c(Achunk*Achunksize,Afeatures))-1)+1
  Bind <- ((Bchunk-1)*Bchunksize):(min(c(Bchunk*Bchunksize,Bfeatures))-1)+1
return(list(Achunkind=Aind,Bchunkind=Bind))
}
TestTrain <- function(testsize,kcount,A,B) {
  testind <- replicate(kcount,sample(0:(nrow(A)-1),testsize,replace = F))
  trainind <- apply(testind,2,function(x) (0:(nrow(A)-1))[-(x+1)])
  testl <- apply(testind,2,MatSplit,A=A,B=B)
  trainl <- apply(trainind,2,MatSplit,A=A,B=B)
  return(list(testl,trainl))
}

Adim <- h5getdim(h5file,"EXP")
Bdim <- h5getdim(h5file,"SNP")
if(Adim[1]!=Bdim[1]){
  stop(paste0("Number of cases differs between Amat(",Adim[1],") and Bmat(",Bdim[1],")"))
}


#First we'll generate the boostrap matrix
for(i in chunkstart:chunkstop) {
  pointmadout <- paste0( outfile, "_",i,"_point_mad.txt")
  pointrmseout<- paste0( outfile, "_",i,"_point_rmse.txt")
  bootmadout <- paste0(  outfile, "_",i,"_boot_mad.txt")
  bootrmseout <- paste0( outfile,"_",i,"_boot_rmse.txt")

  print(paste0("Chunk ",i," of ",chunkstop))
  chunkl <- genchunks(i,Achunksize,Adim[2],Bchunksize,Bdim[2])
  matA <- h5read(file = h5file,name="EXP",index=list(1:(Adim[1]),chunkl$Achunkind))
  Anames <- h5read(file=h5file,name="EXPsymbol",index=list(chunkl$Achunkind))
  Achrom <- h5read(file=h5file,name="EXPchrom",index=list(chunkl$Achunkind))
  Astart <- h5read(file=h5file,name="EXPstart",index=list(chunkl$Achunkind))
  Astop <- h5read(file=h5file,name="EXPstop",index=list(chunkl$Achunkind))
  matB <- h5read(file=h5file,name="SNP",index=list(1:(Bdim[1]),chunkl$Bchunkind))
  Bnames <- h5read(file=h5file,name="SNPrsid",index=list(chunkl$Bchunkind))
  Bchrom <- h5read(file=h5file,name="SNPchrom",index=list(chunkl$Bchunkind))
  Bpos <- h5read(file=h5file,name="SNPpos",index=list(chunkl$Bchunkind))
  SNPdf <- data.frame(rsid=Bnames,Chrom=Bchrom,Pos=Bpos,stringsAsFactors = F)
  Genedf <- data.frame(Symbol=Anames,Chrom=Achrom,Start=Astart,Stop=Astop,stringsAsFactors = F)
  fullboot <- GenBoot(samplesize=nrow(matA),bootstrapnumber=bsi)
  print("Running Overall Bootstrap")
  fullbootmeds <- BeQTL(matA,matB,fullboot)
  fullPoint <- PointCor(matA,matB)
  ttl <- TestTrain(testsize,kfoldcount,matA,matB)
  testl <- ttl[[1]]
  trainl <- ttl[[2]]

  Bootmat <- GenBoot(samplesize = nrow(trainl[[1]][["A"]]),bootstrapnumber = bsi)

  print("Running kfold bootstrap")
  bootk <- sapply(trainl,function(x){
   BeQTL(x[["A"]],x[["B"]],Bootmat)
  })
  print("Generating Point Estimates")
  pointk <- sapply(trainl,function(x){
    PointCor(x[["A"]],x[["B"]])
  })
  print("Generating Test Estimates")
  testk <- sapply(testl,function(x){
    PointCor(x[["A"]],x[["B"]])
  })
  print("Computing MAD")
  bootmad <- MAD(bootk,testk,T,c(ncol(matA),ncol(matB)))
  pointmad <- MAD(pointk,testk,T,c(ncol(matA),ncol(matB)))
  print("Computing RMSE")
  bootrmse <- RMSE(bootk,testk,c(ncol(matA),ncol(matB)))
  pointrmse <- RMSE(pointk,testk,c(ncol(matA),ncol(matB)))

  bootresmad <- SumRes(fullbootmeds,bootmad,SNPdf,Genedf,nrow(matA),qt(0.95,df = 171))
  bootresrmse <- SumRes(fullbootmeds,bootrmse,SNPdf,Genedf,nrow(matA),qt(0.95,df = 171))
  pointresmad <- SumRes(fullPoint,pointmad,SNPdf,Genedf,nrow(matA),qt(0.95,df=171))
  pointresrmse <- SumRes(fullPoint,pointrmse,SNPdf,Genedf,nrow(matA),qt(0.95,df=171))
  print ("Writing Results")

  write.table(bootresmad,bootmadout,col.names=T,row.names=F,sep="\t",quote=F)
  write.table(bootresrmse,bootrmseout,col.names=T,row.names=F,sep="\t",quote=F)

  write.table(pointresmad,pointmadout,col.names=T,row.names=F,sep="\t",quote=F)
  write.table(pointresrmse,pointrmseout,col.names=T,row.names=F,sep="\t",quote=F)

}





