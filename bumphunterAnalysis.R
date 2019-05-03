library(bumphunter)
load("Covariates.RData")
imposeFilter <- function(vec, range1, range2, thres1, thres2){
  Na1 <- length(which(vec[range1]=="NaN"))
  Na2 <- length(which(vec[range2]=="NaN"))
  if (Na1 <= thres1 & Na2 <= thres2) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}
covarX <- as.matrix(data.frame(intercept = 1, disease = covar$Disease[covar$Male==0], age = covar$Age[covar$Male==0]))
GenomeSSCres <- list()
library(doParallel)
for (i in 1:22) {
  FName <- paste0("MethylationBetaValues_chr",i,".RData")
  load(FName)
  tmp <- datMeth[apply(datMeth, 1, function(x) imposeFilter(x,c(1,2,12:18),c(4,6,7,9),3,1)),]
  print(dim(tmp))
  
  pos <- unlist(strsplit(rownames(tmp),"-"))[seq(2,length(unlist(strsplit(rownames(tmp),"-"))),2)]
  chr <- rep(paste0("chr",i), length(pos))
  pos <- as.numeric(pos)
  cl <- clusterMaker(chr, pos, maxGap = 200)
  print(length(unique(cl)))
  meth <- as.matrix(tmp)
  meth <- meth[,covar$Male == 0]
  registerDoParallel(cores = 16)
  SSCres <- bumphunter(meth,covarX,chr,pos,cl,cutoff=0.2,coef=2,nullMethod="bootstrap",B=250)
  GenomeSSCres[[i]] <- SSCres$table
}

load("MethylationBetaValues_chrX.RData")
tmp <- datMeth[apply(datMeth, 1, function(x) imposeFilter(x,c(1,2,12:18),c(4,6,7,9),3,1)),]
print(dim(tmp))

pos <- unlist(strsplit(rownames(tmp),"-"))[seq(2,length(unlist(strsplit(rownames(tmp),"-"))),2)]
chr <- rep(paste0("chr",i), length(pos))
pos <- as.numeric(pos)
cl <- clusterMaker(chr, pos, maxGap = 200)
print(length(unique(cl)))
meth <- as.matrix(tmp)
meth <- meth[,covar$Male == 0]
registerDoParallel(cores = 16)
SSCres <- bumphunter(meth,covarX,chr,pos,cl,cutoff=0.2,coef=2,nullMethod="bootstrap",B=250)
GenomeSSCres[[23]] <- SSCres$table