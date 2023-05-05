##################################################################################
### Simple example with three populations having 10 SNPs of different minor allele frequencies
### Each population has 20 subjects
##################################################################################


source('SparsePCA_ftn.r')

n1<-n2<-n3<-20; n<-n1+n2+n3; d<-100
grp.idx<-rep(1:3,c(n1,n2,n3))

p1<-c(0.041, 0.650, 0.040, 0.083, 0.625, 0.017, 0.050, 0.600, 0.110, 0.045) # First 10 SNPs' MAF for group 1
p2<-c(0.666, 0.052, 0.080, 0.050, 0.041, 0.083, 0.050, 0.050, 0.610, 0.083) # First 10 SNPs' MAF for group 2
p3<-c(0.045, 0.033, 0.690, 0.667, 0.045, 0.600, 0.600, 0.050, 0.020, 0.625) # First 10 SNPs' MAF for group 3
p<-runif(d-10,min=0,max=0.5) # d-10 SNPs' MAF for all groups

grp1<-matrix(0,nr=n1,nc=10)
grp2<-matrix(0,nr=n2,nc=10)
grp3<-matrix(0,nr=n3,nc=10)

for( j in 1:10 ){
	grp1[,j]<-rowSums(matrix(rbinom(n1*2,1,p1[j]),nc=2,byrow=TRUE))
	grp2[,j]<-rowSums(matrix(rbinom(n2*2,1,p2[j]),nc=2,byrow=TRUE))
	grp3[,j]<-rowSums(matrix(rbinom(n3*2,1,p3[j]),nc=2,byrow=TRUE))
}
X.sim <- rbind( grp1, grp2, grp3 )

X.ran <- matrix(0, nr=n, nc=d-10)
for( j in 1:(d-10) ) X.ran[,j]<-rowSums(matrix(rbinom(n*2,1,p[j]),nc=2,byrow=TRUE))

X<-cbind(X.sim,X.ran)
X<-t(t(X)-colMeans(X))

## standard PCA
sing<-svd(X)
standard.A<-X%*%sing$v; colnames(standard.A)<-paste('PCscore',1:ncol(standard.A),sep='') # PC scores
standard.V<-sing$v; colnames(standard.V)<-paste('PC',1:ncol(standard.V),sep='') # PC
plot(as.data.frame(standard.A[,1:5]),col=grp.idx)
par(mfrow=c(2,1))
plot(standard.V[,1],type='h',xlab='SNP index',ylab='PC1')
plot(standard.V[,2],type='h',xlab='SNP index',ylab='PC2')


## sparse PCA with lasso or adaptive lasso penalty
method<-'alasso' # adaptive lasso penalty. For lasso, use method<-'lasso'
lambda<-2^c(-Inf,(-8):4)
k<-2	# total number of PCs to be estimated

tX<-X
Sparse.U<-NULL; Sparse.V<-NULL; Sparse.s<-NULL; BIC<-NULL; Sparse.nzr<-NULL; Sparse.lam<-NULL
for(a in 1:k){
  cat('The ', a, 'th principal component is being extracted....\n', sep='')
  res<-S.PCA(tX,lambda,method, prt=TRUE)
  min.bic.idx<-which(res$bic==min(res$bic))
  if (length(min.bic.idx)>1) min.bic.idx<-min.bic.idx[1]
  lambda.1<-seq(lambda[min.bic.idx-1],lambda[min.bic.idx+1],length=20)
  res.1<-S.PCA(tX,lambda.1,method,prt=TRUE)
	Sparse.U<-cbind(Sparse.U, res.1$u)
	Sparse.V<-cbind(Sparse.V, res.1$v)
	Sparse.s<-append(Sparse.s, res.1$s)
	BIC<-cbind(BIC, res.1$bic)
	Sparse.nzr<-cbind(Sparse.nzr, res.1$nzr)
	Sparse.lam<-cbind(Sparse.lam, res.1$lambda)
	tX <- tX-res.1$s*outer(res.1$u,res.1$v)
}
Sparse.A<-t(t(Sparse.U)*Sparse.s); colnames(Sparse.A)<-paste('PCscore',1:ncol(Sparse.A),sep='') # PC scores
colnames(Sparse.V)<-paste('PC',1:ncol(Sparse.V),sep='')
plot(as.data.frame(Sparse.A[,1:k]),col=grp.idx)
par(mfrow=c(2,1))
plot(Sparse.V[,1],type='h',xlab='SNP index',ylab='PC1')
plot(Sparse.V[,2],type='h',xlab='SNP index',ylab='PC2')
