rq <- function(data,scale='sd')
{
x<-data
#standardize the data default: by standard deviation and square root of n
x_standard<-scale(x,center=TRUE,scale=TRUE)/sqrt(nrow(x))
#no standardisation; use original data
 if (scale == 'data') { 
x_standard<-as.matrix(x)
}

#standardize by standard deviation and square root of n
if (scale == 'sd') { 
x_standard<-scale(x,center=TRUE,scale=TRUE)/sqrt(nrow(x))

}
#standardize by standard deviation
if (scale == 'normal') { 
x_standard<-scale(x,center=TRUE,scale=TRUE)

}

#standardize by centring
if (scale == 'centre'||scale == 'center'||scale=="coordinate") { 
x_standard<-scale(x,center=TRUE,scale=FALSE)

}


#standardize by square root of n
if (scale =='n') { 
x_standard<-scale(x,center=TRUE,scale=FALSE)/sqrt(nrow(x))
}

#develop correlation matrix
corr_matrix<-cor(x_standard)

if(scale=="data"||scale=="n"||scale=="centre"||scale=="center"
||scale=="coordinate"||scale=="mds"){
#corr_matrix<-t(x_standard)%*%x_standard
corr_matrix<-crossprod(x_standard)
}

if(scale=="coordinate"){
#corr_matrix<-t(x_standard)%*%x_standard
#corr_matrix<-x_standard %o% x_standard
}


#develop eigen vector and values
eigen<-eigen(corr_matrix)
eigen_vector<-eigen$vec
eigen_val<-eigen$val

#diagonal matrix
diagonal_matrix<- diag(sqrt(eigen_val))

#develop R mode loadings
R_loading<-eigen_vector%*%diagonal_matrix
colnames(R_loading) <- colnames(R_loading, do.NULL = FALSE, prefix = "Factor")
#develop Q-mode loadings
Q_mode_loading<-x_standard%*%eigen_vector
colnames(Q_mode_loading) <- colnames(Q_mode_loading, do.NULL = FALSE, prefix = "Factor")
#define one axis for all loadings
all_loadings=rbind(R_loading,Q_mode_loading)

#compute scores
rscores<-x_standard%*%R_loading
rownames(rscores) <- rownames(data)
#qscores<-corr_matrix%*%eigen_vector
qscores<-crossprod(x_standard,Q_mode_loading)
combined.scores<-rbind(rscores,qscores)

#develop row and column names for the loadings
colnames=colnames(x)
nrow=nrow(x)
total_obs=length(colnames)+nrow
listnames=""
i=1
#variables names
variables=""
varsymbols=""
varsymbol=""
while (i <= total_obs) {
   if (i<=(length(colnames))) {
listnames=paste(listnames,colnames[i],sep=",")
variables=paste(variables,colnames[i],sep=",")
}
else {
number=i-length(colnames)
obs_start=NULL
if(!is.null(obs_start))
{
if(is.numeric(obs_start))
{
year=(number-1)+obs_start
listnames=paste(listnames,year,sep=",")

}
else{
year=(number-1)+obs_start
listnames=paste(listnames,year,sep=",obs_start")

}

}
else{
listnames=paste(listnames,number,sep=",")
}

}
i=i+1
 }
rownames=listnames
rownames=strsplit(listnames,",")
rownames1=rownames[[1]]
rownames=rownames1[2:length(rownames1)]
variables=strsplit(variables,',')
variables=variables[[1]]
variables=variables[2:length(variables)]
rownames(R_loading) <- c(variables)
all_loadings=rbind(R_loading,Q_mode_loading)
if(!is.null(obs_start)){
rownames(all_loadings) <- c(rownames)
#if(is.numeric(obs_start))
#{
#total_rows=obs_start:length(data)
#rownames(all_loadings) <- c(rownames)

#}
}

rownames(Q_mode_loading) <- rownames(Q_mode_loading, do.NULL = FALSE, prefix = "")
#check for validity of the model
#scale.matrix<-Q_mode_loading%*%diag((sqrt(eigen_val)^(-0.5)))%*%t(R_loading)
pca=eigen_vector
colnames(pca) <- colnames(pca, do.NULL = FALSE, prefix = "PC")
rownames(pca) <- c(variables)
pca.scores=eigen_vector*data
rownames(pca.scores) <- rownames(data)


list(correlation=corr_matrix,eigen.vector=eigen_vector,eigen.value=eigen_val,
diagonal.matrix=diagonal_matrix,r.loading=R_loading,q.loading=Q_mode_loading,
rloading=R_loading,qloading=Q_mode_loading,rloadings=R_loading,qloadings=Q_mode_loading,
combined.loadings=all_loadings,r.scores=rscores,q.scores=qscores,rscores=rscores,qscores=qscores,combined.scores=combined.scores,data=x,
rownames=rownames,variables=variables,mds=combined.scores,coordinates=Q_mode_loading,x.standard=x_standard,loadings=all_loadings,scores=combined.scores,
pca.loadings=pca,pca=pca,pca.scores=pca.scores,pcascores=pca.scores)
}

#generic function
qrfactor<-function(data,scale="sd") UseMethod ("qrfactor")

#default function
qrfactor.default<-function(data,scale="sd")
{
#convert file data to matrix
x<-as.matrix(data)

factor<-rq(data,scale)

factor$call<-match.call()

class(factor)<-"qrfactor"
factor
}

#summary function
summary.qrfactor<-function(object,...)
{
x<-object
cat("Call:\n")
print(x$call)
cat("\nCall: correlation matrix\n")
print(x$correlation)

cat("\n eigen  value\n")
print(x$eigen.value)

cat("\n R-loadings\n")
print(x$r.loading)

cat("\n Q-loadings\n")
print(x$q.loading)


}

#print function
print.qrfactor<-function(x,...)
{
cat("Call:\n")
print(x$call)

cat("\n Eigen values\n")
print(x$eigen.value)

cat("\nR-mode loadings\n")
print(x$r.loading)

cat("\n Q-mode loadings\n")
print(x$q.loading)

cat("\n Combine loadings; R mode loadings first\n")
print(x$combined.loading)

cat("\n R-mode scores\n")
print(x$r.scores)

cat("\n Q-mode scores\n")
print(x$q.scores)

cat("\n combined scores\n")
print(x$combined.scores)

}

plot.qrfactor<-function(x,factors=c(1,2),type="loading",plot="",...)
{
dev.new()
graphics.off()
par(cex="0.7",cex.lab="0.9")

xlab=paste("Factor ",factors[1])
ylab=paste("Factor ",factors[2])
if(type=="mds"){
xlab=paste("MDS axis ",factors[1])
ylab=paste("MDS axis ",factors[2])
main="Multidimensional Scaling"

}
if(type=="pca2"||type=="coordinate"||type=="coord"){
xlab=paste("Principal coordinate ",factors[1])
ylab=paste("Principal coordinate ",factors[2])
main="Principal Coordinate Analysis"
}

if(type=="ca"||type=="correspondence"){
xlab=paste("Correspondence axis ",factors[1])
ylab=paste("Correspondence axis ",factors[2])
main="Correspondence Analysis"

}


if(type=="scores"|type=="score"||type=="mds"){
#plot(x$combined.scores[,factors[1]],x$combined.scores[,factors[2]],xlab=xlab,ylab=ylab, main="R- and Q- Mode FA")
# text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=4)  
 #text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$r.score), cex=0.9, pos=4, col="red")  

#plot(x$q.score[,factors[1]],x$q.score[,factors[2]],xlab=xlab,ylab=ylab,main="Scores:R- and Q- Mode FA",col="red")
# text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=2, col="red")  
#points(x$r.score[,factors[1]],x$r.score[,factors[2]],pch=2, main="R- and Q- Mode FA")
 #text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$r.score), cex=0.9, pos=2)  

plot(x$scores[,factors[1]],x$scores[,factors[2]],xlab=xlab,ylab=ylab,main="Scores:R- and Q- Mode FA",col="black")
 text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=2, col="red")  
#points(x$r.score[,factors[1]],x$r.score[,factors[2]],pch=2, main="R- and Q- Mode FA")
 text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$data), cex=0.9, pos=2, col="blue")  

if(plot=="r"){
plot(x$r.score[,factors[1]],x$r.score[,factors[2]],xlab=xlab,ylab=ylab, main="Scores:R-Mode FA")
 text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$data), cex=0.9, pos=4, col="red")  
}
if(plot=="q"){
plot(x$q.score[,factors[1]],x$q.score[,factors[2]],xlab=xlab,ylab=ylab,main="Scores:Q- Mode FA",col="black")
text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=4, col="black")  
}
if(plot=="qr"|plot=="rq"){
par(mfrow=c(1,2)) 
plot(x$r.score[,factors[1]],x$r.score[,factors[2]],xlab=xlab,ylab=ylab, main="Scores:R-Mode FA",col="red")
 text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$data), cex=0.9, pos=4, col="red")  
plot(x$q.score[,factors[1]],x$q.score[,factors[2]],xlab=xlab,ylab=ylab,main="Q- Mode FA")
 text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=4, col="black")  
}

if(plot=="all"){
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(x$scores[,factors[1]],x$scores[,factors[2]],xlab=xlab,ylab=ylab,main="Scores:R- and Q- Mode FA",col="black")
 text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=2, col="red")  
#points(x$r.score[,factors[1]],x$r.score[,factors[2]],pch=2, main="R- and Q- Mode FA")
 text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$data), cex=0.9, pos=2, col="blue")  
plot(x$r.score[,factors[1]],x$r.score[,factors[2]],xlab=xlab,ylab=ylab, main="Scores:R-Mode FA")
 text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$data), cex=0.9, pos=4, col="blue")  
plot(x$q.score[,factors[1]],x$q.score[,factors[2]],xlab=xlab,ylab=ylab,main="Scores:Q- Mode FA")
 text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=4, col="red")  
}

}
else{
if(type=="pcaloadings"||type=="pca"||type=="eigen.vectors"||type=="eigenvectors"){
plot(x$pca.loadings[,factors[1]],x$pca.loadings[,factors[2]],xlab=xlab,ylab=ylab, main="PCA Loadings")
 text(x$pca.loadings[,factors[1]],x$pca.loadings[,factors[2]], row.names( x$pca.loadings), cex=0.9, pos=4)  

}else{
plot(x$loadings[,factors[1]],x$loadings[,factors[2]],xlab=xlab,ylab=ylab,main="Loadings:R- and Q- Mode FA") 
#xlim=c((min(x$r.loading[,factors[1]])-(0.2*max(x$r.loading[,factors[1]]))):(max(x$r.loading[,factors[1]])+(0.2*max(x$r.loading[,factors[1]])))),
#ylim=c(min(x$r.loading[,factors[2]]):max(x$r.loading[,factors[2]])),xaxs="r",yaxs="r",las=1)
 text(x$r.loading[,factors[1]],x$r.loading[,factors[2]], row.names( x$r.loading), cex=0.9, pos=3, col="red")  
#points(x$q.loading[,factors[1]],x$q.loading[,factors[2]],pch=2, main="R- and Q- Mode FA")
text(x$q.loading[,factors[1]],x$q.loading[,factors[2]], row.names( x$q.loading), cex=0.9, pos=2,col="blue")  
#xlim=c((min(x$q.loading[,factors[1]])-(0.2*max(x$q.loading[,factors[1]]))):(max(x$q.loading[,factors[1]])+(0.2*max(x$q.loading[,factors[1]])))),
#ylim=c(min(x$q.loading[,factors[2]]):max(x$q.loading[,factors[2]])),xaxs="r",yaxs="r",las=1)

if(plot=="r"){
plot(x$r.loading[,factors[1]],x$r.loading[,factors[2]],xlab=xlab,ylab=ylab, main="Loadings:R-Mode FA",col="red")
 text(x$r.loading[,factors[1]],x$r.loading[,factors[2]], row.names( x$r.loading), cex=0.9, pos=4, col="red")  
}
if(plot=="q"||type=="coordinate"||plot=="coordinate"||type=="coord"||plot=="coord"){
plot(x$q.loading[,factors[1]],x$q.loading[,factors[2]],xlab=xlab,ylab=ylab,main="Q- Mode FA",col="black")
text(x$q.loading[,factors[1]],x$q.loading[,factors[2]], row.names( x$q.loading), cex=0.9, pos=4, col="black")  
}
if(plot=="qr"|plot=="rq"){
par(mfrow=c(1,2)) 
plot(x$r.loading[,factors[1]],x$r.loading[,factors[2]],xlab=xlab,ylab=ylab, main="R-Mode FA",col="red")
 text(x$r.loading[,factors[1]],x$r.loading[,factors[2]], row.names( x$r.loading), cex=0.9, pos=4, col="red")  
plot(x$q.loading[,factors[1]],x$q.loading[,factors[2]],xlab=xlab,ylab=ylab,main="Q- Mode FA",col="black")
 text(x$q.loading[,factors[1]],x$q.loading[,factors[2]], row.names( x$q.loading), cex=0.9, pos=4, col="black")  
}

if(plot=="all"){
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(x$loadings[,factors[1]],x$loadings[,factors[2]],xlab=xlab,ylab=ylab,main="Loadings:R- and Q- Mode FA",col="black")
 text(x$r.loading[,factors[1]],x$r.loading[,factors[2]], row.names( x$r.loading), cex=0.9, pos=2, col="red")  
#points(x$q.loading[,factors[1]],x$q.loading[,factors[2]],pch=2, main="R- and Q- Mode FA")
 text(x$q.loading[,factors[1]],x$q.loading[,factors[2]], row.names( x$q.loading), cex=0.9, pos=2,col="blue")  
plot(x$r.loading[,factors[1]],x$r.loading[,factors[2]],xlab=xlab,ylab=ylab, main="R-Mode FA",col="red")
 text(x$r.loading[,factors[1]],x$r.loading[,factors[2]], row.names( x$r.loading), cex=0.9, pos=4, col="red")  
plot(x$q.loading[,factors[1]],x$q.loading[,factors[2]],xlab=xlab,ylab=ylab,main="Q- Mode FA",col="black")
 text(x$q.loading[,factors[1]],x$q.loading[,factors[2]], row.names( x$q.loading), cex=0.9, pos=4, col="blue")  
}

}
}
}
