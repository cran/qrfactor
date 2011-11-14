rq <- function(data,obs_start=NULL,mod_type='sd')
{
x<-data
#standardize the data default: by standard deviation 
x_standard<-scale(x,center=TRUE,scale=TRUE)/sqrt(nrow(x))
#no standardisation; use original data
 if (mod_type == 'data') { 
x_standard<-as.matrix(x)
}

#standardize by standard deviation
if (mod_type == 'sd') { 
x_standard<-scale(x,center=TRUE,scale=TRUE)/sqrt(nrow(x))

}

#standardize by square root of n
if (mod_type =='n') { 
x_standard<-scale(x,center=TRUE,scale=FALSE)/sqrt(nrow(x))
}

#develop correlation matrix
corr_matrix<-cor(x_standard)

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
qscores<-corr_matrix%*%eigen_vector
combined.scores<-rbind(rscores,qscores)

#develop row and column names for the loadings
colnames=colnames(x)
nrow=nrow(x)
total_obs=length(colnames)+nrow
listnames=""
i=1
#variables names
variables=""
while (i <= total_obs) {
   if (i<=(length(colnames))) {
listnames=paste(listnames,colnames[i],sep=",")
variables=paste(variables,colnames[i],sep=",")
}
else {
number=i-length(colnames)
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
rownames(all_loadings) <- c(rownames)

variables=strsplit(variables,',')
variables=variables[[1]]
variables=variables[2:length(variables)]
rownames(R_loading) <- c(variables)

rownames(Q_mode_loading) <- rownames(Q_mode_loading, do.NULL = FALSE, prefix = "")

#check for validity of the model
#scale.matrix<-Q_mode_loading%*%diag((sqrt(eigen_val)^(-0.5)))%*%t(R_loading)


list(correlation=corr_matrix,eigen.vector=eigen_vector,eigen.value=eigen_val,
diagonal.matrix=diagonal_matrix,r.loading=R_loading,q.loading=Q_mode_loading,
combined.loadings=all_loadings,r.scores=rscores,q.scores=qscores,combined.scores=combined.scores,data=x,
rownames=rownames,variables=variables,x.standard=x_standard,loadings=all_loadings,scores=combined.scores)
}

#generic function
qrfactor<-function(data,obs_start=NULL,mod_type="sd") UseMethod ("qrfactor")

#default function
qrfactor.default<-function(data,obs_start=NULL,mod_type="sd")
{
#convert file data to matrix
x<-as.matrix(data)

factor<-rq(data,obs_start,mod_type)

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

plot.qrfactor<-function(x,...)
{
plot(x$combined.loadings[,1],x$combined.loadings[,2],xlab="Factor 1",ylab="Factor 2")
 text(x$combined.loadings[,1],x$combined.loadings[,2], row.names( x$combined.loadings), cex=0.6, pos=4, col="red")  
}

