\name{qrfactor}
\alias{qrfactor}
\alias{qrfactor.default}
\alias{print.qrfactor}
\alias{summary.qrfactor}
\alias{plot.qrfactor}
\title{R package for Q and R mode Factor analysis }
\description{The qrfactor package simultaneously estimates Q and R mode factor analysis loadings and scores.}
\usage{
qrfactor(data,obs_start=NULL)
\method{qrfactor}{default}(data,obs_start)
\method{print}{qrfactor}(x,...)
\method{summary}{qrfactor}(object,...)
\method{plot}{qrfactor}(x,...)
}
\arguments{
\item{data}{ a numeric design matrix for the model. All records be numeric. It also accepts continous data. Avoid using categorical and characters }
\item{obs_start}{ a numeric record where the row names start. For example, in time series analysis a record can start from 1960. Insert 1 if your record  tarts from observation 1. Set it to obs_start = NULL if you do not want specific row names assigned to loadings and scores }
\item{x}{ an object of class \code{"qrfactor"}, i.e., a fitted model. }
\item{object}{ an object of class \code{"qrfactor"}, i.e., a fitted model. }
\item{...}{any other R parameters can be added}
}
\value{
An object of class  basically a list including  elements
\item{data}{ Original data for the model. All records be numeric. It also accepts continous data}
\item{x.standard}{ it is the scale matrix: scale matrix = scale(data,center=TRUE,scale=TRUE)/sqrt(nrow(data))}
\item{correlation}{ correlation matrix for the data }
\item{eigen.vector}{ eigen vector of correlation matrix of the data}
\item{diagonal.matrix}{ diagonal matrix of eigen vector }
\item{r.loading}{ R-mode loadings}
\item{q.loading}{ Q-mode loadings }
\item{loadings}{ combined loadings of R and Q on the same axis}
\item{combined.loadings}{ combined loadings of R and Q on the same axis}
\item{r.scores}{ computed R-mode scores }
\item{q.scores}{ computed Q-mode scores }
\item{scores}{ combined R-mode and Q-mode scores on the same axis }
\item{combined.scores}{ combined R-mode and Q-mode scores on the same axis }
\item{rownames}{ row names of the loadings}
\item{variables}{ variables names of the loadings, of the original data}
}
\author{George Owusu}
\examples{
#data must be numeric, avoid using categorical and characters
data(UScereal, package="MASS")
#create object with observation number starting with nothing. 
variables=c("calories","protein","sodium","carbo","sugars","potassium")
data=UScereal[variables]
mod1 <- qrfactor(data,obs_start=NULL)
#print object
mod1
#extract the loadings
loadings=mod1$loadings
#summary of the loadings
summary(mod1)
#plot the the first two axes of the loadings
plot(mod1)


#2nd model with observations or records starting from zero
mod2<-qrfactor(data,obs_start=0)

#print object 2
mod2

#extract the loadings
loadings=mod2$loadings

#summary of the loadings
summary(mod2)
#plot the the first two axes of the loadings
plot(mod2)

#ploting other axes for mod2 for axes 3 and 4
plot(mod2$loadings[,3],mod2$loadings[,4])

}
\references{
Davis, J. C. (2002) Statistics and Data Analysis in Geology.  Third Edition. John Wiley & Sons
}
\keyword{Multivariate}