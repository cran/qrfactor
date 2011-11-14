\name{qrfactor}
\alias{qrfactor}
\alias{qrfactor.default}
\alias{print.qrfactor}
\alias{summary.qrfactor}
\alias{plot.qrfactor}
\title{R package for Q and R mode Factor analysis }
\description{The qrfactor package simultaneously estimates Q and R mode factor analysis loadings and scores.}
\usage{
qrfactor(data,obs_start=NULL,mod_type='sd')
\method{qrfactor}{default}(data,obs_start=NULL,mod_type='sd')
\method{print}{qrfactor}(x,...)
\method{summary}{qrfactor}(object,...)
\method{plot}{qrfactor}(x,...)
}
\arguments{
\item{data}{A numeric design matrix for the model. All records must be numeric; it also accepts continuous data. Avoid using categorical variables and characters}
\item{obs_start}{ a numeric record where the row number starts. For example, in time series analysis a record can start from 1960. Insert 1 if your record  starts from observation 1. Set it to obs_start = NULL if you do not want specific row numbers assigned to loadings and scores }
\item{mod_type}{ Standardisation method that you want to use. Set it to "data" if you do not desire data transformation; set it to "sd" if you want the data to be standardised after centring it; set it to "n" if you want to divide the centred data by square root of the number of observations. The default is "sd" i.e. data.standard<-scale(data,center=TRUE,scale=TRUE)/sqrt(nrow(data))
   }
\item{x}{ an object of class \code{"qrfactor"}, i.e., a fitted model. }
\item{object}{ an object of class \code{"qrfactor"}, i.e., a fitted model. }
\item{...}{any other R parameters can be added}
}
\value{
Objects of the class that basically list its elements
\item{data}{ Original data for the model. All records must be numeric. It also accepts continous data}
\item{x.standard}{ it is the scale matrix of the original data}
\item{correlation}{ The correlation matrix for the data }
\item{eigen.value}{ eigen value of correlation matrix of the data}
\item{eigen.vector}{ eigen vector of correlation matrix of the data}
\item{diagonal.matrix}{ diagonal matrix of eigen vector }
\item{r.loading}{ R-mode loadings}
\item{q.loading}{ Q-mode loadings }
\item{loadings}{ combined loadings of R and Q on the same axis}
\item{q.scores}{ computed Q-mode scores }
\item{scores}{ combined R-mode and Q-mode scores on the same axis }
\item{rownames}{ row names of the loadings}
\item{variables}{ variables names of the loadings, of the original data}
}
\author{George Owusu}
\examples{
#data must be numeric, avoid using categorical and characters
data(UScereal, package="MASS")
variables=c("calories","protein","sodium","carbo","sugars","potassium")
data=UScereal[variables]
#create object with observation number starting with nothing. 
mod1 <- qrfactor(data)
#print object
mod1
#extract the loadings
loadings=mod1$loadings
#summary of the loadings
summary(mod1)
#plot the first two axes of the loadings
plot(mod1)


#2nd model with observations or records starting from zero and standardising
mod2<-qrfactor(data,obs_start=0,mod_type='sd')

#print object 2
mod2

#extract the loadings
loadings=mod2$loadings

#summary of the loadings
summary(mod2)
#plot the the first two axes of the loadings
plot(mod2)

#plotting other axes for mod2 for axes 3 and 4
plot(mod2$loadings[,3],mod2$loadings[,4])

}
\references{
Davis, J. C. (2002) Statistics and Data Analysis in Geology.  Third Edition. John Wiley & Sons
}
\keyword{Multivariate}