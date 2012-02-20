\name{qrfactor}
\alias{qrfactor}
\alias{qrfactor.default}
\alias{print.qrfactor}
\alias{summary.qrfactor}
\alias{plot.qrfactor}
\title{R package for Q and R mode Factor analysis }
\description{The qrfactor package simultaneously estimates Q and R mode factor analysis loadings and scores.}
\usage{
qrfactor(data,scale='sd')
\method{qrfactor}{default}(data,scale='sd')
\method{print}{qrfactor}(x,...)
\method{summary}{qrfactor}(object,...)
\method{plot}{qrfactor}(x,factors=c(1,2),type="loading",plot="",...)
}
\arguments{
\item{data}{A numeric design matrix for the model. All records must be numeric; it also accepts continuous data. Avoid using categorical variables and characters}
\item{scale}{ Standardisation method that you want to use. Set it to "data" if you do not desire data transformation; set it to "sd" if you want the data to be standardised after centring it; set it to "n" if you want to divide the centred data by square root of the number of observations. Normal scaling of subtracting the centring the data and dividing it by standard dviation is "msd". Set it to "msd" if you want to estiamte PCA. The default is "sd" 
   }
\item{x}{ an object of class \code{"qrfactor"}, i.e., a fitted model. }
\item{object}{ an object of class \code{"qrfactor"}, i.e., a fitted model. }
\item{type}{The type of results one wants to plot. It takes "scores", "loadings", pca or eigenvectors. The deault is loadings. }
\item{plot}{The type of plots one desires. It takes "all" for all the 3 plots or "q" for q plot  or "r" for r plot or 'qr' for both q and r plots}
\item{factors}{list of factors one wants to plot. The default is factors=c(1,2). Please do not forget "c" in the list. }
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
\item{pca}{ pca loadings}
\item{pcascores}{ PCA scores }
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
mod2<-qrfactor(data,scale='sd')

#print object 2
mod2

#extract the loadings
loadings=mod2$loadings

#summary of the loadings
summary(mod2)
#plot the the first two axes of the loadings
plot(mod2)

#plotting other axes for mod2 for axes 3 and 4
plot(mod2,factors=c(3,4))
plot(mod2,factors=c(1,3))
#plotting scores fro axes 2 and 3
plot(mod2,factors=c(3,4),type='scores')

}
\references{
Davis, J. C. (2002) Statistics and Data Analysis in Geology.  Third Edition. John Wiley & Sons
}
\keyword{Multivariate}