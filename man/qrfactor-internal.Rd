\name{rq}
\alias{rq}
\title{Internal Q and R mode Factor analysis function}
\description{An internal function to simultaneoulsy estimates of Q and R mode factor analysis loadings and scores. It is for internal use only}
\usage{
rq(data,obs_start,mod_type)
}
\arguments{
\item{data}{a numeric design matrix for the model. }
\item{obs_start}{ a numeric record where the row names start.  }
\item{mod_type}{ Standardisation method that you want to use. Set it to "data" if you do not desire data transformation; set it to "sd" if you want the data to be standardised after centring it; set it to "n" if you want to divide the centred data by square root of the number of observations. The default is "sd" i.e. data.standard<-scale(data,center=TRUE,scale=TRUE)/sqrt(nrow(data))
   }
}
\value{
An object of class  basically a list including  elements
\item{data}{ Original data for the model. All records must be numeric. It also accepts continous data}
\item{x.standard}{ it is the scale matrix}
\item{correlation}{ The correlation matrix for the data }
\item{eigen.value}{ eigen value of correlation matrix of the data}
\item{eigen.vector}{ eigen vector of correlation matrix of the data}
\item{diagonal.matrix}{ diagonal matrix of eigen vector }
\item{r.loading}{ R-mode loadings}
\item{q.loading}{ Q-mode loadings }
\item{loadings}{ combined loadings of R and Q on the same axis}
\item{combined.loadings}{ combined loadings of R and Q on the same axis}
\item{q.scores}{ computed Q-mode scores }
\item{scores}{ combined R-mode and Q-mode scores on the same axis }
\item{rownames}{ row names of the loadings}
\item{variables}{ variables names of the loadings, of the original data}
}
\author{George Owusu}
\keyword{Multivariate}