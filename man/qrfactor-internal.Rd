\name{rq}
\alias{rq}
\title{Internal Q and R mode Factor analysis function}
\description{internal function to simultaneoulsy estimates of Q and R mode factor analysis loadings and scores. It is for internal use only}
\usage{
rq(data,obs_start)
}
\arguments{
\item{data}{a numeric design matrix for the model. }
\item{obs_start}{ a numeric record where the row names start.  }
}
\value{
An object of class  basically a list including  elements
\item{data}{ original data for the model. It must be numeric}
\item{x.standard}{ it is the scale matrix: scale matrix = scale(x,center=TRUE,scale=TRUE)/sqrt(nrow(x))}
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
\keyword{Multivariate}