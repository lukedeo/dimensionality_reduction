\name{boxcox}
\alias{boxcox}
\docType{package}
\title{
Generalized Multidimensional Scaling.
}
\description{
Multidimensional Scaling using the Boxcox family of stress functions
}
\synopsis{
}
\usage{
boxcox(D, A, d = 2, X1 = NULL, \dots)	
}
\arguments{
	\item{D}{ is a (not necessarily symmetric) distance matrix.}
	\item{A}{ is an adjacency matrix associated with \code{D} -- can be k-nearest neighbors or \eqn{\varepsilon}{epsilon} neighborhoods.}
	\item{X1}{is an initial configutation of the embedding, default = \code{NULL}}
	\item{d}{ is the desired dimension of the embedding, default = \code{2}}
	\item{\dots}{Any of the following additional arguments:}
	\item{}{\describe{
		\item{\code{cmds_start}}{Whether or not we begin by performing Classical MDS -- if we do, it becomes our initial guess at an optimal embedding. default = \code{TRUE} (\emph{recommended})}
		\item{\code{lambda}}{default = \code{TRUE}}
		\item{\code{mu}}{default = \code{TRUE}}
		\item{\code{nu}}{default = \code{FALSE}}
		\item{\code{tau}}{default = \code{TRUE}}
		\item{\code{niter}}{How many iterations for the optimization procedure, default = \code{500}}
		\item{\code{sample_rate}}{Specifies how often we update the adjusted metacriterion for preservation. Also decides how often we sample for the (metrically) best embedding, default = \code{100}}
		\item{\code{bfgs}}{If set to \code{TRUE}, use \emph{full} BFGS (Quasi-Newton) optimization. Recommended only for small iteration counts or getting out of local minima from previous embeddings, default = \code{FALSE}}
		\item{\code{adaptive}}{If \code{FALSE} and \code{bfgs = FALSE}, then uses Armijo Backtracking Line-Search to determine a stepsize, default = \code{TRUE}(\emph{highly recommended})}
		\item{\code{scale_out}}{Should outputs be scaled to zero mean and unit variance, default = \code{TRUE}}
	}}
}
\details{
The implementation of this algorithm and the Local Continuity Metacriterion follow the specification of the Boxcox Stress functions a la Chen and Buja (2012). 
}
\value{
	\item{Y}{ The final embedding produced by the algorithm.}
	\item{best}{ The best embedding as produced through sampling during optimization, judged by the Local Continuity Metacriterion.}
	\item{description}{ A qualitative description of the parameters used.}
}
\references{
Chen, L. and Buja, A. (2012) \emph{Stress Functions for Nonlinear Dimension Reduction, Proximity Analysis, and Graph Drawing}, Journal of Machine Learning Research.
}

     

