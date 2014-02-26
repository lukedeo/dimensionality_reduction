\name{local_linear_embedding}
\alias{local_linear_embedding}
\docType{package}
\title{
Locally Linear Embedding
}
\description{
Locally Linear Embedding (LLE) is a non-linear, spectral kernel method for embedding high dimensional data into lower dimensional spaces.
}
\synopsis{
}
\usage{
local_linear_embedding(X, k, d, verbose)
}
\arguments{
	\item{X}{ is a numeric data matrix (n x p)}
	\item{k}{ number of neighbors to consider for neighborhood construction, default = \code{6}}
	\item{d}{ is the desired dimension of the embedding, default = \code{2}}
	\item{verbose}{ provides a walkthrough of the status of the algorithm, default = \code{FALSE}}
}
\details{
The implementation of this algorithm follows the original publication from S.T. Roweis, L.K. Saul (2000), except using an SVD decomposition rather than carry out a sparse formation of a Gram matrix.
}
\value{
	\item{Y}{ The final embedding produced by the algorithm.}
	\item{description}{ A qualitative description of the parameters and method used.}
}
\references{
S.T. Roweis, L.K. Saul (2000) \emph{Nonlinear dimensionality reduction by locally linear embedding}, Science.
}

     

