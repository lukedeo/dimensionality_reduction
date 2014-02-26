\name{diffusion_map}
\alias{diffusion_map}
\docType{package}
\title{
Diffusion Maps
}
\description{
Locally Linear Embedding (LLE) is a non-linear, spectral kernel method for embedding high dimensional data into lower dimensional spaces.
}
\synopsis{
}
\usage{
diffusion_map(X, d, t, sigma, verbose)
}
\arguments{
	\item{X}{ is a numeric data matrix (n x p)}
	\item{d}{ is the desired dimension of the embedding, default = \code{2}}
	\item{t}{ number of timesteps to run the diffusion process (no functionality yet), default = \code{2.0}}
	\item{t}{ number of timesteps to run the diffusion process (no functionality yet), default = \code{2.0}}

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

     

