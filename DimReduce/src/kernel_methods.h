#include <RcppArmadillo.h>
#include "fastPdist.h"
#include "neighbor_graph.h"





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


    // graph <- graph>0
    // for(i in 1:n) {
    //     one_vec <- matrix(1, sum(graph[i, ]), 1)
    //     Z <- t(X[graph[i, ], ]) - c(X[i, ])
    //     C <- t(Z) %*% Z
    //     diag(C) <- diag(C) + 0.001 * tr(C)
    //     w <- solve_singular(C, one_vec)
    //     w <- w / sum(w)
    //     if(sum(is.nan(w))>0 | sum(is.infinite(w))>0){
    //         print(i)
    //         break
    //     }
    //     if(positive)
    //     {
    //         w[w<0] <- 0
    //         W[i, graph[i, ]] <- -w
    //     }
    //     else{
    //         W[i, graph[i, ]] <- -w
    //     }
    // }
    // diag(W) <- 1
    // decomp <- svd(W)
    // if(!symmetric){sym="non-"}
    // else{sym=""}
    // selection = c((n - d):(n-1))
    // Y <- decomp$v[, selection]
    // list(Y = Y,
    //      X = X,
    //      k = k,
    //      description = paste("LLE, k = ", k, ", ", sym, "symmetrized adjacency.", sep = ""),
    //      d = d) 
//Standard Implementation of Linear Embedding
SEXP local_linear_embedding(arma::mat X, int k, int d = 2)
{
	unsigned int n = X.n_rows;
	unsigned int p = X.n_cols;
	arma::mat Adj = fastPdist(X, X);
	Adj = neighbor_graph(Adj, k);
	arma::mat W(n, n), U(n, n), V(n, n);

	arma::vec ones_vec(n), s();
	ones_vec.ones();

	for (int i = 0; i < n; ++i)
	{
		arma::mat Z = X.rows(Adj.row(i)).t() - X.row(i);
		amra::mat C = Z.t() * Z.
		C.diag() += 0.001 * arma::trace(C);
		arma::vec w = arma::solve(C, ones_vec);
		W.row(i).cols(Adj.row()) = -w;
	}
	W.diag().ones();

	svd(U, s, V, X, method = "dc");

	arma::mat Y = V.cols((n - d), (n - 1));

	fast_scale(Y);
    return Rcpp::List::create(
        Rcpp::Named("embedding") = (Y)
    );






}