solve_singular <- function(x, y) {
    # solves:  x %*% b = y
    d = svd(x)
    # min-norm solution
    b.min = d$v %*% diag(1/d$d, length(d$d)) %*% t(d$u) %*% y[1:ncol(x)]
    return(b.min)
}


M_times_v <- function(v, extra = NULL){
    W <- extra
    return((v - W %*% v) - t(W) %*% (v - W %*% v))
}


# computes the d major eigenvectors from a local reconstruction matrix W
# avoids explicit computation of W^T %*% W using ARPACK
null_space <- function(W, d) { 
    options = list(n=nrow(W), nev=d, ncv=n, sym=TRUE, which="LM", maxiter=200)
    eig_vectors <- arpack(M_times_v, extra = W, sym = TRUE, options = options)
    return(eig_vectors$vectors)
#     
#     n <- nrow(M)
#     eig_appx = eigen(M, symmetric)
#     to_select = c((n - d):(n-1))
#     eig_appx$vectors[, to_select]
}

tr <- function(mat) {
    if(ncol(mat) == nrow(mat))
    {
        trace = 0
        for (i in 1:nrow(mat)) {
            trace = trace + mat[i, i]
        }
        return (trace)
    }
    else{
        cat("\nInvalid Dims\n")
    }
}

LLE <- function(data, d, k) {
    x = data
    cat("\nForming Euclidean Neighborhoods...")
    L = as.matrix((nnwhich(data, k = c(1:k))))
    cat("Done.")
    n = nrow(data)
    
#     W = Matrix(0, nrow = n, ncol = n, sparse = TRUE) # this is very, very sparse.
    W = matrix(0, nrow = n, ncol = n) # this is very, very sparse.
    
    G = matrix(0, k, k) # this is not.
    ones = matrix(1, k, 1)
    cat("\nBuilding Local Reconstructions...\n")
    for(i in 1:n) 
    {
        for(j in 1:k) 
        {
            ind_1 = L[i, j]
            for(l in 1:k) 
            {
                ind_2 = L[i, l]
                G[j, l] = t((x[i, ] - x[ind_1, ])) %*% (x[i, ] - x[ind_2, ])
            }
        }    
        w = solve_singular(G, ones)
        w = w / sum(w)
        W[i, L[i, ]] = -w
        cat(sprintf("\r%.2f%% complete.", (i / n*100)))
    }
#     W = -W
    for(i in 1:n) {
        W[i, i] = W[i, i] + 1
    }
#     return(eig_appx$vectors[, to_select])
    list(Y = null_space(W, d),
         X = x,
         k = k,
         d = d) 
}

lin_embed <- function(x, d, k, scale = TRUE,...) UseMethod("lin_embed")

lin_embed.default <- function(x, d, k, scale = TRUE,...){
    if(scale)
    {
        x <- scale(as.matrix(x))
    }
    else
    {
        x <- as.matrix(x)
    }
    reduction <- LLE(x, d, k)
    reduction$call <- match.call()
    class(reduction) <- "lin_embed"
    reduction
}

lin_embed.formula <- function(formula, data=list(), d, k, scale = TRUE,...)
{
    mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    
    est <- lin_embed.default(x, d, k, scale = TRUE,...)
    est$call <- match.call()
    est$formula <- formula
    est
}






