library(spatstat)
solve.sing <- function(x, y) {
    # solves:  x %*% b = y
    d = svd(x)
    # min-norm solution
    b.min = d$v %*% diag(1/d$d, length(d$d)) %*% t(d$u) %*% y[1:ncol(x)]
    return(b.min)
}

null_space <- function(M, k, symmetric = TRUE) {
    n <- nrow(M)
    eig_appx = eigen(M, symmetric)
    to_select = c((n - d):(n-1))
    eig_appx$vectors[, to_select]
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

LLE_ <- function(data, d, k) {
    x = data
    cat("\nForming Euclidean Neighborhoods...")
    L = as.matrix((nnwhich(data, k = c(1:k))))
    cat("Done.")
    n = nrow(data)
    W = matrix(0, n, n)
    G = matrix(0, k, k)
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
        w = solve.sing(G, ones)
        w = w / sum(w)
        W[i, L[i, ]] = w
        cat(sprintf("\r%.2f%% complete.", (i / n*100)))
    }
    W = -W
    for(i in 1:n) {
        W[i, i] = W[i, i] + 1
    }
    cat("\nForming Objective Functional...")
    M = t(W)%*%(W)
    cat("Done.")
    cat("\nComputing eigen-decomposition...")
    eig_appx = eigen(M, symmetric=TRUE)
    cat("Done.\n")
    to_select = c((n - d):(n-1))
#     return(eig_appx$vectors[, to_select])
    list(Y = eig_appx$vectors[, to_select],
         X = x,
         k = k,
         d = d) 
}


hessian.LLE_ <- function(data, d, k) {
    x = data
    cat("\nForming Euclidean Neighborhoods...")
    L = as.matrix((nnwhich(data, k = c(1:k))))
    cat("Done.")
    n = nrow(data)
    W = matrix(0, n, n)
    G = matrix(0, k, k)
    M = matrix(0, k, n)
    
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
        w = solve.sing(G, ones)
        w = w / sum(w)
        W[i, L[i, ]] = w
        cat(sprintf("\r%.2f%% complete.", (i / n*100)))
    }
    W = -W
    for(i in 1:n) {
        W[i, i] = W[i, i] + 1
    }
    cat("\nForming Objective Functional...")
    M = t(W)%*%(W)
    cat("Done.")
    cat("\nComputing eigen-decomposition...")
    eig_appx = eigen(M, symmetric=TRUE)
    cat("Done.\n")
    to_select = c((n - d):(n-1))
    #     return(eig_appx$vectors[, to_select])
    list(Y = eig_appx$vectors[, to_select],
         X = x,
         k = k,
         d = d) 
}






lin_embed <- function(x, d, k, scale = TRUE...) UseMethod("lin_embed")

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




















