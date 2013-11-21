pdist <- function(A,B) {
    an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
    bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
    m = nrow(A)
    n = nrow(B)
    tmp = matrix(rep(an, n), nrow=m)
    tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
    sqrt( tmp - 2 * tcrossprod(A,B) )
}
solve_singular <- function(x, y) {
    # solves:  x %*% b = y
    d = svd(x)
    # min-norm solution
    b.min = d$v %*% diag(1/d$d, length(d$d)) %*% t(d$u) %*% y[1:ncol(x)]
    return(b.min)
}

compare_reductions <- function(...){
    candidates <- list(...)
    kvals <- c(1:10, seq(15, 30, 5), 40, 50)
    LC <- matrix(0, length(candidates), length(kvals))
    idx <- 1
    cols <- rainbow(length(candidates))
    for(manifold in candidates){
        LC[idx, ] <- LC_metacriterion(manifold, kvals)
        idx <- idx + 1
    }
    idx <- 1
    minval <- 1.14 * min(LC)
    maxval <- 1.14 * max(LC)
    
    plot(kvals, LC[1, ], type = "b", xlim = c(1, 50), ylim = c(minval, maxval), col = cols[1], ylab = "Adjusted Criterion", xlab = "Neighborhood Size", log = "x", main = "Methods Comparison")
    for(manifold in candidates){
        if(idx != 1){
            points(kvals, LC[idx, ], type = "b", col = cols[idx])
        }
        idx <- idx + 1    
    } 
    names <- NULL
    coresp_colors <- NULL
    idx <- 1
    lty <- NULL
    for(manifold in candidates){
        names <- cbind(names, manifold$description)
        coresp_colors <- cbind(coresp_colors, cols[idx])
        lty<-cbind(lty, 1)
        idx <- idx + 1
    }
    legend("topright", cex = 0.8, legend=names, col=coresp_colors, lty=lty)
}

random_projection <- function(X, d, type = 1){
    X_new <- t(X)
    dim <- nrow(X_new)
    R <- randmat(dim, d, type)
    Y <- t(1 / sqrt(d) * t(R) %*% X_new)
    if(type==1){proj="Gaussian"}
    if(type==2){proj="+/- 1 Uniform"}
    list(Y = Y,
         X = X_new,
         d = d,
         description = paste("Random Projection,", proj))
}

diffmap <- function(X, d, t = 1.0, sigma = -1.0){
    n <- nrow(X)
    dist <- pdist(X, X)
    if(sigma == -1)
    {
        sigma <- (max(dist))
    }
    K <- exp(-(dist * dist) / (2 * sigma ^ 2))

    for(i in 1:n)
    {
        K[i, ] <- K[i, ] / sum(K[i, ])
    }
    
    K <- K %^% t
    
    eigs <- eigen(K)
    
    Y <- Re(eigs$vectors[, 2:(d+1)])
    
    list(Y = Y,
         X = X,
         sigma = sigma,
         t = t,
         description = paste("Diffusion Map, t = ", t, ", sigma = ", sigma, sep = ""),
         d = d) 
}




