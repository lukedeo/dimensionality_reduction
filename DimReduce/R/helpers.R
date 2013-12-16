meta_criterion <- function(original_graph, reduced_graph)
{
    n <- nrow(original_graph)

    M_list <- rep(0, n)

    for (i in 1:n) {

       M_k <- length(intersect(which(reduced_graph[i, ] > 0), which(original_graph[i, ] > 0)))
       M_k <- M_k / sum(original_graph[i, ] > 0) - (1 / (n - 1)) * sum(reduced_graph[i, ] > 0)
       M_list[i] <- M_k
    }
    return((M_list))
}

knn_criterion <- function(D_original, D_reduced, k) {
    if(length(k) == 1)
    {
        original_graph <- neighbor_graph(D_original, k, 0)
        reduced_graph <- neighbor_graph(D_reduced, k, 0)
        return(mean(meta_criterion(original_graph, reduced_graph)))
    }
    else
    {
        ks <- length(k)
        mks <- rep(0, ks)
        for(i in 1:ks)
        {
            mks[i] <- knn_criterion(D_original, D_reduced, k[i])
        }
        return(mks)
    }

}

pdist <- function(A,B) {
    an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
    bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
    m = nrow(A)
    n = nrow(B)
    tmp = matrix(rep(an, n), nrow=m)
    tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
    sqrt( tmp - 2 * tcrossprod(A,B) )
}

#minimum norm solution to linear system
solve_singular <- function(x, y) {
    # solves:  x %*% b = y
    d = svd(x)
    # min-norm solution
    b.min = d$v %*% diag(1/d$d, length(d$d)) %*% t(d$u) %*% y[1:ncol(x)]
    return(b.min)
}

swiss_roll <- function(n = 400, noisy = FALSE) {
    r <- seq(0, 1, length.out = n)
    t <- (3*pi/2)*(1+2*r)
    x <- t*cos(t) + ifelse(noisy, rnorm(n, 0, .5), 0.0)
    y <- t*sin(t) + ifelse(noisy, rnorm(n, 0, .5), 0.0)
    z <- 20 * runif(n, 0, 1)
    return(scale(as.matrix(cbind(x, y, z))))
}

s_curve <- function(n = 400, noisy = FALSE)
{
    X <-  matrix(0, ncol = 3, nrow = n)
    s <- runif(n)
    t <- seq(0, 2, length.out=n)
    X[, 1] <- -cos(1.5 * pi * t)
    X[, 2] <- s
    X[t <= 1, 3] <- 2 * (-sin(1.5 * pi * t[t <= 1]))  + rnorm(n, 0, .3)[t <= 1]
    X[t > 1, 3] <- 2 * (2 + sin(1.5 * pi * t[t > 1])) + rnorm(n, 0, .3)[t > 1]
    return(scale(X))
}




compare_reductions <- function(X_orig, ...){
    candidates <- list(...)
    D_orig <- fastPdist(X_orig, X_orig)
    kvals <- c(1:10, seq(15, 30, 5), 40, 50)
    LC <- matrix(0, length(candidates), length(kvals))
    idx <- 1
    cols <- rainbow(length(candidates))
    for(manifold in candidates){
        Dnew <- manifold$Y
        Dnew <- fastPdist(Dnew, Dnew)
        LC[idx, ] <- knn_criterion(D_orig, Dnew, kvals)
        idx <- idx + 1
    }
    idx <- 1
    minval <- min(LC)
    minval <- ifelse(minval < 0, - 1.14 * abs(minval) , .84 * minval)
    
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




