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
dot <- function(x, y){
    return(sum(x * y))
}
mat.sqrt<-function(A){
    ei<-eigen(A)
    d<-ei$values
    d<-(d+abs(d))/2
    d2<-sqrt(d)
    ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
    return(ans)
}
mgs <- function(A){
    n <- ncol(A)
    A <- as.matrix(A)
    R <- matrix(0, nrow = n, ncol = n)
    for(i in 1:(n-1))
    {
        R[i,i] = norm(as.matrix(A[ ,i]))
        A[,i] = A[,i]/R[i,i] 
        for(j in (i+1):n)
        {
            R[i,j] = t(A[,i]) %*% A[,j]
            A[,j] = A[,j] - R[i,j] * A[,i]
        }
    }
    R[n,n] = norm(as.matrix(A[,n]))
    A[,n] = A[,n]/R[n,n]
    Q = A
    return(Q)
}
diag_inverse <- function(A){
    for(i in 1:nrow(A))
    {
        A[i, i] <- 1 / A[i, i]
    }
    return(A)
}
null_space <- function(M, d, symmetric = TRUE) {
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
adjacency <- function(neighbor_list){
    L <- neighbor_list
    n <- nrow(L)
    A <- matrix(0, n, n)
    for(i in 1:n){
        A[i, L[i, ]] <- 1
    }
    A <- (A + t(A)) > 0
    return(A)
}
randmat <- function(m, n, type = 1){
    if(type == 1) {
        A <- matrix(rnorm(m * n), m, n)
        return(A)
    }
    if(type == 2) {
        A <- matrix(sample(c(1, -1), m * n, replace = TRUE), m, n)        
        return(A)
    }
    else
    {
        cat("Error: Matrix must be of type 1 or 2.")
    }
}
LC_metacriterion <- function(embedding, k){
    if(length(k) == 1){
        X <- embedding$X
        Y <- embedding$Y
        n <- nrow(X)
        N_orig <- as.matrix((nnwhich(X, k = c(1:k))))
        N_embed <- as.matrix((nnwhich(Y, k = c(1:k))))
        N_union <- N_orig == N_embed
        N_union <- rowSums(N_union)
        N_k <- mean(N_union)
        M_k <- N_k / k
        M_k_adjusted <- M_k - (k / (n - 1))
        return(M_k_adjusted)
    }
    else{
        adjusted_vals <- k
        idx <- 1
        for(value in k){
            adjusted_vals[idx] <- LC_metacriterion(embedding, value)
            idx <- idx + 1
        }
        return(adjusted_vals)
    }
}



compare_reduction <- function(...){
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






mat_sqrt <- function(A){
    e <- eigen(A)
    V <- e$vectors
    B <- V %*% diag(sqrt(e$values)) %*% t(V)
    return(B)
}

RANDOM_PROJECTION <- function(X, d, type = 1){
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

LEIGENMAP <- function(X, d, k, heat = 2.0){
    n <- nrow(X)
    dist <- pdist(X, X)
    dist_copy <- dist
    graph <- matrix(0, n, n)
    for(i in 1:n){
        dist[i, i] <- Inf
        for(j in 1:k){
            idx <- which.min(dist[i, ])
            graph[i, idx] <- 1.0
            graph[idx, i] <- graph[i, idx]
            dist[i, idx] <- Inf
        }
    }
    if(heat != 0){
        non_zero <- graph>0
        graph[non_zero] = graph[non_zero] * exp((-dist_copy[non_zero]) / heat)
    }
    weight <- diag(rowSums(graph))
    laplacian <- weight - graph
    
    eigs <- eigen(solve(weight) %*% laplacian)
    to_select <- c((n - d):(n-1))
    Y <- eigs$vectors[, to_select] 
    list(Y = Y,
         X = X,
         k = k,
         description = paste("Laplacian Eigenmap, k = ", k, ", heat = ", heat, sep = ""),
         d = d) 
}

DIFFMAP <- function(X, d, t = 1.0, sigma = -1.0){
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


LLE <- function(X, d, k, symmetric = FALSE, positive = FALSE) {
    n <- nrow(X)
    dist <- pdist(X, X)
    dist_copy <- dist
    graph <- matrix(0, n, n)
    W <- matrix(0, n, n)
    for(i in 1:n) {
        dist[i, i] <- Inf
        for(j in 1:k){
            idx <- which.min(dist[i, ])
            graph[i, idx] <- 1.0
            if(symmetric){
                graph[idx, i] <- graph[i, idx]   
            }
            dist[i, idx] <- Inf
        }
    }
    graph <- graph>0
    for(i in 1:n) {
        one_vec <- matrix(1, sum(graph[i, ]), 1)
        Z <- t(X[graph[i, ], ]) - c(X[i, ])
        C <- t(Z) %*% Z
        diag(C) <- diag(C) + 0.001 * tr(C)
        w <- solve_singular(C, one_vec)
        w <- w / sum(w)
        if(sum(is.nan(w))>0 | sum(is.infinite(w))>0){
            print(i)
            break
        }
        if(positive)
        {
            w[w<0] <- 0
            W[i, graph[i, ]] <- -w
        }
        else{
            W[i, graph[i, ]] <- -w
        }
    }
    diag(W) <- 1
    decomp <- svd(W)
    if(!symmetric){sym="non-"}
    else{sym=""}
    selection = c((n - d):(n-1))
    Y <- decomp$v[, selection]
    list(Y = Y,
         X = X,
         k = k,
         description = paste("LLE, k = ", k, ", ", sym, "symmetrized adjacency.", sep = ""),
         d = d) 
}

manifold <- function(x, d, k, method = "standard", scale = TRUE, heat = 1, sigma = -1, t = 3, symmetric = FALSE,...) UseMethod("manifold")

manifold.default <- function(x, d, k = 3, method = "normal", scale = TRUE, heat = 1, sigma = -1, t = 3, symmetric = FALSE,...){
    methods <- c("hessian", "standard", "normal", "laplacian", "projection")
    if(scale)
    {
        x <- scale(as.matrix(x))
    }
    else
    {
        x <- as.matrix(x)
    }
    if((method == "normal") | (method == "standard"))
    {
        reduction <- LLE(x, d, k, symmetric=symmetric)
    }
    if(method == "laplacian")
    {
        reduction <- LEIGENMAP(x, d, k, heat)
    }
    if(method == "diffusion")
    {
        reduction  <- DIFFMAP(x, d, t, sigma)
    }
    if(method == "projection")
    {
        reduction <- RANDOM_PROJECTION(X=X, d=d)
    }
    reduction$call <- match.call()
    class(reduction) <- "manifold"
    reduction
}

manifold.formula <- function(formula, data=list(), d, k, scale = TRUE,...){
    mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    est <- manifold.default(x, d, k, scale = TRUE,...)
    est$call <- match.call()
    est$formula <- formula
    est
}






