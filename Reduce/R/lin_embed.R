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

mat_sqrt <- function(A){
    e <- eigen(A)
    V <- e$vectors
    B <- V %*% diag(sqrt(e$values)) %*% t(V)
    return(B)
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
         d = d) 
}

DIFFMAP <- function(X, d, t = 1.0, sigma = 1.0){
    n <- nrow(X)
    sumX <- rowSums(X ^ 2)
    
    K <- X %*% t(X)
    
    for(i in 1:n){
        K[i, ] <- (sumX - 2 * K[i, ]) / (2 * sigma ^ 2)
    }
    for(i in 1:n){
        K[i, ] <- exp(-(K[i, ] + sumX[i]))
    }
    
    
    
    
    
    
}

HLLE <- function(X, d, k){
    r <- ((d + 2) * (d + 1)) / 2
    if(k < r){
        stop("Error: formation of Hessian required k >= ((d + 2)(d + 1)) / 2.")
    }
    
    n <- nrow(X)
    dist <- pdist(X, X)
    diag(dist) <- Inf
    neighborhood <- matrix(0, n, n)
    for(i in 1:n){
        neighborhood[i, ] <- sort.int(dist[i, ], index.return=TRUE)$ix
    }
    
    dp <- d * (d + 1) / 2
    
    W <- matrix(0, n * dp, n)
    G <- matrix(0, n, n)
    for(i in 1:n){
        tmp_ind <- neighborhood[i, 1:k]
        thissrc <- X[tmp_ind, ]
        thissrc <- t(thissrc - colMeans(thisX))
        decomp <- svd(thissrc)
        
        U <- decomp$u
        Vpr <- decomp$v
        D <- decomp$d
    
        V <- Vpr[, 1:d]
        
        Yi <- matrix(1, k, r)
        Yi[, 2:(d+1)] <- V
        ct <- d + 1
        for(mm in 1:d)
        {
            nn <- c(1:(d + 1 - mm))
            Yi[, ct+nn] <- V[, nn] * V[, mm:d]
            ct <- ct + d - mm + 1
        }
        Yt <- mgs(Yi)
        Pii <- t(Yt[, (d+2):(ncol(Yt))])
        
        for(j in 1:dp)
        {
            if(sum(Pii[j, ]) > 0.0001)
            {
                tpp = Pii[j, ] / sum(Pii[j, ])
            }
            else
            {
                tpp = Pii[j, ]
            }
            G[tmp_ind, tmp_ind] <- G[tmp_ind, tmp_ind] + tpp %*% t(tpp)
        }
        
    }

    eigs <- eigen(G)    
    
    Y <- eigs$vectors
    to_select <- c((n - d):(n-1))
    Y <- (Y[, to_select]) * sqrt(n)
    
    plot(t(Y), pch=19, col=rainbow(N, start=0, end = .7))
    
    
    
    
    
    
    
    
    
    
    
    
    
    n <- nrow(X)
    k_orig <- k
    L <- as.matrix((nnwhich(X, k = c(1:k))))
    A <- adjacency(L)
    G <- matrix(0, n, n)
    mse <- rep(0, n)
    k <- rep(0, n)
    for(i in 1:n)
    {
        idx <- which(A[i, ])
        k[i] <- length(idx)
        thisX <- t(X[idx, ])
        thisX <- thisX - rowMeans(thisX)
        decomp <- svd(thisX)
        vals <- (decomp$d)
        D <- diag(vals)
        mse[i] <- sum(vals[(d+1):length(vals)])
        Vpr <- decomp$v
        V <- Vpr[, 1:d]
        Yi <- matrix(1, k[i], r)
        Yi[, 2:(d+1)] <- V
        ct <- d + 1
        for(mm in 1:d)
        {
            nn <- c(1:(d + 1 - mm))
            Yi[, ct+nn] <- V[, nn] * V[, mm:d]
            ct <- ct + d - mm + 1
        }
#         Yt <- qr(Yi)$qr
        Yt <- mgs(Yi)
        Pi <- Yt[, (d+2):(ncol(Yt))]
        G[idx, idx] <- G[idx, idx] + Pi %*% t(Pi)
    }
    eigs <- eigen(G, symmetric=TRUE)
    Y <- eigs$vectors
    to_select <- c((n - d):(n-1))
    Y <- t(Y[, to_select]) * sqrt(n)
    R <- t(Y) %*% Y
    R2 <- mat.sqrt(R)
    Y <- Y %*% R2
    plot(t(Y), pch=19, col=rainbow(N, start=0, end = .7))
    
    
    
}

LLE <- function(data, d, k) {
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
        w = solve_singular(G, ones)
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
#     return(eig_appx$vectors[, to_select])
    list(Y = null_space(M, d),
         X = x,
         k = k,
         d = d) 
}

manifold <- function(x, d, k, method = "standard", scale = TRUE, heat = 2.0,...) UseMethod("manifold")

manifold.default <- function(x, d, k, method = "normal", scale = TRUE, heat = 2.0,...){
    methods <- c("hessian", "standard", "normal", "laplacian")
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
        reduction <- LLE(x, d, k)
    }
    if(method == "laplacian")
    {
        reduction <- LEIGENMAP(x, d, k, heat)
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






