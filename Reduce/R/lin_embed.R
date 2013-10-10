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

LEIGENMAP <- function(X, d, k){
    n <- nrow(X)
    L <- as.matrix((nnwhich(X, k = c(1:k))))
    A <- adjacency(L)
    D <- diag(rowSums(A))
    L <- D - A
    L <- diag_inverse(D) %*% L
    eigs <- eigen(L)
    to_select <- c((n - d):(n-1))
    Y <- eigs$vectors[, to_select] 
    list(Y = Y,
         X = X,
         k = k,
         d = d) 
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

manifold <- function(x, d, k, method = "standard", scale = TRUE,...) UseMethod("manifold")

manifold.default <- function(x, d, k, method = "normal", scale = TRUE,...){
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
        reduction <- LEIGENMAP(x, d, k)
    }
    reduction$call <- match.call()
    class(reduction) <- "manifold"
    reduction
}

manifold.formula <- function(formula, data=list(), d, k, scale = TRUE,...)
{
    mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    est <- manifold.default(x, d, k, scale = TRUE,...)
    est$call <- match.call()
    est$formula <- formula
    est
}






