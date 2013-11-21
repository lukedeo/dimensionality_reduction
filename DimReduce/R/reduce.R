embedding <- function(X = NULL, A = NULL, D = NULL, d = 2,...) UseMethod("embedding")

embedding.default <- function(X = NULL, A = NULL, D = NULL, d = 2, method = "lle", scale = TRUE, ...){
    if(scale)
    {
        X <- scale(as.matrix(x))
    }
    else
    {
        x <- as.matrix(x)
    }
    if(method == "lle")
    {
        reduction <- local_linear_embedding(X = X, d = d, ...)
    }
    if(method == "laplacian")
    {
        reduction <- laplacian_eigenmap(X = X, d = d, ...)
    }
    if(method == "diffusion")
    {
        reduction  <- diffusion_map(X = X, d = d, ...)
    }
    if(method == "projection")
    {
        reduction <- random_projection(X = X, d = d, ...)
    }
    if(method == "mds")
    {
        reduction <- cmds(D = D, d = d)
    }
    if(method == "boxcox")
    {   
        reduction <- boxcox(D = D, A = A, d = d, ...)
    }
    if(method == "isomap")
    {
        reduction <- isomap(X = X, d = d, ...)
    }
    reduction$call <- match.call()
    class(reduction) <- "embedding"
    reduction
}

plot.embedding <- function(x,...)
{
    if(ncol(x$Y) == 2)
    {
        plot(x$Y, main = x$description, xlab = "Feature 1", ylab = "Feature 2", ...);
    }
    else
    {
        cat("Printing not defined for an embedding of dimension > 2.")
    }
    
}
