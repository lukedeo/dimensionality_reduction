reduce <- function(x, d,...) UseMethod("reduce")

reduce.default <- function(x, d, method = "lle", scale = TRUE, ...){
    if(scale)
    {
        x <- scale(as.matrix(x))
    }
    else
    {
        x <- as.matrix(x)
    }
    if(method == "lle")
    {
        reduction <- local_linear_embedding(x, d, ...)
    }
    if(method == "laplacian")
    {
        reduction <- laplacian_eigenmap(x, d, ...)
    }
    if(method == "diffusion")
    {
        reduction  <- diffusion_map(x, d, ...)
    }
    if(method == "projection")
    {
        reduction <- random_projection(x, d, ...)
    }
    # if(method == "mds")
    # {
    #     reduction <- cmds()
    # }
    reduction$call <- match.call()
    class(reduction) <- "manifold"
    reduction
}
