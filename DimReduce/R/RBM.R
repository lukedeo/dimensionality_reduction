#### Code for normal RBM


sigmoid <- function(x) 
{
    return(1 / (1 + exp(-x)))
}
identity <- function(x) 
{
    return(x)
}

bernoulli <- function(x) 
{
    return(runif(length(x)) < x)
}

RBM <- setRefClass("RBM", 
    fields = list(
        weights = "matrix", 
        hidden_bias = "matrix", 
        visible_bias = "matrix", 
        visible_form = "function",
        num_hidden = "numeric",
        num_visible = "numeric"
    ),
    methods = list(
        initialize = function(visible, hidden, binary = TRUE, scale = 0.01)
        {
            num_hidden <<- hidden
            num_visible <<- visible

            weights <<- matrix(rnorm(visible * hidden, 0, scale), hidden, visible)
            hidden_bias <<- matrix(rnorm(hidden, 0, scale), 1, hidden)
            visible_bias <<- matrix(rnorm(visible, 0, scale), 1, visible)
            if(binary)
            {
                visible_form <<- sigmoid
            }
            else
            {
                visible_form <<- identity
            }
        },
        reset = function(visible = NULL, hidden = NULL, binary = TRUE, scale = 0.01)
        {
            if(is.null(visible)) visible <- num_visible
            if(is.null(hidden)) hidden <- num_hidden

            if(!is.numeric(hidden) | !is.numeric(visible)) error("Dimensions have to be integers.")

            num_hidden <<- hidden
            num_visible <<- visible

            weights <<- matrix(rnorm(visible * hidden, 0, scale), hidden, visible)
            hidden_bias <<- matrix(rnorm(hidden, 0, scale), 1, hidden)
            visible_bias <<- matrix(rnorm(visible, 0, scale), 1, visible)
            if(binary)
            {
                visible_form <<- sigmoid
            }
            else
            {
                visible_form <<- identity
            }
        },
        expected_hidden = function(visible, bias = 0)
        {
            # matrix(m$hidden_bias, 10, 3)
            if(!is.null(nrow(visible)))
            {
                return(t(sigmoid(weights %*% t(visible))) + t(matrix(hidden_bias, num_hidden, nrow(visible))) + bias)
            }
            else
            {
                return(t(sigmoid(weights %*% t(visible))) + hidden_bias + bias)    
            }
            
        },
        expected_visible = function(hidden, bias = 0)
        {
            if(!is.null(nrow(hidden)))
            {
                return(visible_form((hidden %*% weights) + t(matrix(visible_bias, num_visible, nrow(hidden))) + bias))
            }
            else
            {
                return(visible_form((hidden %*% weights) + visible_bias + bias))
            }
        }
        # train = function(X, learning_rate = 0.1) 
        # {
        #     'Replaces the range [i, j] of the
        #     object by value.
        #     '
        #     backup <- list(i, j, data[i,j])
        #     data[i,j] <<- value
        #     edits <<- c(edits, list(backup))
        #     invisible(value)
        # },
        # undo = function() 
        # {
        #     'Undoes the last edit() operation
        #     and update the edits field accordingly.
        #     '
        #     prev <- edits
        #     if(length(prev)) prev <- prev[[length(prev)]]
        #     else stop("No more edits to undo")
        #     edit(prev[[1]], prev[[2]], prev[[3]])
        #     ## trim the edits list
        #     length(edits) <<- length(edits) - 2
        #     invisible(prev)
        # },
        # show = function() 
        # {
        #     'Method for automatically printing matrix editors'
        #     cat("Reference matrix editor object of class",
        #     classLabel(class(.self)), "\n")
        #     cat("Data: \n")
        #     methods::show(data)
        #     cat("Undo list is of length", length(edits), "\n")
        # }
    )
)