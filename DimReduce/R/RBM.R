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
    if(is.matrix(x))
    {
        return(matrix(runif(length(x)), nrow(x), ncol(x)) < x)
    }
    return(runif(length(x), nrow(x), ncol(x)) < x)   
}

RBM <- setRefClass("RBM", 
    fields = list(
        weights = "matrix", 
        hidden_bias = "matrix", 
        visible_bias = "matrix", 
        visible_form = "function",
        num_hidden = "numeric",
        num_visible = "numeric",
        sampling_sequence = "list"
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
                return(t(sigmoid(weights %*% (visible))) + hidden_bias + bias)    
            }  #transpose (technically) here ~~~~^
            
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
        },
        pass_between = function(visible, num_passes = 1)
        {
            passes <- 1
            sampling_sequence <<- list()
            while(passes <= num_passes)
            {
                hidden <- expected_hidden(visible)
                sampling_sequence[[passes]] <<- list(hidden = hidden, visible = visible)
                visible <- expected_visible(bernoulli(hidden))
                passes <- passes + 1
            }
        },
        reconstruct = function(visible, num_passes = 1)
        {
            pass_between(visible, num_passes)
            return(sampling_sequence[[num_passes]])
        },
        calculate_gradients = function(visible_batch)
        {
            gradients <- list()
            pass_between(visible_batch, 2)

            v_0 <- sampling_sequence[[1]]$visible
            h_0 <- sampling_sequence[[1]]$hidden

            v_1 <- sampling_sequence[[2]]$visible
            h_1 <- sampling_sequence[[2]]$hidden

            gradients$grad_W <- ((t(h_0) %*% v_0) - (t(h_1) %*% v_1)) / length(visible_batch[, 1])
            gradients$grad_v <- colMeans(v_0 - v_1)
            gradients$grad_h <- colMeans(h_0 - h_1)
            return(gradients)
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