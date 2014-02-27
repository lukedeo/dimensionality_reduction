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
            if(!is.null(nrow(visible)))
            {
                return(t(sigmoid(weights %*% t(visible))) + t(matrix(hidden_bias, num_hidden, nrow(visible))) + bias)
            }
            else
            {
                return(t(sigmoid(weights %*% (visible))) + hidden_bias + bias)    
            } #transpose (technically) here ~~~~^
            
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
        calculate_gradients = function(visible_batch, sparsity = NULL)
        {
            gradients <- list()
            pass_between(visible_batch, 2)

            v_0 <- sampling_sequence[[1]]$visible
            h_0 <- sampling_sequence[[1]]$hidden

            v_1 <- sampling_sequence[[2]]$visible
            h_1 <- sampling_sequence[[2]]$hidden

            n_samples <- ifelse(is.matrix(visible_batch), nrow(visible_batch), 1)

            gradients$grad_W <- ((t(h_0) %*% v_0) - (t(h_1) %*% v_1)) / n_samples
            gradients$grad_v <- colMeans(v_0 - v_1)
            if(is.null(sparsity))
            {
                gradients$grad_h <- colMeans(h_0 - h_1)
            }
            else
            {
                gradients$grad_h <- sparsity - colMeans(h_0)
            }
            return(gradients)
        },
        apply_gradients = function(gradients, learning_rate = 0.1, momentum = 0.0, l2_regularizer = 0)
        {
            gradients$grad_W <- gradients$grad_W * (1 - momentum)
            gradients$grad_h <- gradients$grad_h * (1 - momentum)
            gradients$grad_v <- gradients$grad_v * (1 - momentum)





            gradients$grad_W <- gradients$grad_W + momentum * (gradients$grad_W - l2_regularizer * weights)

            # don't apply l2 regularizer to biases
            gradients$grad_h <- gradients$grad_h + momentum * (gradients$grad_h)
            gradients$grad_v <- gradients$grad_v + momentum * (gradients$grad_v)

            weights <<- weights + learning_rate * gradients$grad_W
            hidden_bias <<- hidden_bias + learning_rate * gradients$grad_h
            visible_bias <<- visible_bias + learning_rate * gradients$grad_v
        },
        learn = function(X, learning_rate = 0.1, momentum = 0, l2_regularizer = 0, sparsity = NULL)
        {
            apply_gradients(calculate_gradients(X, sparsity), learning_rate, momentum, l2_regularizer)
        }
    )
)