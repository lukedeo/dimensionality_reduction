autoencoder <- function(X, ...) UseMethod("autoencoder")

autoencoder.default <- function(X, n_hidden, decoder = "linear", 
	variant = "denoising", epochs = 1, batch = 1, learning = 0.02, 
	momentum = 0.9, regularization = 0.001, noise = 0.02, optimizer = "sgd")
{
	L <- NULL
	if(length(n_hidden) == 1)
	{
		if (variant == "denoising")
		{
			cat("Learning a denoising autoencoder...")
			L <- learn_denoising_autoencoder(X, n_hidden, decoder[1], epochs, batch, learning, momentum, regularization, noise)
		}
		else if (variant == "normal")
		{
			if (optimizer == "sgd")
			{
				cat("Learning a standard autoencoder trained with SGD...")
				L <- learn_autoencoder(X, n_hidden, decoder[1], epochs, batch, learning, momentum, regularization)
			}
			else if (optimizer == "bfgs")
			{
				cat("Learning a standard autoencoder trained with S-LM-BFGS...")
				L <- learn_bfgs_autoencoder(X, n_hidden, decoder[1], epochs, batch, learning, momentum, regularization)
			}
			else
			{
				stop("Unrecognized optimizer, needs to be one of \"sgd\" or \"bfgs\"")
			}
		}
		else
		{
			stop("Unrecognized variant, needs to be one of \"normal\" or \"denoising\"")
		}
		class(L) <- "autoencoder"
	}
	else
	{
		decoders <- decoder
		if(length(decoder) == 1)
		{
			decoders <- c(rep(decoder, length(n_hidden)))
		}
		if(length(decoders) != length(n_hidden))
		{
			stop("vectors of hidden units and of types of units dont match in length.")
		}
		if(variant == "denoising")
		{
			cat("Learning a deep denoising autoencoder...\n")
			L <- learn_deep_autoencoder(X, n_hidden, decoders, epochs, batch, learning, 
			momentum, regularization, TRUE, noise)
		}
		else
		{
			cat("Learning a deep autoencoder...\n")
			L <- learn_deep_autoencoder(X, n_hidden, decoders, epochs, batch, learning, 
			momentum, regularization, FALSE, noise)
		}

	
		class(L) <- c("autoencoder", "deep_autoencoder")
	}
	L
}

update.autoencoder <- function(autoenc, X, variant = "denoising", epochs = 1, batch = 1, learning = 0.02, 
	momentum = 0.9, regularization = 0.001, noise = 0.02, optimizer = "sgd")
{
	L <- NULL
	if("deep_autoencoder" %in% class(autoenc))
	{
		denoising <- ifelse(variant == "denoising", TRUE, FALSE)
		L <- continue_learn_deep_autoencoder(autoenc, X, epochs, batch, learning, 
			momentum, regularization, denoising, noise)
		class(L) <- c("autoencoder", "deep_autoencoder")
	}
	else
	{
		if (variant == "denoising")
		{
			cat("Updating a denoising autoencoder...")
			L <- continue_learn_denoising_autoencoder(autoenc, X, epochs, batch, learning, momentum, regularization, noise)
		}
		else if (variant == "normal")
		{
			if (optimizer == "sgd")
			{
				cat("Updating a standard autoencoder trained with SGD...")
				L <- continue_learn_autoencoder(autoenc, X, epochs, batch, learning, momentum, regularization)
			}
			else if (optimizer == "bfgs")
			{
				cat("Updating a standard autoencoder trained with S-LM-BFGS...")
				L <- continue_learn_bfgs_autoencoder(autoenc, X, epochs, batch, learning, momentum, regularization)
			}
			else
			{
				stop("Unrecognized optimizer, needs to be one of \"sgd\" or \"bfgs\"")
			}
		}
		else
		{
			stop("Unrecognized variant, needs to be one of \"normal\" or \"denoising\"")
		}
		class(L) <- "autoencoder"
	}	
	L	
}

encode <- function(autoenc, X)
{
	if("deep_autoencoder" %in% class(autoenc))
	{
		return(encode_deep_autoencoder(autoenc, X))
	}
	return(encode_autoencoder(autoenc, X))
}

decode <- function(autoenc, Y)
{
	if("deep_autoencoder" %in% class(autoenc))
	{
		return(decode_deep_autoencoder(autoenc, Y))
	}
	return(decode_autoencoder(autoenc, Y))
}

reconstruct <- function(autoenc, X)
{
	if("deep_autoencoder" %in% class(autoenc))
	{
		return(reconstruct_deep_autoencoder(autoenc, X))
	}
	return(reconstruct_autoencoder(autoenc, X))
}

