
load_mnist()

sel <- sample(1:train$n, 1000, replace=F)

X <- train$x / max(train$x)

mnist.dae <- learn_denoising_autoencoder(X=X[sel, ], n_hidden=500, activation_type="sigmoid", epochs=1, batch = 1)
mnist.dae.small <- learn_denoising_autoencoder(X=X[sel, ], n_hidden=10, activation_type="sigmoid", epochs=4, batch = 1)

mnist.DAE <- autoencoder(X=X[sel, ], n_hidden=500, 
                         decoder="sigmoid", variant="denoising")


samp <- c(sample(which(train$y == 0), 1),
          sample(which(train$y == 1), 1),
          sample(which(train$y == 2), 1),
          sample(which(train$y == 3), 1),
          sample(which(train$y == 4), 1))

samp <- c(sample(which(train$y == 5), 1),
          sample(which(train$y == 6), 1),
          sample(which(train$y == 7), 1),
          sample(which(train$y == 8), 1),
          sample(which(train$y == 9), 1))

orig <- X[samp, ]
reco <- reconstruct_autoencoder(mnist.DAE, X[samp, ])

grid.arrange(show_digit(orig[1, ]),
             show_digit(orig[2, ]),
             show_digit(orig[3, ]),
             show_digit(orig[4, ]),
             show_digit(orig[5, ]),
             show_digit(reco[1, ]),
             show_digit(reco[2, ]),
             show_digit(reco[3, ]),
             show_digit(reco[4, ]),
             show_digit(reco[5, ]),
             nrow=2)




show_digit((mnist.dae$encoder$W[idx[1], ]))

filt <- abs(matrix(mnist.dae$encoder$W[idx[1], ], 28, 28)[,28:1])
filt <- filt / sqrt(sum(filt^2))

dataL <-melt(filt)

p <- ggplot(dataL, aes(x=Var1, y=Var2, fill=value))+geom_tile()+scale_fill_gradient(low="white",high="black")
p <- p + theme_bw()+theme(line = element_blank(),
                          text = element_blank(),
                          line = element_blank(),
                          title = element_blank(),
                          plot.background = element_blank(),
                          panel.border = element_blank(),
                          legend.position = "none")
p


idx = sort.int(apply(X=mnist.dae$encoder$W, MARGIN=1, function(x){mean(x)}), decreasing=F, index.return=T)$ix

idx <- sample(1:600, 9, replace=F)

grid.arrange(show_digit((DAE[[1]]$encoder$W[idx[1], ]), T),
             show_digit((DAE[[1]]$encoder$W[idx[2], ]), T),
             show_digit((DAE[[1]]$encoder$W[idx[3], ]), T),
             show_digit((DAE[[1]]$encoder$W[idx[4], ]), T),
             show_digit((DAE[[1]]$encoder$W[idx[5], ]), T),
             show_digit((DAE[[1]]$encoder$W[idx[6], ]), T),
             show_digit((DAE[[1]]$encoder$W[idx[7], ]), T),
             show_digit((DAE[[1]]$encoder$W[idx[8], ]), T),
             show_digit((DAE[[1]]$encoder$W[idx[9], ]), T),
             nrow=3)





sel <- sample(1:train$n, 40000, replace=F)

digit.ae <- autoencoder(X[sel, ], n_hidden=c(400, 300, 100, 80, 40, 26, 15, 7, 2), variant="denoising", epochs=1, batch=1,
            decoder=c(rep("sigmoid", 8), "linear"))

encoded_digits <- reconstruct(digit.ae, X[sel, ])

encoded_numbers <- data.frame(x=encoded_digits[, 1], y=encoded_digits[, 2], Number=factor(train$y[sel]))
p <- ggplot(encoded_numbers, aes(x=x, y=y, color=Number))
p + geom_point(aes(color = Number), alpha = 0.4) + xlab("Dimension 1") + ylab("Dimension 2") + 
    labs(title="Unsupervised Feature Representation\nof MNIST, 9 Layers.") +
    guides(color = guide_legend(override.aes= list(alpha = 1, size = 4.2)))
