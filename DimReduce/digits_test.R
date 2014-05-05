sigmoid <- function(x)
{
    return(1 / (1 + exp(-x)))
}

plot_image <- function(mat)
{
    dataL <- melt(mat)
    
    p <- ggplot(dataL, aes(x=Var1, y=Var2, fill=value))+geom_tile()+scale_fill_gradient(low = "white", high="black")
    p <- p + theme_bw()+theme(line = element_blank(),
                         text = element_blank(),
                         line = element_blank(),
                         title = element_blank(),
                         plot.background = element_blank(),
                         panel.border = element_blank(),
                         legend.position = "none")
    return(p)
}

p <- ggplot()

library(reshape2)

rotate <- function(x) t(apply(x, 2, rev))

im <- rotate(matrix(X[52, ], 16, 16))



dataL <- melt(im)


p <- ggplot(dataL, aes(x=Var1, y=Var2, fill=value))+geom_tile()+scale_fill_gradient(low="white", high="black")
p + theme_bw()+theme(line = element_blank(),
          text = element_blank(),
          line = element_blank(),
          title = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "none")






labels <- cbind(rep("0", 1100), rep("1", 1100), rep("2", 1100), rep("3", 1100)
                , rep("4", 1100), rep("5", 1100), rep("6", 1100), rep("7", 1100)
                , rep("8", 1100), rep("9", 1100))
color_vec <- cbind(rep("blue", 1100), rep("red", 1100), rep("green", 1100), rep("black", 1100)
                   , rep("purple", 1100), rep("grey", 1100), rep("pink", 1100), rep("yellow", 1100)
                   , rep("orange", 1100), rep("brown", 1100))
N <- 2200
samp <- (sample(x=1:nrow(digits), size=N, replace=F))
X <- (digits[samp, ]) / max(digits[samp, ])






digit.autoencoder <- learn_denoiding_autoencoder(X, 100)

digit.autoencoder$encoder$W[1, ]


X_tilde <- reconstruct_autoencoder(autoenc=digit.autoencoder_1, X=X)


im <- rotate(matrix(X_tilde[3, ], 16, 16))




dataL <- melt(im)


p <- ggplot(dataL, aes(x=Var1, y=Var2, fill=value))+geom_tile()+scale_fill_gradient(low = "white", high="black")
p + theme_bw()+theme(line = element_blank(),
                     text = element_blank(),
                     line = element_blank(),
                     title = element_blank(),
                     plot.background = element_blank(),
                     panel.border = element_blank(),
                     legend.position = "none")








data(binaryalphadigits)
digits <- binaryalphadigits
data(faces)
faces <- t(faces)

faces <- faces / max(faces)

faces <- 1 - faces


dim(digits)



data(frey)

frey <- t(frey) 
frey <- frey / max(frey)

frey <- 1 - frey
plot_image(rotate(rotate(matrix(frey[1, ], 20, 28))))

face.autoencoder <- learn_denoising_autoencoder(X=frey[sample(1:nrow(frey), size=1000, replace=T), ], n_hidden=200, "sigmoid", epochs=10, batch=2)
face.autoencoder.std <- learn_autoencoder(X=frey[sample(1:nrow(frey), size=1000, replace=T), ], n_hidden=200, epochs=1, batch=5)
face.autoencoder.bfgs <- learn_bfgs_autoencoder(X=frey[sample(1:nrow(frey), size=1000, replace=T), ], n_hidden=200, epochs=1, batch=3)

faces_pred <- reconstruct_autoencoder(autoenc=face.autoencoder.std, X=frey)



require(gridExtra)

for(i in sample(1:nrow(frey), size=20, replace=F))
{
    act <- plot_image(rotate(rotate(matrix(frey[i, ], 20, 28))))
    rec <- plot_image(rotate(rotate(matrix(faces_pred[i, ], 20, 28))))
    grid.arrange(act, rec, ncol=2)
}

sel <- sample(1:nrow(frey), size=5, replace=F)
sel <- c(500, 400, 800, 1200, 1900)
grid.arrange(plot_image(rotate(rotate(matrix(frey[sel[1], ], 20, 28)))),
             plot_image(rotate(rotate(matrix(frey[sel[2], ], 20, 28)))),
             plot_image(rotate(rotate(matrix(frey[sel[3], ], 20, 28)))),
             plot_image(rotate(rotate(matrix(frey[sel[4], ], 20, 28)))),
             plot_image(rotate(rotate(matrix(frey[sel[5], ], 20, 28)))),
             plot_image(rotate(rotate(matrix(faces_pred[sel[1], ], 20, 28)))),
             plot_image(rotate(rotate(matrix(faces_pred[sel[2], ], 20, 28)))),
             plot_image(rotate(rotate(matrix(faces_pred[sel[3], ], 20, 28)))),
             plot_image(rotate(rotate(matrix(faces_pred[sel[4], ], 20, 28)))),
             plot_image(rotate(rotate(matrix(faces_pred[sel[5], ], 20, 28)))),
             nrow=2)

11.69

filt <- (rotate(rotate(matrix((face.autoencoder$encoder$W[50, ]), 20, 28))))

filt[filt > 0.2] <- 0.2
filt[filt < -0.2] <- -0.2

dataL <- melt(filt)


p <- ggplot(dataL, aes(x=Var1, y=Var2, fill=value))+geom_tile()+scale_fill_gradient(low="white", high="black")
p <- p + theme_bw()+theme(line = element_blank(),
                          text = element_blank(),
                          line = element_blank(),
                          title = element_blank(),
                          plot.background = element_blank(),
                          panel.border = element_blank(),
                          legend.position = "none")
p


idx = sort.int(apply(X=face.autoencoder$encoder$W, MARGIN=1, function(x){sum(x^2)}), decreasing=T, index.return=T)$ix

idx <- sample(1:200, 16, replace=F)

face.autoencoder <- face.autoencoder.std
grid.arrange(show_filter((face.autoencoder$encoder$W[idx[1], ])),
             show_filter((face.autoencoder$encoder$W[idx[2], ])),
             show_filter((face.autoencoder$encoder$W[idx[3], ])),
             show_filter((face.autoencoder$encoder$W[idx[4], ])),
             show_filter((face.autoencoder$encoder$W[idx[5], ])),
             show_filter((face.autoencoder$encoder$W[idx[6], ])),
             show_filter((face.autoencoder$encoder$W[idx[7], ])),
             show_filter((face.autoencoder$encoder$W[idx[8], ])),
             show_filter((face.autoencoder$encoder$W[idx[9], ])),
             show_filter((face.autoencoder$encoder$W[idx[10], ])),
             show_filter((face.autoencoder$encoder$W[idx[11], ])),
             show_filter((face.autoencoder$encoder$W[idx[12], ])),
             show_filter((face.autoencoder$encoder$W[idx[13], ])),
             show_filter((face.autoencoder$encoder$W[idx[14], ])),
             show_filter((face.autoencoder$encoder$W[idx[15], ])),
             show_filter((face.autoencoder$encoder$W[idx[16], ])),
             nrow=4)








data(binaryalphadigits)
digits <- binaryalphadigits


digits <- as.matrix(digits / max(digits))


dataL <- melt(matrix(as.numeric(digits[1, ]), 20, 16))


p <- ggplot(dataL, aes(x=Var1, y=Var2, fill=value))+geom_tile()
p <- p + theme_bw()+theme(line = element_blank(),
                          text = element_blank(),
                          line = element_blank(),
                          title = element_blank(),
                          plot.background = element_blank(),
                          panel.border = element_blank(),
                          legend.position = "none")
p


plot_image((rotate(matrix((digits[40, ]), 20, 16))))

char.autoencoder <- learn_denoising_autoencoder(X=digits[sample(1:nrow(digits), size=1000, replace=T), ], "sigmoid", n_hidden=200, epochs=10, batch=1)
char.autoencoder.std <- learn_autoencoder(X=digits[sample(1:nrow(digits), size=1000, replace=T), ], n_hidden=200, epochs=20, batch=100)

char_pred <- reconstruct_autoencoder(autoenc=char.autoencoder.std, X=digits)

sellet <- sample((10*39+1):nrow(digits), size=5, replace=F)

selnum <- sample(1:(10*39+1), size=5, replace=F)
sel <- sellet
sel <- selnum
grid.arrange(plot_image((rotate(matrix(digits[sel[1], ], 20, 16)))),
             plot_image((rotate(matrix(digits[sel[2], ], 20, 16)))),
             plot_image((rotate(matrix(digits[sel[3], ], 20, 16)))),
             plot_image((rotate(matrix(digits[sel[4], ], 20, 16)))),
             plot_image((rotate(matrix(digits[sel[5], ], 20, 16)))),
             plot_image((rotate(matrix(char_pred[sel[1], ], 20, 16)))),
             plot_image((rotate(matrix(char_pred[sel[2], ], 20, 16)))),
             plot_image((rotate(matrix(char_pred[sel[3], ], 20, 16)))),
             plot_image((rotate(matrix(char_pred[sel[4], ], 20, 16)))),
             plot_image((rotate(matrix(char_pred[sel[5], ], 20, 16)))),
             nrow=2)


char.autoencoder <- char.autoencoder.std
idx <- sample(1:100, 16, replace=F)
grid.arrange(show_filter((char.autoencoder$encoder$W[idx[1], ])),
             show_filter((char.autoencoder$encoder$W[idx[2], ])),
             show_filter((char.autoencoder$encoder$W[idx[3], ])),
             show_filter((char.autoencoder$encoder$W[idx[4], ])),
             show_filter((char.autoencoder$encoder$W[idx[5], ])),
             show_filter((char.autoencoder$encoder$W[idx[6], ])),
             show_filter((char.autoencoder$encoder$W[idx[7], ])),
             show_filter((char.autoencoder$encoder$W[idx[8], ])),
             show_filter((char.autoencoder$encoder$W[idx[9], ])),
             show_filter((char.autoencoder$encoder$W[idx[10], ])),
             show_filter((char.autoencoder$encoder$W[idx[11], ])),
             show_filter((char.autoencoder$encoder$W[idx[12], ])),
             show_filter((char.autoencoder$encoder$W[idx[13], ])),
             show_filter((char.autoencoder$encoder$W[idx[14], ])),
             show_filter((char.autoencoder$encoder$W[idx[15], ])),
             show_filter((char.autoencoder$encoder$W[idx[16], ])),
             nrow=4)

labels <- cbind(rep("0", 39), rep("1", 39), rep("2", 39), rep("3", 39), rep("4", 39), 
                rep("5", 39), rep("6", 39), rep("7", 39), rep("8", 39), rep("9", 39),rep("A", 39),
                rep("B", 39),rep("C", 39),rep("D", 39),rep("E", 39),rep("F", 39),rep("G", 39),rep("H", 39),
                rep("I", 39),rep("J", 39),rep("K", 39),rep("L", 39),rep("M", 39),rep("N", 39),rep("O", 39),
                rep("P", 39),rep("Q", 39),rep("R", 39),rep("S", 39),rep("T", 39),rep("U", 39),rep("V", 39),
                rep("W", 39),rep("X", 39),rep("Y", 39),rep("Z", 39))


char.autoencoder <- learn_denoising_autoencoder(X=digits[sample(1:nrow(digits), size=1000, replace=T), ], n_hidden=100, epochs=4, batch=5)




library(kernlab)

sel <- sample(1:nrow(digits), size=1000, replace=T)
alphanum_model <- ksvm(x=digits[sel, ], y=factor(labels[sel]), cross = 5)

pred_class <- predict(alphanum_model, digits[-sel, ])


sum(pred_class == labels[-sel]) / length(pred_class)




sel <- sample(1:nrow(digits), size=1000, replace=T)
char.autoencoder <- learn_denoising_autoencoder(X=digits[sel, ], n_hidden=100, epochs=4, batch=5)
encoded_digits <- encode_autoencoder(char.autoencoder, digits[sel, ])
char.autoencoder_2 <- learn_denoising_autoencoder(X=encoded_digits, n_hidden=50, epochs=4, batch=5)

encoded_digits <- encode_autoencoder(char.autoencoder, digits[sel, ])
encoded_digits <- encode_autoencoder(char.autoencoder_2, encoded_digits)


alphanum_model_DL <- ksvm(x=encoded_digits, y=factor(labels[sel]), cross = 5)




encoded_digits <- encode_autoencoder(char.autoencoder, digits[-sel, ])
encoded_digits <- encode_autoencoder(char.autoencoder_2, encoded_digits)

pred_class <- predict(alphanum_model_DL, encoded_digits)


sum(pred_class == labels[-sel]) / length(pred_class)






# Reducing dimension ------------------------------------------------------

data(digits)
digits <- t(digits)
digits <- digits / max(digits)
labels <- cbind(rep("0", 1100), rep("1", 1100), rep("2", 1100), rep("3", 1100)
                , rep("4", 1100), rep("5", 1100), rep("6", 1100), rep("7", 1100)
                , rep("8", 1100), rep("9", 1100))
color_vec <- cbind(rep("blue", 1100), rep("red", 1100), rep("green", 1100), rep("black", 1100)
                   , rep("purple", 1100), rep("grey", 1100), rep("pink", 1100), rep("yellow", 1100)
                   , rep("orange", 1100), rep("brown", 1100))
N <- 3000
samp <- (sample(x=1:nrow(digits), size=N, replace=F))
X <- (digits[samp, ]) / max(digits[samp, ])

digit.autoencoder_1 <- learn_denoising_autoencoder(X, 150, "sigmoid", epochs=20, batch=1)
encoded_digits <- encode_autoencoder(digit.autoencoder_1, X)

digit.autoencoder_2 <- learn_denoising_autoencoder(encoded_digits, 120, "sigmoid", epochs=20, batch=1)
encoded_digits <- encode_autoencoder(digit.autoencoder_2, encoded_digits)

digit.autoencoder_3 <- learn_denoising_autoencoder(encoded_digits, 80, "sigmoid", epochs=12, batch=1)
encoded_digits <- encode_autoencoder(digit.autoencoder_3, encoded_digits)

digit.autoencoder_4 <- learn_denoising_autoencoder(encoded_digits, 50, "sigmoid", epochs=12, batch=1)
encoded_digits <- encode_autoencoder(digit.autoencoder_4, encoded_digits)

digit.autoencoder_5 <- learn_denoising_autoencoder(encoded_digits, 40, "sigmoid", epochs=12, batch=1)
encoded_digits <- encode_autoencoder(digit.autoencoder_5, encoded_digits)

digit.autoencoder_6 <- learn_denoising_autoencoder(encoded_digits, 26, "sigmoid", epochs=12, batch=1)
encoded_digits <- encode_autoencoder(digit.autoencoder_6, encoded_digits)

digit.autoencoder_7 <- learn_denoising_autoencoder(encoded_digits, 15, "sigmoid", epochs=12, batch=1)
encoded_digits <- encode_autoencoder(digit.autoencoder_7, encoded_digits)

digit.autoencoder_8 <- learn_denoising_autoencoder(encoded_digits, 7, "sigmoid", epochs=12, batch=1)
encoded_digits <- encode_autoencoder(digit.autoencoder_8, encoded_digits)

digit.autoencoder_9 <- learn_denoising_autoencoder(encoded_digits, 4, "sigmoid", epochs=50, batch=1)
encoded_digits <- encode_autoencoder(digit.autoencoder_9, encoded_digits)

digit.autoencoder_10 <- learn_denoising_autoencoder(encoded_digits, 2, "linear", epochs=50, batch=1)
encoded_digits <- encode_autoencoder(digit.autoencoder_10, encoded_digits)


all <- (digits) / max(digits)

encoded_digits <- encode_autoencoder(digit.autoencoder_1, all)
encoded_digits <- encode_autoencoder(digit.autoencoder_1, X)
encoded_digits <- encode_autoencoder(digit.autoencoder_2, encoded_digits)
encoded_digits <- encode_autoencoder(digit.autoencoder_3, encoded_digits)
encoded_digits <- encode_autoencoder(digit.autoencoder_4, encoded_digits)
encoded_digits <- encode_autoencoder(digit.autoencoder_5, encoded_digits)
encoded_digits <- encode_autoencoder(digit.autoencoder_6, encoded_digits)
encoded_digits <- encode_autoencoder(digit.autoencoder_7, encoded_digits)
encoded_digits <- encode_autoencoder(digit.autoencoder_8, encoded_digits)
encoded_digits <- encode_autoencoder(digit.autoencoder_9, encoded_digits)
encoded_digits <- encode_autoencoder(digit.autoencoder_10, encoded_digits)

qplot(x=encoded_digits[, 1], y=encoded_digits[, 2], color=factor(labels))

encoded_numbers <- data.frame(x=encoded_digits[, 1], y=encoded_digits[, 2], Number=factor(labels))
p <- ggplot(encoded_numbers, aes(x=x, y=y, color=Number))
p + geom_point(aes(color = Number), alpha = 0.6) + xlab("Dimension 1") + ylab("Dimension 2") + 
    labs(title="Unsupervised Feature Representation\nof USPS Handwritten Digits, 10 Layers.") +
    guides(color = guide_legend(override.aes= list(alpha = 1, size = 4.2)))

dim(X)

linemb <- local_linear_embedding(X = X + matrix(rnorm(n=1000*256, mean=0, sd=0.001), 1000, 256), k = 10, d = 2, verbose = TRUE)

encoded_numbers <- data.frame(x=linemb$Y[, 1], y=linemb$Y[, 2], Number=factor(labels[samp]))
p <- ggplot(encoded_numbers, aes(x=x, y=y, color=Number))
p + geom_text(aes(color = Number, label=labels[samp])) + xlab("Dimension 1") + ylab("Dimension 2") + 
    labs(title="Unsupervised Feature Representation\nof USPS Handwritten Digits, LLE.") +
    guides(color = guide_legend(override.aes= list(alpha = 1, size = 4.2)))

