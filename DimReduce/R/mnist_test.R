library(ggplot2)

p <- ggplot(data.frame(matrix(MNIST[1], 28, 28)))


image((matrix(MNIST[1], 28, 28)))

mnist <- as.matrix(MNIST)

im <- as.matrix((matrix(mnist[1, ], 28, 28)))


image(im)

library(ggplot2)
library(ggthemes)
library(reshape2)
library(scales)
library(RColorBrewer)



# faces stuff -------------------------------------------------------------


data(faces)

faces <- as.matrix(t(faces))

faces <- faces / max(faces)

face_gen <- function(idx, data = faces)
{
    im <- matrix(0, 64, 64)
    for(i in 0:63)
    {
        for(j in 0:63)
        {
            im[i+1, 64 - j] <- data[idx, (i *64 + (j)) + 1 ]
        }
    }
    return(im)
}

a <- faces_rbm$expected_visible(faces_rbm$expected_hidden(faces[2, ]))

a <- faces_rbm$reconstruct(faces[2, ], 20)$visible

orig <- ggplot(melt(face_gen(1, a)), aes(x=Var1,y=Var2))
orig <- orig + geom_tile(aes(fill=value)) + 
    scale_fill_gradient(low = "black", high = "white") + 
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())  
orig <- orig + labs(title = "Original Face")
orig

for(i in seq(1, 200, by=10))
{
    orig <- ggplot(melt(face_gen(i, data = faces_rbm$reconstruct(faces[i:(i+1), ], 3)$visible)), aes(x=Var1,y=Var2))
    orig <- orig + geom_tile(aes(fill=value)) + 
        scale_fill_gradient(low = "black", high = "white") + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank())  
    orig <- orig + labs(title = "Original Face")
    multiplot(orig)
}

faces_rbm <- RBM(64^2, 4000, TRUE)

faces_rbm$learn(faces[1:5, ])

for(i in 1:(nrow(faces)))
{
    faces_rbm$learn(faces[sample(x=c(1:nrow(faces)), 15, replace=T), ], 0.06, momentum = 0.1)
    if(i %% 5 == 0){
        cat(sprintf("%i is done.\n", i))
    }
}


reconstructions <- faces_rbm$reconstruct(faces[c(1:100), ], 2)
reconstructions <- faces_rbm$expected_visible(faces_rbm$expected_hidden(faces[1:100, ])))

for(i in c(1:30)))
{
    orig <- ggplot(melt((face_gen(i))), aes(x=Var1,y=Var2))
    orig <- orig + geom_tile(aes(fill=value)) + 
        scale_fill_gradient(low = "black", high = "white") + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank())  
    orig <- orig + labs(title = "Original Face")
    
    learned <- ggplot(melt((face_gen(i, data = reconstructions$visible))), aes(x=Var1,y=Var2))
    learned <- learned + geom_tile(aes(fill=value))+ 
        scale_fill_gradient(low = "black", high = "white") + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank())
    
    learned <- learned + labs(title = "Reconstructed/Learned Face")
    
    multiplot(orig, learned, cols=2)
}






# mnist stuff -------------------------------------------------------------


image_data <- data.frame(x.idx = rep(c(1:28), 28), y.idx = rep(c(1:28), 28), value = rep(0, 784))
ctr <- 1
for(i in 1:28)
{
    for(j in 1:28)
    {
        image_data$value[ctr] <- im[i, j]
        image_data$x.idx[ctr] <- i
        image_data$y.idx[ctr] <- j
        ctr <- ctr + 1
    }
}


im <- matrix(t(t(mnist[2, ])), 28, 28)

image_gen <- function(idx, data = mnist)
{
    im <- matrix(0, 28, 28)
    for(i in 0:27)
    {
        for(j in 0:27)
        {
            im[i+1, 28 - j] <- data[idx, (i *28 + (j)) + 1 ]
        }
    }
    return(im)
}

mnist_bin <- ifelse(mnist > 0.2, 1, 0)
myrbm <- RBM(28^2, 600, binary = TRUE)

for(i in 1:(nrow(mnist)))
{
    myrbm$learn(mnist_bin[sample(1:nrow(mnist), 5, replace=T), ], 0.2, 0.1, 0.0001)
    if(i %% 10 == 0)
    {
        cat(sprintf("%i is done.\n", i))
    }
}
myrbm_2 <- RBM(600, 200, binary = T)
layer_2 <- (myrbm$expected_hidden(mnist_bin))

for(i in 1:(nrow(mnist)))
{
    myrbm_2$learn(layer_2[sample(1:nrow(mnist), 5, replace=T), ], 0.2, 0.1, 0.0001)
    if(i %% 10 == 0)
    {
        cat(sprintf("%i is done.\n", i))
    }
}

myrbm_3 <- RBM(200, 10, binary = T)
layer_3 <- (myrbm_2$expected_hidden(layer_2))

for(i in 1:(nrow(mnist)))
{
    myrbm_3$learn(layer_3[sample(1:nrow(mnist), 5, replace=T), ], 0.2, 0.1, 0.0001)
    if(i %% 10 == 0)
    {
        cat(sprintf("%i is done.\n", i))
    }
}


twolayers <- (myrbm_2$reconstruct(myrbm$expected_hidden(mnist_bin), 2))

reconstructions <- myrbm$reconstruct(mnist_bin, 2)
reconstructions <- myrbm$reconstruct(mnist_bin, 2)
reconstructions <- list(visible = myrbm$expected_visible(myrbm$expected_hidden(mnist)))

for(i in seq(1, 600, 20))
{
    orig <- ggplot(melt((image_gen(i, data = mnist))), aes(x=Var1,y=Var2))
    orig <- orig + geom_tile(aes(fill=value)) + 
                   scale_fill_gradient(low = "white", high = "black") + 
                   theme(axis.line=element_blank(),
                         axis.text.x=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks=element_blank(),
                         axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         legend.position="none",
                         panel.background=element_blank(),panel.border=element_blank(),
                         panel.grid.major=element_blank(),
                         panel.grid.minor=element_blank(),
                         plot.background=element_blank())  
    orig <- orig + labs(title = "Original Digit")
    
    learned <- ggplot(melt((image_gen(i, data = reconstructions$visible))), aes(x=Var1,y=Var2))
    learned <- learned + geom_tile(aes(fill=value))+ 
                         scale_fill_gradient(low = "white", high = "black") + 
                         theme(axis.line=element_blank(),
                               axis.text.x=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks=element_blank(),
                               axis.title.x=element_blank(),
                               axis.title.y=element_blank(),
                               legend.position="none",
                               panel.background=element_blank(),
                               panel.border=element_blank(),
                               panel.grid.major=element_blank(),
                               panel.grid.minor=element_blank(),
                               plot.background=element_blank())
    
    learned <- learned + labs(title = "Reconstructed/Learned Digit")
    
    multiplot(orig, learned, cols=2)
}


randomat <-  myrbm_3$expected_hidden(myrbm_2$expected_hidden(myrbm$expected_hidden(mnist_bin[100, ])))


randim <- myrbm$expected_visible(myrbm_2$expected_visible(myrbm_3$expected_visible(randomat)))
orig <- ggplot(melt((image_gen(100, data = mnist))), aes(x=Var1,y=Var2))
orig <- orig + geom_tile(aes(fill=value)) + 
    scale_fill_gradient(low = "white", high = "black") + 
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())  
orig <- orig + labs(title = "Random Sampled Digit")
orig


randomat <-  myrbm_3$expected_hidden(myrbm_2$expected_hidden(myrbm$expected_hidden(mnist_bin)))


randim <- myrbm$expected_visible(myrbm_2$expected_visible(myrbm_3$expected_visible(randomat)))

for(i in seq(12, 60, 2))
{
    orig <- ggplot(melt((image_gen(i, data = mnist))), aes(x=Var1,y=Var2))
    orig <- orig + geom_tile(aes(fill=value)) + 
        scale_fill_gradient(low = "white", high = "black") + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank())  
    orig <- orig + labs(title = "Original Digit")
    
    learned <- ggplot(melt((image_gen(i, data = randim))), aes(x=Var1,y=Var2))
    learned <- learned + geom_tile(aes(fill=value))+ 
        scale_fill_gradient(low = "white", high = "black") + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank())
    
    learned <- learned + labs(title = "Reconstructed/Learned Digit")
    
    multiplot(orig, learned, cols=2)
}



# Digits dataset ----------------------------------------------------------

data(digits)
digits <- t(digits)

labels <- cbind(rep("1", 1100), rep("2", 1100), rep("3", 1100)
                , rep("4", 1100), rep("5", 1100), rep("6", 1100), rep("7", 1100)
                , rep("8", 1100), rep("9", 1100), rep("0", 1100))
color_vec <- cbind(rep("blue", 1100), rep("red", 1100), rep("green", 1100), rep("black", 1100)
                   , rep("purple", 1100), rep("grey", 1100), rep("pink", 1100), rep("yellow", 1100)
                   , rep("orange", 1100), rep("brown", 1100))

digits <- digits / max(digits)

N <- 600
samp <- sort(sample(x=1:nrow(digits), size=N, replace=F))
X <- (digits[samp, ]) / max(digits[samp, ])
labs <- labels[samp]

image_gen <- function(idx, data = digits)
{
    im <- matrix(0, 16, 16)
    for(i in 0:15)
    {
        for(j in 0:15)
        {
            im[i+1, 16 - j] <- data[idx, (i *16 + (j)) + 1 ]
        }
    }
    return(im)
}
orig <- ggplot(melt((image_gen(200, data = mnist_bin))), aes(x=Var1,y=Var2))
orig <- orig + geom_tile(aes(fill=value)) + 
    scale_fill_gradient(low = "white", high = "black") + 
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())  
orig <- orig + labs(title = "Original Digit")
orig

mnist_bin <- ifelse(mnist > 0.2, 1, 0)
mnist_bin <- mnist
myrbm <- RBM(16^2,150, binary = T)

for(k in 1:30){
for(i in 1:(nrow(mnist)))
{
    myrbm$learn(mnist_bin[sample(1:nrow(mnist), 5, replace=T), ], 0.2, 0.1, 0.0001)
    if(i %% 10 == 0)
    {
        cat(sprintf("%i is done.\n", i))
    }
}
}
myrbm_2 <- RBM(150, 70, binary = F)
layer_2 <- (myrbm$expected_hidden(mnist_bin))

for(k in 1:30){
for(i in 1:(nrow(mnist)))
{
    myrbm_2$learn(layer_2[sample(1:nrow(mnist), 5, replace=T), ], 0.2, 0.1, 0.0001)
    if(i %% 10 == 0)
    {
        cat(sprintf("%i is done.\n", i))
    }
}
}
myrbm_3 <- RBM(70, 20, binary = T)
layer_3 <- (myrbm_2$expected_hidden(layer_2))

for(k in 1:30){
for(i in 1:(nrow(mnist)))
{
    myrbm_3$learn(layer_3[sample(1:nrow(mnist), 5, replace=T), ], 0.2, 0.15, 0.00001)
    if(i %% 10 == 0)
    {
        cat(sprintf("%i is done.\n", i))
    }
}
}
myrbm_4 <- RBM(20, 10, binary = T)
layer_4 <- (myrbm_3$expected_hidden(layer_3))

for(k in 1:30){
for(i in 1:(nrow(mnist)))
{
    myrbm_4$learn(layer_4[sample(1:nrow(mnist), 5, replace=T), ], 0.1, 0.1, 0.0001)
    if(i %% 10 == 0)
    {
        cat(sprintf("%i is done.\n", i))
    }
}
}
myrbm_5 <- RBM(10, 2, binary = F)
layer_5 <- (myrbm_4$expected_hidden(layer_4))

for(k in 1:10){
for(i in 1:(nrow(mnist)))
{
    myrbm_5$learn(layer_5[sample(1:nrow(mnist), 5, replace=T), ], 0.1, 0.1, 0.0001)
    if(i %% 10 == 0)
    {
        cat(sprintf("%i is done.\n", i))
    }
}
}
randomat <- myrbm_5$expected_hidden(myrbm_4$expected_hidden(myrbm_3$expected_hidden(myrbm_2$expected_hidden(myrbm$expected_hidden(mnist_bin)))))
reduced_data <- data.frame(digit = labs, DIM_1 = randomat[, 1], DIM_2 = randomat[, 2])

qplot(data=reduced_data,x=DIM_1,y=DIM_2,color=digit)

randomat <-  myrbm_3$expected_hidden(myrbm_2$expected_hidden(myrbm$expected_hidden(mnist_bin)))



randomat <-  myrbm_4$expected_hidden(myrbm_3$expected_hidden(myrbm_2$expected_hidden(myrbm$expected_hidden(mnist_bin))))


randim <- myrbm$expected_visible(myrbm_2$expected_visible(myrbm_3$expected_visible(myrbm_4$expected_visible(randomat))))

for(i in seq(1, 600, 50))
{
    orig <- ggplot(melt((image_gen(i, data = mnist))), aes(x=Var1,y=Var2))
    orig <- orig + geom_tile(aes(fill=value)) + 
        scale_fill_gradient(low = "white", high = "black") + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank())  
    orig <- orig + labs(title = "Original Digit")
    
    learned <- ggplot(melt((image_gen(i, data = randim))), aes(x=Var1,y=Var2))
    learned <- learned + geom_tile(aes(fill=value))+ 
        scale_fill_gradient(low = "white", high = "black") + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank())
    
    learned <- learned + labs(title = "Reconstructed/Learned Digit")
    
    multiplot(orig, learned, cols=2)
}

















# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}