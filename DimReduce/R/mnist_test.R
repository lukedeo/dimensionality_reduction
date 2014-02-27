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
myrbm <- RBM(28^2, 400)

for(i in 1:(nrow(mnist) - 3))
{
    myrbm$learn(ifelse(mnist[i:(1+3), ] > 0.2, 1, 0))
}


reconstructions <- myrbm$reconstruct(mnist, 2)

for(i in 1:50)
{
    orig <- ggplot(melt((image_gen(i, data = mnist))), aes(x=Var1,y=Var2))
    orig <- orig + geom_tile(aes(fill=value)) 
                 + scale_fill_gradient(low = "white", high = "black") 
                 + theme(axis.line=element_blank(),
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
    learned <- learned + geom_tile(aes(fill=value))
                       + scale_fill_gradient(low = "white", high = "black") 
                       + theme(axis.line=element_blank(),
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