source("../Reduce/R/lin_embed.R")
library("rgl")
library("spatstat")

#dynamic eigenproblem

#online linear embedding....how would it work?

#samuel kou


# for given methods, plot versus K'

vars <- c("significance3d", "nVTX", "nTracksAtVtx", "nSingleTracks", "nTracks", "deltaphi", "deltaeta", "meanTrackRapidity", "meanTrackPtRel", "minTrackRapidity", "minTrackPtRel", "maxTrackRapidity", "maxTrackPtRel", "maxSecondaryVertexRho", "maxSecondaryVertexZ", "subMaxSecondaryVertexRho", "subMaxSecondaryVertexZ", "SVInfoPlus_mass", "SVInfoPlus_energyfrac", "SVInfoPlus_normdist", "n_SVInfoPlus_gt_jet", "n_SVInfoPlus_gt_svx","n_SVInfoPlus_2t")

Do <- pdist(X, X)
dist <- pdist(X, X)


# Generate Datasets -------------------------------------------------------
#Generate the Swiss-Roll
N = 10
N = 150
r = seq(0, 1, length.out=N)
t = (3*pi/2)*(1+2*r)
x = t*cos(t) #+ rnorm(N, 0, .5) #add noise
y = t*sin(t) #+ rnorm(N, 0, .5) #add noise
# z <- c(10 * runif(N/2, 0, 1), 20 * runif(N/2, 0, 1))
z <- 20 * runif(N, 0, 1)
data = scale(as.matrix(cbind(x, y, z)))
X = data
X <- scale(X)
Do <- pdist(X, X)
Daux <- apply(Do,2,sort)[k+1,]
Inb <- ifelse(Do>Daux, 0, 1)
X_init <- X[, c(1, 3)]
ng.i.digits <- ng_image_array_gray('USPS Handwritten Digits',p.digits,16,16,invert = TRUE,img_in_row = FALSE)

N = 1000

X <-  matrix(0, ncol = 3, nrow = N)
s <- runif(N)
t <- seq(0, 2, length.out=N)
X[, 1] <- -cos(1.5 * pi * t)
X[, 2] <- s
X[t <= 1, 3] <- 2 * (-sin(1.5 * pi * t[t <= 1]))  + rnorm(N, 0, .3)[t <= 1]
X[t > 1, 3] <- 2 * (2 + sin(1.5 * pi * t[t > 1])) + rnorm(N, 0, .3)[t > 1]
X <- scale(X)
Do <- pdist(X, X)
Daux <- apply(Do,2,sort)[k+1,]
Inb <- ifelse(Do>Daux, 0, 1)
X_init <- X[, c(1, 3)]

plot3d(X[, 1], X[, 2], X[, 3], col = rainbow(N, start=0, end = .7), pch=19)


data(digits)

X <- swiss_roll(n=N)

data(frey)


digits <- t(digits)

labels <- cbind(rep("0", 1100), rep("1", 1100), rep("2", 1100), rep("3", 1100)
                , rep("4", 1100), rep("5", 1100), rep("6", 1100), rep("7", 1100)
                , rep("8", 1100), rep("9", 1100))
color_vec <- cbind(rep("blue", 1100), rep("red", 1100), rep("green", 1100), rep("black", 1100)
                , rep("purple", 1100), rep("grey", 1100), rep("pink", 1100), rep("yellow", 1100)
                , rep("orange", 1100), rep("brown", 1100))
N <- 600
samp <- sort(sample(x=1:nrow(digits), size=N, replace=F))
X <- (digits[samp, ]) / max(digits[samp, ])
labs <- labels[samp]
col_samp <- color_vec[samp]
Do <- fastPdist(X, X)
Daux <- apply(Do,2,sort)[k+1,]
Inb <- ifelse(Do>Daux, 0, 1)
X_init <- X[, c(1, 3)]

library(RnavGraphImageData)
data(frey)

frey <- t(frey)

frey_scaled <- scale(frey)
N <- 1000
X <- frey_scaled[sort(sample(x=1:nrow(frey), size=N, replace=F)), ]
Do <- fastPdist(X, X)
Daux <- apply(Do,2,sort)[k+1,]
Inb <- ifelse(Do>Daux, 0, 1)
X_init <- X[, c(1, 3)]

Xnew <- swiss_roll(n=N, noisy=F)
X <- swiss_roll(n=N)
X <- s_curve(n=N)
Do <- fastPdist(X, X)
Daux <- apply(Do,2,sort)[k+1,]
Inb <- ifelse(Do>Daux, 0, 1)
X_init <- X[, c(1, 3)]
plot3d(Xnew[, 1], Xnew[, 2], Xnew[, 3], col = rainbow(N, start=0, end = .7), pch=19)


frey_lle <- local_linear_embedding(X, k = 10, d = 2, verbose = TRUE)

plot(frey_lle$Y, pch=19, col=rainbow(N, start=0, end = .7), xlab = "Feature 1", ylab = "Feature 2",main = frey_lle$description)
locle10 <- local_linear_embedding(X, k = 10, d = 2, verbose = TRUE)
locle15 <- local_linear_embedding(X, k = 15, d = 2, verbose = TRUE)
locle20 <- local_linear_embedding(X, k = 20, d = 2, verbose = TRUE)
compare_reductions(X_orig=X, locle10, locle15, locle20)

plot(locle15$Y, pch=19, col=rainbow(N, start=0, end = .7), xlab = "Feature 1", ylab = "Feature 2",main = locle15$description)





Y=frey_lle$Y
Dy <- fastPdist(Y, Y)
plot3d(Y[, 1], Y[, 2], Y[, 3], col = rainbow(N, start=0, end = .7), pch=19, main = frey_lle$description)

# New LLE -----------------------------------------------------------------

linemb <- local_linear_embedding(X = X, k = 20, d = 3, verbose = TRUE)
plot(linemb$Y, pch=19, col=rainbow(N, start=0, end = .7), main = linemb$description)
plot(decomp$v[,c(198, 199)], pch=19, col=rainbow(N, start=0, end = .7))

Y=linemb$Y
plot3d(X[, 1], X[, 2], X[, 3], col = rainbow(N, start=0, end = .7), pch=19)



# Isomap ------------------------------------------------------------------
iso <- embedding(X=X, k = 10, d = 2, verbose=TRUE, method="isomap", mode = "classical", weighted = FALSE)

iso10 <- isomap(X=X, k = 10,  d = 2, verbose=TRUE, weighted=TRUE)
iso15 <- isomap(X=X, k = 15,  d = 2, verbose=TRUE, weighted=TRUE)
iso20 <- isomap(X=X, k = 20,  d = 2, verbose=TRUE, weighted=TRUE)
compare_reductions(X_orig=X, iso10, iso15, iso20)
plot(iso15$Y, col=rainbow(N, start=0, end = .7), main = iso15$description, pch = 19)
Y=iso$Y
plot3d(Y[, 1], Y[, 2], Y[, 3], col = rainbow(N, start=0, end = .7), pch=19, main = frey_lle$description)




iso <- embedding(X=X, k = 6, d = 2, method = "lle")
plot(iso,  col=rainbow(N, start=0, end = .7))


A <- fastPdist(X, X)
A <- weighted_neighbor_graph(A, k)
G <- graph.adjacency(A, mode = "undirected", weighted=TRUE)
A <- shortest.paths(G, mode="all")
if(sum(!is.finite(A)) > 0)
{
    A[!is.finite(A)] <- 1.1* max(A[is.finite(A)])
}


Do <- fastPdist(X, X)

k = 10
Daux <- apply(Do,2,sort)[k+1,]
Inb <- ifelse(Do>Daux, 0, 1)
lmds <- boxcox(D=Do, A=Inb, d = 2, sample_rate = 10, niter = 500, cmds_start=T, verbose = T,scale_out=F)
plot(lmds$best, col=rainbow(N, start=0, end = .7), main = lmds$description, pch = 19)



library(kernlab)

X <- s_curve(N, noisy=F)

data <- data.frame(out=lmds$Y, X)

svm_coord1 <- ksvm(out.1~.-out.2, data = data, type = "eps-svr")
svm_coord2 <- ksvm(out.2~.-out.1, data = data, type = "eps-svr")
plot(predict(svm_coord1, data), predict(svm_coord2, data), col = rainbow(N, start=0, end = .7), pch=19, main = "Embedding of LMDS, k = 10 training data using kSVM")


Ntest <- 1000
X_new <- s_curve(Ntest)
X_new[, 1] <- X_new[, 1] + rnorm(n=Ntest, 0, 0.1)
data_new <- data.frame(out.1=rep(0, Ntest), out.2=rep(0, Ntest), X_new)
plot(predict(svm_coord1, data_new), predict(svm_coord2, data_new), col = rainbow(Ntest, start=0, end = .7), pch=19, , main = "Embedding of noisy new data \nusing k = 10 LMDS-estimated kSVM")

plot3d(X[, 1], X[, 2], X[, 3], col = rainbow(Ntest, start=0, end = .7), pch=19)

predemb <- list()
predemb$description <- "Embedding generated by estimated kSVM\non out-of-sample data"
predemb$Y <- cbind(predict(svm_coord1, data_new), predict(svm_coord2, data_new))

Dnew <- fastPdist(X_new, X_new)

k = 10
Da <- apply(Dnew,2,sort)[k+1,]
Adj <- ifelse(Dnew>Da, 0, 1)

lmds_spec <- boxcox(D=Dnew, A=Adj, d = 2, tau = 1, lambda=1, sample_rate = 10, niter = 200, cmds_start=F, verbose = T,scale_out=F)





lmds_frey_clust <- boxcox(D=Do, A=Inb, d = 2, tau = 1, lambda=1, sample_rate = 10, niter = 150, cmds_start=F, verbose = TRUE, bfgs=F, scale_out=F)
lmds_frey <- boxcox(D=Do, A=Inb, d = 2, tau = 1, lambda=1, cmds_start=F, sample_rate = 1, niter = 10, verbose = TRUE, bfgs=T, X1=lmds_frey$Y, scale_out=F)

k = 5
Do <- fastPdist(X, X)
Daux <- apply(Do,2,sort)[k+1,]
Inb <- ifelse(Do>Daux, 0, 1)
lmds_k5 <- boxcox(D=Do, A=Inb, d = 2, tau = 1, lambda=1, cmds_start=F, sample_rate = 10, niter = 300, verbose = TRUE, bfgs=F, scale_out=F)

k = 10
Do <- fastPdist(X, X)
Daux <- apply(Do,2,sort)[k+1,]
Inb <- ifelse(Do>Daux, 0, 1)
lmds_k10 <- boxcox(D=Do, A=Inb, d = 2, tau = 1, lambda=1, cmds_start=F, sample_rate = 10, niter = 300, verbose = TRUE, bfgs=F, scale_out=F)

k = 15
Do <- fastPdist(X, X)
Daux <- apply(Do,2,sort)[k+1,]
Inb <- ifelse(Do>Daux, 0, 1)
lmds_k15 <- boxcox(D=Do, A=Inb, d = 2, tau = 1, lambda=1, cmds_start=F, sample_rate = 10, niter = 300, verbose = TRUE, bfgs=F, scale_out=F)
lmds = lmds_k5 
plot(lmds$Y, col=rainbow(N, start=0, end = .7), main = lmds$description, pch = 19)
lmds = lmds_k10
plot(lmds$Y, col=rainbow(N, start=0, end = .7), main = lmds$description, pch = 19)
lmds = lmds_k15
plot(lmds$Y, col=rainbow(N, start=0, end = .7), main = lmds$description, pch = 19)

compare_reductions(X_orig=X, lmds_k5, lmds_k10, lmds_k15)

Y = lmds_frey_geodesic$Y
Y=lmds_frey_clust$best

plot3d(Y[, 1], Y[, 2], Y[, 3], col = rainbow(N, start=0, end = .7), pch=19, main = lmds_frey$description)



rand_prj <- random_projection(X=X, d=2, type=2)
plot(rand_prj$Y, col=rainbow(N, start=0, end = .7), main = rand_prj$description, pch = 19)
Y = rand_prj$Y
plot3d(Y[, 1], Y[, 2], Y[, 3], col = rainbow(N, start=0, end = .7), pch=19, main = rand_prj$description)


cms <- cmds(D=Do, d=2, verbose=T)
plot(cms$Y, col=rainbow(N, start=0, end = .7), main = cms$description, pch = 19)




plot(princomp(iso$Y)$scores[, 1:2], pch=19, col=rainbow(N, start=0, end = .7), main =iso$description)
plot(iso$Y, pch=19, col=rainbow(N, start=0, end = .7), main =iso$description)
plot(A, pch=19, col=rainbow(N, start=0, end = .7), main =A$description)


# Doing LMDS --------------------------------------------------------------


cm1 <- boxcox(dist=Do, Adj=Inb, X1 = X_init, random.start=0, d = 2, niter=500)
cm2 <- boxcox(dist=Do, Adj=Inb, X1 = cm1$X, random.start=0, d = 2, niter=500)
cm <- BOXCOX(D=Do, A=Inb, X1=X_init, cmds_start=0, random_start=0, d = 2)



plot(cm$embedding, pch=19, col=rainbow(N, start=0, end = .7))
plot(cm$best, pch=19, col=rainbow(N, start=0, end = .7))


plot(cm2$X, pch=19, col=rainbow(N, start=0, end = .7))


cm_bfgs <- BOXCOX(D=Do, A=Inb, X1=X_init, cmds_start=1, random_start=0, d = 2, sample_rate=2, niter=30, bfgs=1, tau=1)

cm_bfgs <- BOXCOX(D=Do, A=Inb, X1=X_init, cmds_start=1, random_start=0, d = 2, sample_rate=2, niter=40, bfgs=1, tau=1)

cm_gd <- BOXCOX(D=Do, A=Inb, X1=cm_bfgs$embedding, cmds_start=0, random_start=0, d = 2, sample_rate=50, niter=700, bfgs=0, tau=1)

cm_gd <- BOXCOX(D=Do, A=Inb, X1=X_init, cmds_start=1, random_start=0, d = 2, sample_rate=50, niter=1000, bfgs=0, tau=1)
cm_gd_bt <- BOXCOX(D=Do, A=Inb, X1=X_init, cmds_start=1, random_start=0, d = 2, sample_rate=2, niter=1000, bfgs=0, tau=1, adaptive=0)


plot(X_init, pch=19, col=rainbow(N, start=0, end = .7))


cm_gd <- boxcox(D=Do, A=Inb, d=3, tau = .2, niter=120)


#adaptive lambda

lmds_pre <- boxcox(D=Do, A=Inb, d = 2, tau = .6, niter = 5, bfgs = TRUE, verbose = TRUE, sample_rate=1)

lmds <- boxcox(D=Do, A=Inb, d = 2, tau = 1, sample = 20, niter = 500, cmds_start=TRUE, verbose = TRUE)

t1 <- boxcox(D=Do, A=Inb, d = 3, tau = 6, sample_rate = 10, niter = 300, cmds_start=T, verbose = TRUE)
t1 <- boxcox(D=Do, A=Inb, d = 3, tau = 1, sample_rate = 10, niter = 300, X1=t1$best,cmds_start=F, verbose = TRUE)
t1 <- boxcox(D=Do, A=Inb, d = 3, tau = .7, sample_rate = 10, niter = 300,X1=t1$Y,cmds_start=F, verbose = TRUE)
t1 <- boxcox(D=Do, A=Inb, d = 3, tau = .3, lambda = 2,cmds_start=F, sample_rate= 10, niter = 300, X1=t1$Y, verbose = TRUE)






f
t.5 <- boxcox(D=Do, A=Inb, d = 3, tau = .5, sample = 10, niter = 300, cmds_start=T, verbose = TRUE)
t.1 <- boxcox(D=Do, A=Inb, d = 3, tau = .1, sample = 10, niter = 300, cmds_start=T, verbose = TRUE)



lmds = t1
lmds = t.5
lmds = t.1

plot(lmds$Y, pch=19, col=rainbow(N, start=0, end = .7), main = lmds$description)
plot(lmds$best, pch=19, col=rainbow(N, start=0, end = .7), main = lmds$description)


Y=lmds$best
plot3d(Y[, 1], Y[, 2], Y[, 3], col = rainbow(N, start=0, end = .7), pch=19, size = 5)

plot(cm_gd_bt$embedding, pch=19, col=rainbow(N, start=0, end = .7))
plot(cm_gd_bt$best, pch=19, col=rainbow(N, start=0, end = .7))

plot(cm_bfgs$embedding, pch=19, col=rainbow(N, start=0, end = .7))
plot(cm_bfgs$best, pch=19, col=rainbow(N, start=0, end = .7))



cm_bfgs <- BOXCOX(D=Do, A=Inb, X1=cm_gd$embedding, cmds_start=0, random_start=0, d = 2, sample_rate=1, niter=30, bfgs=1, tau=1)

plot(emb$Y, pch=19, col=rainbow(N, start=0, end = .7), main = emb$description)


plot3d(x, y, z, col = rainbow(N, start=0, end = .7), pch=19)



colors <- rainbow(n, start=0, end = .7)
plot3d(X[, 1], X[, 2], X[, 3], col = rainbow(N, start=0, end = .7), pch=19)

data = scale(as.matrix(cbind(x, y, z)))

X = data

data <- data.frame(x=x, y=y, z=z)


net <- nnet(X, X, size=2)



X <- scale(X)

lle_local_quark <- manifold(X_p, d=2, k=10, method="normal")

cmds <- list(X=X, Y=cmdscale(dists, 2), description="CMDS")



# Diffusion Map -----------------------------------------------------------


diffmap_local <- manifold(X, 2, sigma=0.3, t=2, method = "diffusion")
diffmap_semiglobal <- manifold(X, 2, sigma=0.4, t=2, method = "diffusion")
diffmap_global <- manifold(X, 2, sigma=0.9, t=2, method = "diffusion")

diffmap_local_t <- manifold(X, 2, sigma=0.1, t=5, method = "diffusion")
diffmap_semiglobal_t <- manifold(X, 2, sigma=0.4, t=5, method = "diffusion")
diffmap_global_t <- manifold(X, 2, sigma=0.9, t=5, method = "diffusion")

diffmap_local_t <- DIFFMAP(X, 2, sigma = 0.3, t=2)
plot(diffmap_local_t$Y, pch=19, col=rainbow(N, start=0, end = .7))

diffmap_local <- diffusion_map(X, 3, t=2, sigma=-1, verbose=T)
plot(diffmap_local$Y, pch=19, col=rainbow(N, start=0, end = .7))
Y=diffmap_local$Y
plot3d(Y[, 1], Y[, 2], Y[, 3], col = rainbow(N, start=0, end = .7), pch=19, main = diffmap_local$description)


plot(diffmap_local_t$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(diffmap_semiglobal$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(diffmap_global_t$Y, pch=19, col=rainbow(N, start=0, end = .7))


# LLe ---------------------------------------------------------------------



lle_local <- manifold(X, d=2, k=5, method="normal")
lle_semiglobal <- manifold(X, d=2, k=10, method="normal")
lle_semiglobal2 <- manifold(X, d=2, k=11, method="normal")
lle_semiglobal_sym <- manifold(X, d=2, k=8, method="normal")
lle_global <- manifold(X, d=2, k=30, method="normal")
proj <- manifold(X, d=2, method="projection")

plot(lle_local$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(lle_semiglobal$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(lle_global$Y, pch=19, col=rainbow(N, start=0, end = .7))



# Laplacian ---------------------------------------------------------------


laplacian5 <- manifold(X, 2, 4, method="laplacian", heat=0)
plot(laplacian5$Y, pch=19, col=rainbow(N, start=0, end = .7))

laplacian5 <- laplacian_eigenmap(X=X,d=3, k=11, heat=100.0, verbose=T)
plot(laplacian5$Y, pch=19, col=rainbow(N, start=0, end = .7))

Y=laplacian5$Y
plot3d(Y[, 1], Y[, 2], Y[, 3], col = rainbow(N, start=0, end = .7), pch=19, main = laplacian5$description)


old <- laplacian4$Y
new <- laplacian4$Y
plot(double(dec$vectors[, 298:297]), pch=19, col=rainbow(N, start=0, end = .7))



dec4 <- eigen(laplacian4$L)
dec5 <- eigen(laplacian5$L)



new_data <- manifold(X, 2, 6, method="laplacian")
lle <- manifold(X, d=2, k=10, method="normal") # not using formulas
# new_data <- manifold(~x + y + z - 1, data, 2, 8) #working with formulas
plot(diffmap$Y, pch=19, col=rainbow(N, start=0, end = .7))
rd <- LEIGENMAP(X, 2, 10)
plot(laplacian4$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(laplacian6$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(laplacian10$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(rd$Y, pch=19, col=rainbow(N, start=0, end = .7))

plot((cm$vectors[, c(1,2)]), pch=19, col=rainbow(N, start=0, end = .7))

laplacian4 <- manifold(X, 2, 4, method="laplacian", heat=0)
laplacian6 <- manifold(X, 2, 6, method="laplacian", heat=2)
laplacian10 <- manifold(X, 2, 10, method="laplacian", heat=4)

kvals <- c(1:8, seq(10, 50, 5))


kprin <- kpca(X)
Y <- kprin@pcv[, 1:2]
ksv <- list(X=X, Y=Y, description="Kernel PCA")



# scurve=X
# sroll=X

X=scale(sroll)
X=scale(scurve)

n = 3
N = 1000


X <-  matrix(0, ncol = 3, nrow = N)
s <- runif(N)

t <- seq(0, 2, length.out=N)



X[, 1] <- -cos(1.5 * pi * t)
X[, 2] <- s
X[t <= 1, 3] <- 2 * (-sin(1.5 * pi * t[t <= 1]))  #+ rnorm(N, 0, .3)[t <= 1]
X[t > 1, 3] <- 2 * (2 + sin(1.5 * pi * t[t > 1]))# + rnorm(N, 0, .3)[t > 1]


idx <- sort(X[, 3], index.return = TRUE)$ix

X <- X[idx, ]









pca <- princomp(X)

dists <- pdist(X, X)
Y <- cmdscale(dists, k=2)



Y <- manifold(X, d=2, k=11, method="normal")$Y



plot((Y), pch=19, col=rainbow(N, start=0, end = .7))



Y <- pca$scores[, 1:2]







plot_2D_projection <- function(manifold, dataset_name){
    label <- paste(dataset_name, ", ",manifold$description, sep = "")
    pdf(paste0(label, ".pdf"))
    plot(manifold$Y, pch=19, col=rainbow(N, start=0, end = .7), xlab = "Component 1", ylab = "Component 2", main = label)
    dev.off()
}

for(set in c("Swiss Roll", "S-Curve")) {
    if(set == "Swiss Roll")
    {
        X <- sroll
    }
    if(set == "S-Curve")
    {
        X <- scurve
    }
    manifold_t <- manifold(X, d=2, k=5, method="normal")
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, d=2, k=8, method="normal")
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, d=2, k=11, method="normal")
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, d=2, k=30, method="normal")
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, 2, 4, method="laplacian", heat=0)
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, 2, 6, method="laplacian", heat=2)
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, 2, 10, method="laplacian", heat=4)
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, 2, sigma=0.3, t=2, method = "diffusion")
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, 2, sigma=0.4, t=2, method = "diffusion")
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, 2, sigma=0.9, t=2, method = "diffusion")
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, 2, sigma=0.3, t=5, method = "diffusion")
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, 2, sigma=0.4, t=5, method = "diffusion")
    plot_2D_projection(manifold_t, set)
    manifold_t <- manifold(X, 2, sigma=0.9, t=5, method = "diffusion")
    plot_2D_projection(manifold_t, set)
}

X <- s_curve(n=1000) 
X <- swiss_roll(n=1000)

lle1 <- local_linear_embedding(X=X, k=20, d=2, verbose=T)
lle2 <- local_linear_embedding(X=X, k=15, d=2, verbose=T)
lle3 <- local_linear_embedding(X=X, k=10, d=2, verbose=T)
lle4 <- local_linear_embedding(X=X, k=6, d=2, verbose=T)

iso1 <- isomap(X=X, k=20, d = 2, verbose=T)
iso2 <- isomap(X=X, k=15, d = 2, verbose=T)
iso3 <- isomap(X=X, k=10, d = 2, verbose=T)
iso4 <- isomap(X=X, k=5, d = 2, verbose=T)

candidates <- list(lle1, lle2, lle3, lle4, iso1, iso2, iso3, iso4)
X_orig <- X
D_orig <- fastPdist(X_orig, X_orig)
kvals <- c(1:10, seq(15, 30, 5), 40, 50, 100)

LC <- data.frame(k=rep(0, length(candidates) * length(kvals)), crit=rep(0, length(candidates) * length(kvals)))


len_k <- length(kvals)
mult = 0
names <- rep("", length(candidates) * length(kvals))
for(manifold in candidates){
    Dnew <- manifold$Y
    Dnew <- fastPdist(Dnew, Dnew)
    idx <- seq(1 + mult * len_k, by=1, to= (mult + 1) * len_k)
    LC$crit[idx] <- knn_criterion(D_orig, Dnew, kvals)
    LC$k[idx] <- kvals
    names[idx] <- (manifold$description)
    mult <- mult + 1
}
idx <- 1
minval <- min(LC)
minval <- ifelse(minval < 0, - 1.14 * abs(minval) , .84 * minval)
maxval <- 1.14 * max(LC)

LC$method <- names
LC$method <- factor(names)

library(scales)
gg <- qplot(x=k,y=crit,color=method, data = LC)+geom_line()

gg <- gg +scale_x_continuous(trans = 'log10',
                             breaks = trans_breaks('log10', function(x) 10^x),
                             labels = trans_format('log10', math_format(10^.x)))


gg <- gg+ xlab(label="Neighborhood size, k")
gg <- gg + ylab("Local MetaCriterion")
gg <- gg+labs(title = "Local MetaCriterion for Considered Methods", cex = 1.3)
gg + theme(plot.title = element_text(size = rel(1.5))) + theme(text=element_text(family="Trebuchet MS")) 
gg + scale_colour_tableau()

gg + theme_solarized(light = FALSE) + scale_colour_solarized("red")





