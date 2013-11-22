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
N = 400
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


# New LLE -----------------------------------------------------------------

linemb <- local_linear_embedding(X=X, k = 11, d = 2, verbose = TRUE)
plot(linemb$Y, pch=19, col=rainbow(N, start=0, end = .7), main = linemb$description)
plot(decomp$v[,c(198, 199)], pch=19, col=rainbow(N, start=0, end = .7))


# Isomap ------------------------------------------------------------------
iso <- embedding(X=X, k = 6, d = 2, verbose=TRUE, method="isomap", mode = "classical", weighted = FALSE)
iso <- isomap(X=X, k = 10, d = 2, verbose=TRUE, weighted=TRUE)
plot(iso$Y, pch=19, col=rainbow(N, start=0, end = .7), main = iso$description)




iso <- embedding(X=X, k = 6, d = 2, method = "lle")
plot(iso,  col=rainbow(N, start=0, end = .7))












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

lmds_pre <- boxcox(D=Do, A=Inb, d = 2, tau = .7, niter = 20, bfgs = TRUE, verbose = TRUE, sample_rate=1)

lmds <- boxcox(D=Do, A=Inb, X1=lmds_pre$Y, d = 2, tau = 1, sample = 20, niter = 600, cmds_start=FALSE, verbose = TRUE)

lmds_cg <- boxcox(D=Do, A=Inb, d = 2, tau = .5, sample = 20, niter = 600, cmds_start=FALSE, verbose = TRUE)




plot(lmds$Y, pch=19, col=rainbow(N, start=0, end = .7), main = lmds$description)
plot(lmds$best, pch=19, col=rainbow(N, start=0, end = .7), main = lmds$description)


plot(cm_gd$best, pch=19, col=rainbow(N, start=0, end = .7))
Y=lmds$Y
plot3d(Y[, 1], Y[, 2], Y[, 3], col = rainbow(N, start=0, end = .7), pch=19)

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

diffmap_local_t <- manifold(X, 2, sigma=0.3, t=5, method = "diffusion")
diffmap_semiglobal_t <- manifold(X, 2, sigma=0.4, t=5, method = "diffusion")
diffmap_global_t <- manifold(X, 2, sigma=0.9, t=5, method = "diffusion")

diffmap_local_t <- DIFFMAP(X, 2, sigma = 0.3, t=2)
plot(diffmap_local_t$Y, pch=19, col=rainbow(N, start=0, end = .7))

diffmap_local <- diffusion_map(X, 2, t=2, sigma=0.8)
plot(diffmap_local$Y, pch=19, col=rainbow(N, start=0, end = .7))

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

laplacian4 <- laplacian_eigenmap(X=X, 2, 7, heat=.5)
plot(laplacian4$Y, pch=19, col=rainbow(N, start=0, end = .7))

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










