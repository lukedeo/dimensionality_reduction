source("Reduce/R/lin_embed.R")
library("rgl")
library("spatstat")

#dynamic eigenproblem

#online linear embedding....how would it work?

#samuel kou


# for given methods, plot versus K'





dists <- pdist(X, X)

#Generate the Swiss-Roll
N = 1000
r = seq(0, 1, length.out=N)
t = (3*pi/2)*(1+2*r)
x = t*cos(t) + rnorm(N, 0, .5) #add noise
y = t*sin(t) + rnorm(N, 0, .5) #add noise
# z <- c(10 * runif(N/2, 0, 1), 20 * runif(N/2, 0, 1))
z <- 20 * runif(N, 0, 1)

plot3d(x, y, z, col = rainbow(N, start=0, end = .7), pch=19)



colors <- rainbow(n, start=0, end = .7)
plot3d(X[, 1], X[, 2], X[, 3], col = rainbow(N, start=0, end = .7), pch=19)

data = scale(as.matrix(cbind(x, y, z)))

X = data

data <- data.frame(x=x, y=y, z=z)


X <- scale(X)


diffmap_local <- manifold(X, 2, sigma=0.3, t=2, method = "diffusion")
diffmap_semiglobal <- manifold(X, 2, sigma=0.4, t=2, method = "diffusion")
diffmap_global <- manifold(X, 2, sigma=0.9, t=2, method = "diffusion")

diffmap_local_t <- manifold(X, 2, sigma=0.3, t=5, method = "diffusion")
diffmap_semiglobal_t <- manifold(X, 2, sigma=0.4, t=5, method = "diffusion")
diffmap_global_t <- manifold(X, 2, sigma=0.9, t=5, method = "diffusion")

plot(diffmap_local$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(diffmap_semiglobal$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(diffmap_global$Y, pch=19, col=rainbow(N, start=0, end = .7))


lle_local <- manifold(X, d=2, k=5, method="normal")
lle_semiglobal <- manifold(X, d=2, k=8, method="normal")
lle_global <- manifold(X, d=2, k=16, method="normal")

plot(lle_local$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(lle_semiglobal$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(lle_global$Y, pch=19, col=rainbow(N, start=0, end = .7))




new_data <- manifold(X, 2, 6, method="laplacian")
lle <- manifold(X, d=2, k=10, method="normal") # not using formulas
# new_data <- manifold(~x + y + z - 1, data, 2, 8) #working with formulas
plot(diffmap$Y, pch=19, col=rainbow(N, start=0, end = .7))
rd <- LEIGENMAP(X, 2, 10)
plot(laplacian4$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(laplacian6$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(laplacian10$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(rd$Y, pch=19, col=rainbow(N, start=0, end = .7))

plot((Y), pch=19, col=rainbow(N, start=0, end = .7))

laplacian4 <- manifold(X, 2, 4, method="laplacian", heat=0)
laplacian6 <- manifold(X, 2, 6, method="laplacian", heat=0)
laplacian10 <- manifold(X, 2, 10, method="laplacian", heat=4)

kvals <- c(1:8, seq(10, 50, 5))






n = 3
N = 1000


X <-  matrix(0, ncol = 3, nrow = N)
s <- runif(N)

t <- seq(0, 2, length.out=N)



X[, 1] <- -cos(1.5 * pi * t)
X[, 2] <- s
X[t <= 1, 3] <- 2 * (-sin(1.5 * pi * t[t <= 1]))  #+ rnorm(N, 0, .6)[t <= 1]
X[t > 1, 3] <- 2 * (2 + sin(1.5 * pi * t[t > 1])) #+ rnorm(N, 0, .6)[t > 1]


idx <- sort(X[, 3], index.return = TRUE)$ix

X <- X[idx, ]

