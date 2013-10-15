source("Reduce/R/lin_embed.R")
library("rgl")
library("spatstat")

#dynamic eigenproblem

#online linear embedding....how would it work?

#samuel kou







#Generate the Swiss-Roll
N = 1000
r = seq(0, 1, length.out=N)
t = (3*pi/2)*(1+2*r)
x = t*cos(t) + rnorm(N, 0, 1) #add noise
y = t*sin(t) + rnorm(N, 0, 1) #add noise
# z <- c(10 * runif(N/2, 0, 1), 20 * runif(N/2, 0, 1))
z <- 20 * runif(N, 0, 1)

plot3d(x, y, z, col = rainbow(N, start=0, end = .7), pch=19)



colors <- rainbow(n, start=0, end = .7)
plot3d(X[, 1], X[, 2], X[, 3], col = rainbow(N, start=0, end = .7), pch=19)

data = scale(as.matrix(cbind(x, y, z)))

X = data

data <- data.frame(x=x, y=y, z=z)


X <- scale(X)
new_data <- manifold(X, 2, 10, method="laplacian") # not using formulas
new_data <- manifold(X, 2, 10, method="normal") # not using formulas
# new_data <- manifold(~x + y + z - 1, data, 2, 8) #working with formulas
plot(new_data$Y, pch=19, col=rainbow(N, start=0, end = .7))
rd <- LEIGENMAP(X, 2, 10)
plot(reduction$Y, pch=19, col=rainbow(N, start=0, end = .7))
plot(rd$Y, pch=19, col=rainbow(N, start=0, end = .7))

plot(t(Y), pch=19, col=rainbow(N, start=0, end = .7))






n = 3
N = 10


X <-  matrix(0, ncol = 3, nrow = N)
s <- runif(N)

t <- seq(0, 2, length.out=N)



X[, 1] <- -cos(1.5 * pi * t)
X[, 2] <- s
X[t <= 1, 3] <- 2 * (-sin(1.5 * pi * t[t <= 1])) # + rnorm(N, 0, .2)[t <= 1])
X[t > 1, 3] <- 2 * (2 + sin(1.5 * pi * t[t > 1])) # + rnorm(N, 0, .2)[t > 1])


idx <- sort(X[, 3], index.return = TRUE)$ix

X <- X[idx, ]

