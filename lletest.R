source("Reduce/R/lin_embed.R")

#Generate the Swiss-Roll
N = 1000
r = seq(0, 1, length.out=N)
t = (3*pi/2)*(1+2*r)
x = t*cos( t^(1/5)) + rnorm(N, 0, .2) #add noise
y = t*sin(t) + rnorm(N, 0, .2) #add noise
z <- c(10 * runif(N/2, 0, 1), 20 * runif(N/2, 0, 1))
#z <- 20 * runif(N, 0, 1)

plot3d(x, y, z, col = rainbow(N, start=0, end = .7), pch=19)

data <- lle_scurve_data

plot3d(X[, 1], X[, 2], X[, 3], col = rainbow(N, start=0, end = .7), pch=19)

data = scale(as.matrix(cbind(x, y, z)))

X = data

data <- data.frame(x=x, y=y, z=z)

new_data <- manifold(X, 2, 15, method="laplacian") # not using formulas
new_data <- manifold(X, 2, 6, method="normal") # not using formulas
new_data <- manifold(~x + y + z - 1, data, 2, 8) #working with formulas
plot(new_data$Y, pch=19, col=rainbow(N, start=0, end = .7))
rd <- LEIGENMAP(X, 2, 10)
plot(rd$Y, pch=19, col=rainbow(N, start=0, end = .7))






n = 3
N = 1000


X <-  matrix(0, ncol = 3, nrow = N)
angle <- pi * (1.5 * runif(N/2) - 1)
height <- 5 * runif(N)
sd <- 0.03
X[, 1] <- c(cos(angle), -cos(angle)) + rnorm(1, 0, sd)
X[, 2] <- height + rnorm(1, 0, sd)
X[, 3] <- c(sin(angle), 2 - sin(angle)) + rnorm(1, 0, sd)


idx <- sort(X[, 3], index.return = TRUE)$ix

X <- X[idx, ]

