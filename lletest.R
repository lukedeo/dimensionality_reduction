source("Reduce/R/lin_embed.R")

#Generate the Swiss-Roll
N = 1000
r = seq(0, 1, length.out=N)
t = (3*pi/2)*(1+2*r)
x = t*cos( t) + rnorm(N, 0, 1) #add noise
y = t*sin(t) + rnorm(N, 0, 1) #add noise
#z <- c(10 * runif(N/2, 0, 1), 20 * runif(N/2, 0, 1))
z <- 20 * runif(N, 0, 1)

plot3d(x, y, z, col = rainbow(N, start=0, end = .7), pch=19)

data = scale(as.matrix(cbind(x, y, z)))

data <- data.frame(x=x, y=y, z=z)

new_data <- lin_embed(data, 2, 8) # not using formulas
new_data <- lin_embed(~x + y + z - 1, data, 2, 8) #working with formulas
rd <- LEIGENMAP(data, 2, 3)
plot(rd$Y, pch=19, col=rainbow(N, start=0, end = .7))

