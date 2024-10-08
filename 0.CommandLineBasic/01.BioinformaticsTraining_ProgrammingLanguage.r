# String
"a"
# Vector
a <- c(1, 2, 3)
a[2]
names(a) <- c("a","b","c")
a["b"]
# Matrix
m <- matrix(seq(6), nrow = 2, ncol = 3)
dim(m)
m[1, 2]

# plot()
cars
plot(cars[,1], cars[,2])
plot(cars)
plot(cars, main = "Speed and Stopping Distances of Cars (1920)", xlab = "Speed(mph)", ylab = "Stopping distance (ft)")
# points()
points(cars[,1], cars[,2], type = "o", pch = 15)
# abline()
abline(v = 15, lty = 2)
m <- lm(dist ~ speed, data=cars)
abline(m)
# barplot()
# boxplot()