data <- matrix(AirPassengers, nrow = 12, ncol = 12, byrow = T)
colnames(data) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
rownames(data) <- c("1949", "1950", "1951", "1952", "1953", "1954", "1955", "1956", "1957", "1958", "1959", "1960")
"q" <- data[,1:3]
"tq" <- t(q)
barplot(tq, beside=TRUE, col=c("cyan", "blue", "navy"),ylab="Air Passengers(k)", xlab="year", las=2)
legend(x=1, y=400, legend = c("Jan", "Feb", "Mar"), fill=c("cyan","blue","navy"), bty="n")