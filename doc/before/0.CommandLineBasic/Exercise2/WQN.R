data <- matrix(AirPassengers, nrow=12, ncol=12,byrow=T)
colnames(data) <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
rownames(data) <-c("1949","1950","1951","1952","1953","1954","1955","1956","1957","1958","1959","1960")
data[,c(1,2,3)]
wqn <-data[,c(1,2,3)]
t(wqn)
qn <- t(wqn)
barplot(height=qn,beside=TRUE,xlab="year",ylab="Passenger Numbers(K)",col = c("red","yellow","grey"),las="2")
legend("topleft",c("Jan","Feb","Mar"),col = c("red","yellow","grey"),bty="n",pch=15)

  