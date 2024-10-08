
data <- matrix(AirPassengers,nrow=12,ncol=12,byrow=T)
colnames(data) <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
rownames(data) <- c("1949","1950","1951","1952","1953","1954","1955","1956","1957","1958","1959","1960") 
data2 <- data[,1:3]
data3 <- t(data2)
barplot(data3)
barplot(data3, beside=T,xlab="year",ylab="Passenger Numbers(k)",col=c("red","green","yellow"),las=2)
legend("topleft",legend=rownames(data3), bty="n",col=c("red","green","yellow"),pch=15)