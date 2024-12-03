
data<-matrix(AirPassengers,nrow=12,ncol=12,byrow=T)

colnames(data)<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

rownames(data)<-c("1949","1950","1951","1952","1953","1954","1955","1956","1957","1958","1959","1960")

data_df <- data[,1:3]

barplot(t(data_df),beside = TRUE, las=2,col = c("red","blue","yellow"),xlab="year", ylab="Passenger Numbers(k)")

legend(x = 1, y = 400,legend = c("Jan","Feb","Mar"),bty = "n",fill = c("red","blue","yellow"))
