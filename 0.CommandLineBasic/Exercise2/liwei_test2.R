data <- matrix(AirPassengers, nrow = 12, ncol = 12, byrow = T)
colnames(data)<-  c("Jan","Feb","Mar","Apr","May","Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
rownames(data) <- c("1949","1950","1951","1952","1953","1954","1955","1956","1957","1958","1959","1960") 
new_data <-t(data[,c(1:3)])
barplot(new_data,
        beside =T,las=0,xlab = "year",ylab = "Passenger Numbers(K)",
        col=c('red','yellow','green'),legend=rownames(new_data),
                                      args.legend = list(x = "topleft"))
       

