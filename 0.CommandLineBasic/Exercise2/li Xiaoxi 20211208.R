 data <- matrix(AirPassengers,nrow = 12,ncol = 12,byrow = T)
 colnames(data) <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
 rownames(data) <- c('1949','1950','1951','1952','1953','1954','1955','1956','1957','1958','1959','1960')
 new_data = t(data[,1:3])
 barplot(new_data,beside=TRUE,xlab="year",ylab="Passenger Numbers(k)",las=2,col=c('red','green','yellow'))
 legend(x='topleft',legend=rownames(new_data),fill = c('red','green','yellow'),bty = 'n')
 
 
 