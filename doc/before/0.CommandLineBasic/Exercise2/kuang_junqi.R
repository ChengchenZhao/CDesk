#construct matrix
data = matrix(AirPassengers,nrow = 12,ncol = 12,byrow = T)
colnames(data) = c('Jan',"Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
rownames(data) = c('1949','1950','1951','1952','1953','1954','1955','1956','1957','1958','1959','1960')
dim(data)
head(data)

#abstract data
data.sub = t(data[,1:3])
head(data.sub)
barplot(data.sub,xlab = 'year',ylab = 'Passenger Number(K)', beside = T, 
        col = c('red','green','yellow'),legend.text = row.names(data.sub), args.legend = list(x=6))

