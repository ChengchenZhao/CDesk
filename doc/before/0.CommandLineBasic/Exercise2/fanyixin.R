data<-matrix(AirPassengers,nrow=12,ncol=12,byrow=T)
colnames(data)<-c('Jan','Feb',"Mar",'Apr',"May",'Jun','Jul','Aug','Sep','Otc',"Nov",'Dec')
rownames(data)<-c('1949','1950','1951','1952','1953','1954','1955','1956','1957','1958','1959','1960')
dim(data)

data.a=t(data)
data.b=data.a[1:3,]
barplot(data.b,xlab='year',ylab='Passenger Numbers(k)',beside=T,col=c('red','blue','yellow'),las=2)
legend(x='topleft',legend=rownames(data.b),fill=c('red','blue','yellow'),bty='n')
