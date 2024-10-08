
data <- matrix(AirPassengers, nrow = 12, ncol = 12, byrow = T)
colnames(data) <-c("Jan", "Feb", "Mar", "Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
rownames(data) <-c("1949", "1950", "1951", "1952", "1953","1954", "1955", "1956", "1957", "1958", "1959", "1960")

data2 <- data[,1:3]

barplot(height = matrix (data2, nrow = 3, ncol = 12, byrow = T),names.arg = c("1949", "1950", "1951", "1952", "1953","1954", "1955", "1956", "1957", "1958", "1959", "1960"),
        beside = TRUE, ylim = c(0, 500), las=2,
        xlab = "Year", ylab = "Passenger Numbers(k)", font.lab=2, col = c("red","blue","orange"),)

legend(x = 1, y = 500, bty = "n",legend = c ("Jan","Feb","Mar"), 
       fill = c("red","blue","orange"))




