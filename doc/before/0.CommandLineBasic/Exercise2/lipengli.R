
R version 4.1.2 (2021-11-01) -- "Bird Hippie"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)

R是自由软件，不带任何担保。
在某些条件下你可以将其自由散布。
用'license()'或'licence()'来看散布的详细条件。

R是个合作计划，有许多人为之做出了贡献.
用'contributors()'来看合作者的详细情况
用'citation()'会告诉你如何在出版物中正确地引用R或R程序包。

用'demo()'来看一些示范程序，用'help()'来阅读在线帮助文件，或
用'help.start()'通过HTML浏览器来看帮助文件。
用'q()'退出R.

> ？barplot
错误: unexpected input在"？"里
> ?barplot
starting httpd help server ... done
> ?legend
> ?col
> data<-matrix(AirPassengers,nrow=12,ncol=12,byrow=T)#解释构建12行12列的矩阵
> colnames(data)<-c("Jan","Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
>  rownames(data) <- c("1949", "1950", "1951", "1952", "1953", "1954", "1955", "1956", "1957", "1958", "1959", "1960")
>  barplot(t(data[,1:3]),beside=T,col=c("#CE0013","#16557A","#C7A609"),las=2,xlab="year",ylab="Passenger Numbers (k)")
>  legend("topleft",bty="n",c("Jan", "Feb", "Mar"),col=c("#CE0013","#16557A","#C7A609"),pch=15)
> save.image("C:\\Users\\ABC\\Desktop\\lipengli")
> 
