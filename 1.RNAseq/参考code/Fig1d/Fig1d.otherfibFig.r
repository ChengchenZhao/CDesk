###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol80 <- c("#CE001380","#16557A80","#C7A60980","#87C23280","#64C0AB80","#A14C9480","#15A08C80","#8B7E7580","#1E7CAF80","#EA425F80","#46489A80","#E5003380","#0F231F80","#1187CD80")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)
library(beeswarm)

HalfErrorBar <- function(x, y, upper, lower=upper, length=0.03,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y, angle=90, code=1, length=length, ...)
}

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
data1 <- rbind(c(0.11,0.23,0.52),c(0.33,0.13,0.56))*100
colnames(data1) <- c("0826","0903","0930")
rownames(data1) <- c("control","14h")

data2 <- rbind(c(0.26,0.25),c(0.47,0.33))*100
colnames(data2) <- c("0113","0531")
rownames(data2) <- c("control","14h")


###################################################################################################
################################                PLOT               ################################
###################################################################################################
pdf("otherfibFig_20211019.1.pdf",width=2.5,height=3.6)
par(mar=c(6,6,2,2))
plot_n <- c(length(data1[1,]),length(data1[2,]))
plot_mean <- c(mean(data1[1,]),mean(data1[2,]))
plot_sd <- c(sd(data1[1,]),sd(data1[2,]))
bp <- barplot(plot_mean,col=cccol80[1],border=cccol[1],ylab="Percentage (%)",xlab="",width=1, space=c(0.8,0.2),xlim=c(0,3.8),names.arg=c("control","14h"),main="BDF1_Tail-Fib_NT",las=2,ylim=c(0,100))
box(lwd=2)
beeswarm(at=bp,list(data1[1,],data1[2,]),xaxt="n",spacing = 0.5,col = cccol80[1],pch = 16,cex=0.5,method = "swarm",corral = "wrap",add=T)
# HalfErrorBar(bp,plot_mean, 1.96*plot_sd/sqrt(plot_n),col=cccol[1],lwd=1)
dev.off()

pdf("otherfibFig_20211019.2.pdf",width=2.5,height=3.6)
par(mar=c(6,6,2,2))
plot_n <- c(length(data2[1,]),length(data2[2,]))
plot_mean <- c(mean(data2[1,]),mean(data2[2,]))
plot_sd <- c(sd(data2[1,]),sd(data2[2,]))
bp <- barplot(plot_mean,col=cccol80[1],border=cccol[1],ylab="Percentage (%)",xlab="",width=1, space=c(0.8,0.2),xlim=c(0,3.8),names.arg=c("control","14h"),main="ICR_mef_NT",las=2,ylim=c(0,100))
box(lwd=2)
beeswarm(at=bp,list(data2[1,],data2[2,]),xaxt="n",spacing = 0.5,col = cccol80[1],pch = 16,cex=0.5,method = "swarm",corral = "wrap",add=T)
# HalfErrorBar(bp,plot_mean, 1.96*plot_sd/sqrt(plot_n),col=cccol[1],lwd=1)
dev.off()

t.test(data1[1,],data1[2,])$p.value
# 0.7743993
t.test(data2[1,],data2[2,])$p.value
# 0.284984

t.test(data1[1,],data1[2,],alternative="less")$p.value
# 0.3871996
t.test(data2[1,],data2[2,],alternative="less")$p.value
# 0.142492

t.test(data1[1,],data1[2,],alternative="greater")$p.value
# 0.6128004
t.test(data2[1,],data2[2,],alternative="greater")$p.value
# 0.857508

# p-value of two tailed t-test : 0.7743993
# p-value of one tailed t-test : 0.6128004 (M14h>Ctrl)
# p-value of one tailed t-test : 0.3871996 (Ctrl<M14h)

# p-value of two tailed t-test : 0.284984
# p-value of one tailed t-test : 0.857508 (M14h>Ctrl)
# p-value of one tailed t-test : 0.142492 (Ctrl<M14h)
