###################################################################################################
################################     PARAMETER/FUNCTION DEFINE     ################################
###################################################################################################
# sys_argv <- commandArgs(T)
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")
options(scipen = 200)

red_col <- colorRampPalette(c(cccol[1],cccol[3]), bias=1)(6)[2:5]
blue_col <- colorRampPalette(c(cccol[2],cccol[4]), bias=1)(6)[2:5]
yellow_col <- colorRampPalette(c(cccol[5],cccol[6]), bias=1)(6)[2:5]
grey_col <- colorRampPalette(c(cccol[8],"black"), bias=1)(6)[2:5]

library(beeswarm)

ErrorBar <- function(x, y, upper, lower=upper, length=0.03,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length,lwd=2, ...)
}

###################################################################################################
################################            READ IN DATA           ################################
###################################################################################################
oct4_SCNT <- c(109016,196605,92897,113669,140209,85011,22468,33946,251171,176308,219162)/10000
oct4_14hSCNT <- c(180865,166795,158131,156627,170467,152566,23838,58737,109777,118432,140013,106791,193415,166532,192199,67666,293522,49760,51611,228411,96679,62220,174757,90550)/10000
oct4_Normal <- c(639236,109986,106946,89212,230662,286922,75062,39431,116415,622145,186764,289000,244070,274862,35202,45353,79679,33667,64520,261873,216301,118689,73094,157693,16858,98296,88434,127525,93784,222516,30970,89824,85367,111319,167157)/10000
cdx2_SCNT <- c(252015,993258,426656,524649,611269,227789,280935,274178,605520,124190,327333)/10000
cdx2_14hSCNT <- c(319947,681749,204717,527939,332647,416257,448869,436144,178298,484340,395633,279555,244778,409240,343419,196690,350404,164554,274186,322469,189213,328262,271158,668611)/10000
cdx2_Normal <- c(717941,209591,266518,149481,914040,457941,234388,214891,383621,664795,832397,793498,853611,694905,169373,344884,397351,286845,316298,349628,389963,412059,408629,398419,384233,345334,369336,450899,392102,383805,324392,482801,493210,384703,460093)/10000

cellnum_oct4_SCNT <- c(15,20,14,19,20,7,11,11,18,4,14)
cellnum_oct4_14hSCNT <- c(21,16,12,22,11,12,9,9,7,12,12,10,9,23,25,8,16,12,20,17,8,7,14,10)
cellnum_oct4_Normal <- c(22,26,22,23,19,29,22,27,16,10,17,30,16,26,9,20,23,13,12,7,14,25,16,26,11,18,23,29,23,12,16,17,10,15,18)
cellnum_cdx2_SCNT <- c(23,52,30,31,36,24,25,18,39,12,31)
cellnum_cdx2_14hSCNT <- c(41,46,27,43,31,41,42,41,24,45,37,42,25,42,34,23,36,25,32,42,16,31,29,42)
cellnum_cdx2_Normal <- c(37,35,49,32,63,26,46,54,27,55,70,60,53,55,23,50,58,48,43,57,63,49,30,55,44,25,39,56,60,61,49,62,57,41,64)

cellnumpercentage_oct4_SCNT <- c(0.3947,0.2778,0.3182,0.3800,0.3571,0.2258,0.3056,0.3793,0.3158,0.2500,0.5556,0.3111)*100
cellnumpercentage_oct4_14hSCNT <- c(0.3387,0.2581,0.3077,0.3385,0.2619,0.2264,0.1765,0.1800,0.2258,0.2105,0.2449,0.1923,0.2647,0.3538,0.4237,0.2581,0.3077,0.3243,0.3846,0.2881,0.3333,0.1842,0.3256,0.1176,0.1923)*100
cellnumpercentage_oct4_Normal <- c(0.3729,0.4262,0.3099,0.4182,0.2317,0.5273,0.3235,0.3333,0.3721,0.1538,0.1954,0.3333,0.2319,0.3210,0.2813,0.2857,0.2840,0.2131,0.2182,0.1094,0.1818,0.3378,0.3478,0.3210,0.2000,0.4186,0.3710,0.3412,0.2771,0.1644,0.2462,0.2152,0.1493,0.2679,0.2195)*100
cellnumpercentage_cdx2_SCNT <- c(0.6053,0.7222,0.6818,0.6200,0.6429,0.7742,0.6944,0.6207,0.6842,0.7500,0.4444,0.6889)*100
cellnumpercentage_cdx2_14hSCNT <- c(0.6613,0.7419,0.6923,0.6615,0.7381,0.7736,0.8235,0.8200,0.7742,0.7895,0.7551,0.8077,0.7353,0.6462,0.5763,0.7419,0.6923,0.6757,0.6154,0.7119,0.6667,0.8158,0.6744,0.8824,0.8077)*100
cellnumpercentage_cdx2_Normal <- c(0.6271,0.5738,0.6901,0.5818,0.7683,0.4727,0.6765,0.6667,0.6279,0.8462,0.8046,0.6667,0.7681,0.6790,0.7188,0.7143,0.7160,0.7869,0.7818,0.8906,0.8182,0.6622,0.6522,0.6790,0.8000,0.5814,0.6290,0.6588,0.7229,0.8356,0.7538,0.7848,0.8507,0.7321,0.7805)*100

###################################################################################################
################################                PLOT               ################################
###################################################################################################

pdf("Fig.blastocystIF.20220113.pdf", width=3, height=5)
par(mar=c(9,4,2,2))
plot_list <- list(oct4_SCNT,oct4_14hSCNT,oct4_Normal,cdx2_SCNT,cdx2_14hSCNT,cdx2_Normal)
COL <- c(red_col[1],blue_col[1],grey_col[1],red_col[2],blue_col[2],grey_col[2])
ylim <- range(plot_list,na.rm=T)
dist <- (ylim[2]-ylim[1])/20
beeswarm(plot_list,xaxt="n",col = COL,pch = 16,method = "swarm",corral = "wrap",main="",xlab="",ylab="Mean fluorescence intensity (10^4)",las=2,ylim=c(0,ylim[2]+dist));box(lwd=2)
axis(at=seq(length(plot_list)),side=1,labels=c("oct4_SCNT","oct4_14hSCNT","oct4_Normal","cdx2_SCNT","cdx2_14hSCNT","cdx2_Normal"),las=2)
n <- as.numeric(lapply(plot_list,length))
mean <- as.numeric(lapply(plot_list,mean))
sd <- as.numeric(lapply(plot_list,sd))
points(seq(length(plot_list)),mean,type="p",pch=16,col="black")
ErrorBar(seq(length(plot_list)),mean, qt(0.975, n-1) * sd / sqrt(n),col="black")

plot_list <- list(cellnum_oct4_SCNT,cellnum_oct4_14hSCNT,cellnum_oct4_Normal,cellnum_cdx2_SCNT,cellnum_cdx2_14hSCNT,cellnum_cdx2_Normal)
COL <- c(red_col[1],blue_col[1],grey_col[1],red_col[2],blue_col[2],grey_col[2])
ylim <- range(plot_list,na.rm=T)
dist <- (ylim[2]-ylim[1])/20
beeswarm(plot_list,xaxt="n",col = COL,pch = 16,method = "swarm",corral = "wrap",main="",xlab="",ylab="Cell number",las=2,ylim=c(0,ylim[2]+dist));box(lwd=2)
axis(at=seq(length(plot_list)),side=1,labels=c("cellnum_oct4_SCNT","cellnum_oct4_14hSCNT","cellnum_oct4_Normal","cellnum_cdx2_SCNT","cellnum_cdx2_14hSCNT","cellnum_cdx2_Normal"),las=2)
n <- as.numeric(lapply(plot_list,length))
mean <- as.numeric(lapply(plot_list,mean))
sd <- as.numeric(lapply(plot_list,sd))
points(seq(length(plot_list)),mean,type="p",pch=16,col="black")
ErrorBar(seq(length(plot_list)),mean, qt(0.975, n-1) * sd / sqrt(n),col="black")

plot_list <- list(cellnumpercentage_oct4_SCNT,cellnumpercentage_oct4_14hSCNT,cellnumpercentage_oct4_Normal,cellnumpercentage_cdx2_SCNT,cellnumpercentage_cdx2_14hSCNT,cellnumpercentage_cdx2_Normal)
COL <- c(red_col[1],blue_col[1],grey_col[1],red_col[2],blue_col[2],grey_col[2])
ylim <- range(plot_list,na.rm=T)
dist <- (ylim[2]-ylim[1])/20
beeswarm(plot_list,xaxt="n",col = COL,pch = 16,method = "swarm",corral = "wrap",main="",xlab="",ylab="Cell number percentage",las=2,ylim=c(0,ylim[2]+dist));box(lwd=2)
axis(at=seq(length(plot_list)),side=1,labels=c("cellnumpercentage_oct4_SCNT","cellnumpercentage_oct4_14hSCNT","cellnumpercentage_oct4_Normal","cellnumpercentage_cdx2_SCNT","cellnumpercentage_cdx2_14hSCNT","cellnumpercentage_cdx2_Normal"),las=2)
n <- as.numeric(lapply(plot_list,length))
mean <- as.numeric(lapply(plot_list,mean))
sd <- as.numeric(lapply(plot_list,sd))
points(seq(length(plot_list)),mean,type="p",pch=16,col="black")
ErrorBar(seq(length(plot_list)),mean, qt(0.975, n-1) * sd / sqrt(n),col="black")
dev.off()

p_value <- c()
p_value <- c(p_value,t.test(oct4_SCNT,oct4_14hSCNT)$p.value)
p_value <- c(p_value,t.test(oct4_Normal,oct4_14hSCNT)$p.value)
p_value <- c(p_value,t.test(oct4_SCNT,oct4_Normal)$p.value)
p_value <- c(p_value,t.test(cdx2_SCNT,cdx2_14hSCNT)$p.value)
p_value <- c(p_value,t.test(cdx2_Normal,cdx2_14hSCNT)$p.value)
p_value <- c(p_value,t.test(cdx2_SCNT,cdx2_Normal)$p.value)
p_value <- c(p_value,t.test(cellnum_oct4_SCNT,cellnum_oct4_14hSCNT)$p.value)
p_value <- c(p_value,t.test(cellnum_oct4_Normal,cellnum_oct4_14hSCNT)$p.value)
p_value <- c(p_value,t.test(cellnum_oct4_SCNT,cellnum_oct4_Normal)$p.value)
p_value <- c(p_value,t.test(cellnum_cdx2_SCNT,cellnum_cdx2_14hSCNT)$p.value)
p_value <- c(p_value,t.test(cellnum_cdx2_Normal,cellnum_cdx2_14hSCNT)$p.value)
p_value <- c(p_value,t.test(cellnum_cdx2_SCNT,cellnum_cdx2_Normal)$p.value)
p_value <- c(p_value,t.test(cellnumpercentage_oct4_SCNT,cellnumpercentage_oct4_14hSCNT)$p.value)
p_value <- c(p_value,t.test(cellnumpercentage_oct4_Normal,cellnumpercentage_oct4_14hSCNT)$p.value)
p_value <- c(p_value,t.test(cellnumpercentage_oct4_SCNT,cellnumpercentage_oct4_Normal)$p.value)
p_value <- c(p_value,t.test(cellnumpercentage_cdx2_SCNT,cellnumpercentage_cdx2_14hSCNT)$p.value)
p_value <- c(p_value,t.test(cellnumpercentage_cdx2_Normal,cellnumpercentage_cdx2_14hSCNT)$p.value)
p_value <- c(p_value,t.test(cellnumpercentage_cdx2_SCNT,cellnumpercentage_cdx2_Normal)$p.value)
options(scipen = 1)
print(p_value)

mark <- symnum(p_value, cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","-"))
print(mark)

# p_value <- c()
# p_value <- c(p_value,t.test(oct4_SCNT,oct4_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(oct4_Normal,oct4_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(oct4_SCNT,oct4_Normal,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cdx2_SCNT,cdx2_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cdx2_Normal,cdx2_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cdx2_SCNT,cdx2_Normal,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnum_oct4_SCNT,cellnum_oct4_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnum_oct4_Normal,cellnum_oct4_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnum_oct4_SCNT,cellnum_oct4_Normal,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnum_cdx2_SCNT,cellnum_cdx2_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnum_cdx2_Normal,cellnum_cdx2_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnum_cdx2_SCNT,cellnum_cdx2_Normal,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_oct4_SCNT,cellnumpercentage_oct4_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_oct4_Normal,cellnumpercentage_oct4_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_oct4_SCNT,cellnumpercentage_oct4_Normal,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_cdx2_SCNT,cellnumpercentage_cdx2_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_cdx2_Normal,cellnumpercentage_cdx2_14hSCNT,alternative="less")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_cdx2_SCNT,cellnumpercentage_cdx2_Normal,alternative="less")$p.value)
# options(scipen = 1)
# print(p_value)

# mark <- symnum(p_value, cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","-"))
# print(mark)

# p_value <- c()
# p_value <- c(p_value,t.test(oct4_SCNT,oct4_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(oct4_Normal,oct4_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(oct4_SCNT,oct4_Normal,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cdx2_SCNT,cdx2_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cdx2_Normal,cdx2_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cdx2_SCNT,cdx2_Normal,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnum_oct4_SCNT,cellnum_oct4_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnum_oct4_Normal,cellnum_oct4_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnum_oct4_SCNT,cellnum_oct4_Normal,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnum_cdx2_SCNT,cellnum_cdx2_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnum_cdx2_Normal,cellnum_cdx2_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnum_cdx2_SCNT,cellnum_cdx2_Normal,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_oct4_SCNT,cellnumpercentage_oct4_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_oct4_Normal,cellnumpercentage_oct4_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_oct4_SCNT,cellnumpercentage_oct4_Normal,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_cdx2_SCNT,cellnumpercentage_cdx2_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_cdx2_Normal,cellnumpercentage_cdx2_14hSCNT,alternative="greater")$p.value)
# p_value <- c(p_value,t.test(cellnumpercentage_cdx2_SCNT,cellnumpercentage_cdx2_Normal,alternative="greater")$p.value)
# options(scipen = 1)
# print(p_value)

# mark <- symnum(p_value, cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","-"))
# print(mark)