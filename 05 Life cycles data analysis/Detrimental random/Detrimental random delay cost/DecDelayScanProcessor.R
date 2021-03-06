setwd('/Volumes/OSX/Users/pichugin/Projects/Fragmentation with costs/Local simulations/Classes abundance and groups size/Detrimental random delay cost 2020/')

Write2File=2


# Size of figure
pdfSizeX=5
pdfSizeY=5

## margin parameters
Bottom=3.5
Left=4
Top=0.5
Right=0.5

#palette generated by iwanthue
# CLR_bin <- "#29beb8"
# CLR_equal <- "#d9aa5a"
# CLR_other <- "#67f7ff"
# CLR_unicell <- "#68c4d5"
CLR_bin <- "#60a86a"
CLR_equal <- "#d9aa5a"
CLR_other <- "#67f7ff"
# CLR_unicell <- "#68c4d5"
CLR_unicell <- "#60a86a"


Data=read.csv("Combined_Delay_Scan_Detrimental.txt")
Data=Data[,3:dim(Data)[2]]

Unicell_ID = read.csv(file = "Unicell4dec_ID.dat")
Binaries_ID = read.csv(file = "Binaries4dec_ID.dat")
Equal_ID = read.csv(  file = "Equal4dec_ID.dat")
Other_ID = read.csv(   file = "Other4dec_ID.dat")


ClassPartitions=matrix(0, nrow = 4, ncol = dim(Data)[2])

for(Risk in 1:dim(Data)[2]){
  print(Risk)
  
  Slice <- Data[, Risk]
  
  Binary_num <- 0
  Equal_Num <- 0 
  Other_Num <- 0
  Unicell_Num <- 0
  
  for(Entry in Slice){
    if(Entry %in% Binaries_ID[,1]){
      Binary_num <- Binary_num + 1
    }
    
    if(Entry %in% Equal_ID[,1]){
      Equal_Num <- Equal_Num + 1
    }
    
    if(Entry %in% Other_ID[,1]){
      Other_Num <- Other_Num + 1
    }
    
    if(Entry %in% Unicell_ID[,1]){
      Unicell_Num <- Unicell_Num + 1
    }
  }
  
  Total_Num <- Binary_num + Equal_Num + Other_Num + Unicell_Num
  ClassPartitions[1, Risk] <- Unicell_Num/Total_Num
  ClassPartitions[2, Risk] <- Binary_num/Total_Num
  ClassPartitions[3, Risk] <- Equal_Num/Total_Num
  ClassPartitions[4, Risk] <- Other_Num/Total_Num
  
}


if(Write2File==1){svg(file = "Dec_ClassesVsDelay_v8.svg", height=pdfSizeY, width=pdfSizeX)}
if(Write2File==2){cairo_pdf(file="Dec_ClassesVsDelay_v8.pdf", height=pdfSizeY, width=pdfSizeX)}
# svg(file="Dec_ClassesVsDelay_v5.svg", height=pdfSizeY, width=pdfSizeX)

par(pty="s", mgp=c(2.4,0.6,0.0), xpd=T, mar=c(Bottom,Left, Top, Right), cex.axis=1.4, las = 1)
CLR_list=c(CLR_unicell, CLR_bin, CLR_equal, CLR_other)

MyBP2<-barplot(ClassPartitions[,], col=CLR_list, axes=F, ann=F, xlab="Fragmentation delay", ylab="Fraction of environments", cex.lab=1.6, xaxs="i", yaxs="i", border=NA, space=0)
SQ=seq(0,1, length.out = 6)
SQT=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")
axis(2, at=SQ, labels=SQT)

SQX=seq(0, 20, length.out=5)
PosX=SQX/length(SQX)*dim(ClassPartitions)[2]
axis(1, at=SQX*10, labels=SQX)

# par(pty="s", mgp=c(2.0,0.6,0.0), xpd=T, mar=c(Bottom,Left, Top, Right), cex.axis=1.1, bty="o")
# CLR_list=c(CLR_bin, CLR_equal, CLR_seed, CLR_other)
# 
# # MyBP2<-barplot(ClassPartitions[,], col=CLR_list, axes=F, ann=F, xlab="Fragmentation delay", ylab="Fraction of fitness landscapes", cex.lab=1.6, xaxs="i", yaxs="i", border=NA, space=0)
# Xarray <- seq(0, 20, length.out = dim(Data)[2])
# plot( c(0,20), c(0,1), type = "l", col="white", lwd = 1, ylim=c(0,1), xaxs="i", yaxs="i", ann=F)
# lines(Xarray, ClassPartitions[1,], col=CLR_list[1], lwd=5)
# # lines(Xarray, ClassPartitions[2,], col=CLR_list[2], lwd=5)
# lines(Xarray, ClassPartitions[3,], col=CLR_list[2], lwd=5)
# # lines(Xarray, ClassPartitions[4,], col=CLR_list[4], lwd=5)
# 
# title(xlab="Fragmentation delay", ylab="Fraction of environments", cex.lab=1.6)


if(Write2File>0){dev.off()}



# 
# LC_ID <- read.csv("LC_ID_FFLD_full.dat")
# 
# MeanFragSize <- c()
# StdFragSize <- c()
# MeanOffNum <- c()
# StdOffNum <- c()
# MeanSizeOff <- c()
# StdSizeOff <- c()
# 
# for(Cost in 1:dim(Data)[2]){
#   print(Cost)
#   Slice <- Data[,Cost]
#   FragSize <- c()
#   OffNum <- c()
#   SizeOff <- c()
#   for(Entry in 1:length(Slice)){
#     E_ID <- Slice[Entry]
#     Index <- which(LC_ID$ID == E_ID)
#     FragSize <- c(FragSize, LC_ID$Combined.size[Index] - 1)
#     OffNum <- c(OffNum, LC_ID$Offspring.number[Index])
#     SizeOff <- c(SizeOff, LC_ID$Combined.size[Index] / LC_ID$Offspring.number[Index])
#   }
#   MeanFragSize <- c(MeanFragSize, mean(FragSize))
#   StdFragSize <- c(StdFragSize, sqrt(var(FragSize)))
#   MeanOffNum <- c(MeanOffNum, mean(OffNum))
#   StdOffNum <- c(StdOffNum, sqrt(var(OffNum)))
#   MeanSizeOff <- c(MeanSizeOff, mean(SizeOff))
#   StdSizeOff <- c(StdSizeOff, sqrt(var(SizeOff)))
# }
# 
# 
# 
# 
# CLR_ParentSize <- "#32736d"
# CLR_OffspringNum <- "#f7acba"
# 
# Xseq <- seq(0,20, length.out = length(MeanFragSize))
# SubsetSeq <- seq(1, length(Xseq), 20)
# DX <- 0.4
# 
# if(Write2File==1){svg(file = "StatisticsSizes_Delay_v4.svg", height=pdfSizeY, width=pdfSizeX)}
# if(Write2File==2){cairo_pdf(file="StatisticsSizes_Delay_v4.pdf", height=pdfSizeY, width=pdfSizeX)}
# # svg(file="StatisticsSizes_Delay_v2.svg", height=pdfSizeY, width=pdfSizeX)
# 
# par(pty="s", mgp=c(2.5,0.6,0.0), xpd=T, mar=c(Bottom,Left, Top, Right+2), cex.axis=1.1, tcl = -0.5)
# plot(c(0,20), c(0,10), type="l", col="white", ann=F, xaxt="n", yaxt="n")
# 
# 
# # Ylow <- MeanFragSize-StdFragSize
# # YHigh <- MeanFragSize+StdFragSize
# # Xleft <- Xseq - DX
# # Xright <- Xseq + DX
# # points(Xseq, MeanFragSize, col=CLR_ParentSize, pch=16)
# # for(i in 1:length(SubsetSeq)){
# #   Ind <- SubsetSeq[i]
# #   lines(c(Xseq[Ind], Xseq[Ind]), c(Ylow[Ind], YHigh[Ind]), col=CLR_ParentSize)
# #   lines(c(Xleft[Ind], Xright[Ind]), c(Ylow[Ind], Ylow[Ind]), col=CLR_ParentSize)
# #   lines(c(Xleft[Ind], Xright[Ind]), c(YHigh[Ind], YHigh[Ind]), col=CLR_ParentSize)
# # }
# 
# 
# SubsetSeq <- Xseq
# Ylow <- MeanSizeOff-StdSizeOff
# YHigh <- MeanSizeOff+StdSizeOff
# Xleft <- Xseq - DX
# Xright <- Xseq + DX
# points(Xseq, MeanSizeOff, col=CLR_OffspringNum, pch=16)
# for(i in 1:length(SubsetSeq)){
#   Ind <- SubsetSeq[i]
#   lines(c(Xseq[Ind], Xseq[Ind]), c(Ylow[Ind], YHigh[Ind]), col=CLR_OffspringNum)
#   lines(c(Xleft[Ind], Xright[Ind]), c(Ylow[Ind], Ylow[Ind]), col=CLR_OffspringNum)
#   lines(c(Xleft[Ind], Xright[Ind]), c(YHigh[Ind], YHigh[Ind]), col=CLR_OffspringNum)
# }
# 
# SubsetSeq <- seq(1, length(Xseq), 20)
# Ylow <- MeanFragSize-StdFragSize
# YHigh <- MeanFragSize+StdFragSize
# Xleft <- Xseq - DX
# Xright <- Xseq + DX
# points(Xseq, MeanFragSize, col=CLR_ParentSize, pch=16)
# for(i in 1:length(SubsetSeq)){
#   Ind <- SubsetSeq[i]
#   lines(c(Xseq[Ind], Xseq[Ind]), c(Ylow[Ind], YHigh[Ind]), col=CLR_ParentSize)
#   lines(c(Xleft[Ind], Xright[Ind]), c(Ylow[Ind], Ylow[Ind]), col=CLR_ParentSize)
#   lines(c(Xleft[Ind], Xright[Ind]), c(YHigh[Ind], YHigh[Ind]), col=CLR_ParentSize)
# }
# 
# 
# # 
# # SubsetSeq <- seq(10, length(Xseq), 20)
# # Ylow <- MeanOffNum-StdOffNum
# # YHigh <- MeanOffNum+StdOffNum
# # Xleft <- Xseq - DX
# # Xright <- Xseq + DX
# # points(Xseq, MeanOffNum, col=CLR_OffspringNum, pch=16)
# # for(i in 1:length(SubsetSeq)){
# #   Ind <- SubsetSeq[i]
# #   lines(c(Xseq[Ind], Xseq[Ind]), c(Ylow[Ind], YHigh[Ind]), col=CLR_OffspringNum)
# #   lines(c(Xleft[Ind], Xright[Ind]), c(Ylow[Ind], Ylow[Ind]), col=CLR_OffspringNum)
# #   lines(c(Xleft[Ind], Xright[Ind]), c(YHigh[Ind], YHigh[Ind]), col=CLR_OffspringNum)
# # }
# 
# axis(2, at=c(0,2,4,6,8,10), labels=c(0,2,4,6,8,10))
# par(tcl = 0.5)
# axis(1, at=c(0,5,10,15,20), labels=F)
# 
# # axis(4, at=c(0,5,10,15,20), labels=c(0,5,10,15,20), col.axis=CLR_OffspringNum, col=CLR_OffspringNum)
# # mtext("Offspring number", side=4, line=2, cex=1.6, col=CLR_OffspringNum)
# # mtext("Size at fragmentation", side=2, line=2, cex=1.6, col=CLR_ParentSize)
# mtext("Size", side=2, line=2, cex=1.6)
# # mtext("Fragmentation delay", side=1, line=2, cex=1.6)
# 
# if(Write2File>0){dev.off()}
# # dev.off()

