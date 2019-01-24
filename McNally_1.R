rvlv.129s<- read.delim("~/Desktop/Raw_Data/129s.rvlv.txt", stringsAsFactors=FALSE)


rvlv.129s_2<- subset(rvlv.129s, logFC > 1 | logFC < -1 & PValue < .05)


qa.129s<-read.delim("~/Desktop/Raw_Data/129s.qa.txt", stringsAsFactors=FALSE)

qa.129s_2<- subset(qa.129s, logFC > 1 | logFC < -1 & PValue < .05)

View(qa.129s_2)

write.table(qa.129s_2, "~/Desktop/Output folder/qa.129s.txt", append = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)

lvrv.129wt<-read.delim("~/Desktop/Raw_Data/129wt.qa.txt")

lvrv.129wt_2<-subset(lvrv.129wt, logFC > 1 | logFC < -1 & PValue < .05)


write.table(lvrv.129wt_2, "~/Desktop/Output folder/lvrv.129wt.txt", append = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)

qa.129wt <- read.delim("~/Desktop/Raw_Data/129wt.qa.txt")

qa.129wt_2 <-subset(qa.129wt, logFC > 1 | logFC < -1 & PValue < .05)

write.table(qa.129wt_2, "~/Desktop/Output folder/lvrv.129wt.txt", append = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)

D2s.rvlv <- read.delim("~/Desktop/Raw_Data/D2s.rvlv.txt", stringsAsFactors=FALSE)


D2s_2.rvlv<-subset(D2s.rvlv, logFC > 1 | logFC < -1 & PValue < .05)


write.table(D2s_2.rvlv, "~/Desktop/Output folder/lvrv.D2s.txt", append = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)

d2s.aq <- read.delim("~/Desktop/Raw_Data/d2s.aq.txt", stringsAsFactors=FALSE)

d2s.aq_2<-subset(d2s.aq, logFC > 1 | logFC < -1 & PValue < .05)

write.table(d2s.aq_2, "~/Desktop/Output folder/d2s.aq.txt", append = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)



d2wt.rvlv<-read.delim("~/Desktop/Output folder/mD2WT_LvVsRV.htseq.edgeR.txt", stringsAsFactors=FALSE)

d2wt.rvlv_2<-subset(d2wt.rvlv, logFC > 1 | logFC < -1 & PValue < .05)

write.table(d2wt.rvlv_2, "~/Desktop/Output folder/d2wt.rvlv.txt", append = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)


d2wt.qa<-read.delim("~/Desktop/Output folder/mD2WT_QuadVsAbs.htseq.edgeR.txt", stringsAsFactors=FALSE)


d2wt.qa_2<-subset(d2wt.qa, logFC > 1 | logFC < -1 & PValue < .05)

write.table(d2wt.qa_2, "~/Desktop/Output folder/d2wt.qa.txt", append = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)


# For loop for gathering all data at once
setwd("~/Desktop/Raw_Data")
temp <- list.files(pattern="*.txt")
i<-1
#temp = list.files(pattern="*edgeR.txt")

for (i in 1:length(temp)){
  assign(temp[i], read.delim(temp[i]))
  ShortedName = substr(temp[i],1,nchar(temp[i])-4) 
  WorkingName = na.omit(get(temp[i]))
  FilterUP = subset(WorkingName, adj.p <= .05 & logFC >= 1 & logCPM > 2)
  FilterUP = FilterUP[order(FilterUP$adj.p),]
  #FilterUP=FilterUP[,c(1,16,20,7,8,9,10,11,14,15)]
  FilterDOWN = subset(WorkingName, adj.p <= .05 & logFC <= -1 & logCPM > 2)
  FilterDOWN = FilterDOWN[order(FilterDOWN$adj.p),]
  #FilterDOWN=FilterDOWN[,c(1,16,20,7,8,9,10,11,14,15)]
  #UpName = paste(ShortedName, "UPCANDIDATEGENES.txt", sep = "_")
  #DownName = paste(ShortedName, "DOWNCANDIDATEGENES.txt", sep = "_")
  #write.table(FilterUP,UpName,sep="\t",row.names=FALSE)
  #write.table(FilterDOWN,DownName,sep="\t",row.names=FALSE)
  Together = rbind(FilterUP,FilterDOWN)
  #DuplicatedTogether = subset(Together,duplicated(gene) | duplicated(gene, fromLast=TRUE))
  #SortedDuplicatedTogether = DuplicatedTogether[order(DuplicatedTogether$gene),]
  TogetherName = paste(ShortedName, "LOGCPM2.txt", sep = "_")
  write.table(Together,TogetherName,sep="\t",row.names=FALSE)
}


RNAexp <- list.files(path = "~/Desktop/Raw_Data/",pattern="*.txt") # "*" indicates the end of the file name
RNAexp
i<-1
for (i in 1:length(RNAexp)){
  assign(RNAexp[i], read.delim(RNAexp[i]))
  sn_1 <- substr(RNAexp[i],1,nchar(RNAexp[i])-4) 
  wname <- na.omit(get(RNAexp[i]))
  high.fc<- subset(wname, PValue <= .05 & logFC >= 1)
  high.fc <- high.fc[order(high.fc$PValue),]
  tot.name = paste("~/Desktop/RD_2/", sn_1, ".txt", sep = "_")
  write.table(high.fc,tot.name,sep="\t",row.names=FALSE)
}
?get
x<- 1
for (x in 1:length(RNAexp)){
  assign(RNAexp[x], read.delim(RNAexp[x]))
  sn_2 <- substr(RNAexp[x],1,nchar(RNAexp[x])-4) 
  mname <- na.omit(get(RNAexp[x]))
  low.fc<- subset(mname, PValue <= .05 & logFC <= (-1))
  low.fc <- low.fc[order(low.fc$PValue),]
  tot.name.2 = paste("~/Desktop/RD_3/", sn_2, ".txt", sep = "_")
  write.table(low.fc,tot.name.2,sep="\t",row.names=FALSE)
}




for (x in 1:length(RNAexp)){
  assign(RNAexp[x], read.delim(RNAexp[x]))
  sn_2 <- substr(RNAexp[x],1,nchar(RNAexp[x])-4) 
  mname <- na.omit(get(RNAexp[x]))
  low.fc<- subset(mname, PValue <= .05 & logFC <= (-1))
  low.fc <- low.fc[order(low.fc$PValue),]
  tot.name.2 = paste("~/Desktop/RD_3/", sn_2, ".txt", sep = "_")
  write.table(low.fc,tot.name.2,sep="\t",row.names=FALSE)
}






hello <- (1:8)

sname

for (x in hello){
  print(x)
}


listgenes<-c("Actn3","Atp2a1","Atp2b2","Fgf12","Myl4","Myl7","Myl1","Mylpf","Ryr1","Scn10a","Tnni2","Tnnt3","Mybpc1","Drd2","Mstn","Myh1","Myh2","Myh4","Myoc","Stac3","Oprd1","Angpt1","F5","Slit2","Cacna1d","Ryr3","Cacna2d2","Sln","Ctgf")

lvrv.129wt.sub<- subset(m129WT_LvVsRV.htseq.edgeR, gene %in% listgenes )

View(lvrv.129wt.sub)

Mean.left.vent.exp.129wt<-rowMeans(lvrv.129wt.sub[,7:8])

mean.vent.129wt<- cbind(lvrv.129wt.sub,Mean.left.vent.exp.129wt)  #how to add a list to a column


Mean.right.vent.exp.129wt<-rowMeans(lvrv.129wt.sub[,9:10])


mean.vent.129wt<-cbind(mean.vent.129wt,Mean.right.vent.exp.129wt)

View(mean.vent.129wt)




lvrv.129s.sub<- subset(m129KO_LvVsRV.htseq.edgeR, gene %in% listgenes )

View(lvrv.129s.sub)

Mean.left.vent.exp.129s<-rowMeans(lvrv.129s.sub[,7:9])

View(Mean.left.vent.exp)

mean.vent.129s<- cbind(lvrv.129s.sub,Mean.left.vent.exp.129s)  #how to add a list to a column


Mean.right.vent.exp.129s<-rowMeans(lvrv.129s.sub[,10:12])


mean.vent.129s<-cbind(mean.vent.129s,Mean.right.vent.exp.129s)

View(mean.vent.129s)


View(d2wt.rvlv)
lvrv.d2.sub<- subset(mD2WT_LvVsRV.htseq.edgeR, gene %in% listgenes )

View(lvrv.d2.sub)

Mean.left.vent.exp.d2wt<-rowMeans(lvrv.d2.sub[,7:9])



mean.vent.d2wt<- cbind(lvrv.d2.sub,Mean.left.vent.exp.d2wt)  #how to add a list to a column


Mean.right.vent.exp.d2wt<-rowMeans(lvrv.d2.sub[,10:12])


mean.vent.d2wt<-cbind(mean.vent.d2wt,Mean.right.vent.exp.d2wt)

View(mean.vent.d2wt)



lvrv.d2s.sub<- subset(mD2KO_LvVsRV.htseq.edgeR, gene %in% listgenes )

View(lvrv.d2s.sub)

Mean.left.vent.exp.d2s<-rowMeans(lvrv.d2s.sub[,7:8])

plot(mean.vent.129wt)

mean.vent.d2s<- cbind(lvrv.d2s.sub,Mean.left.vent.exp.d2s)  #how to add a list to a column


Mean.right.vent.exp.d2s<-rowMeans(lvrv.d2s.sub[,9:11])


mean.vent.d2s<-cbind(mean.vent.d2s,Mean.right.vent.exp.d2s)

ncol(mean.vent.129wt)

plot(mean.vent.129wt[,-c(1:20)])

boxplot(mean.vent.129wt[,21],mean.vent.129wt[,22])

summary(mean.vent.129wt[,21])
summary(mean.vent.129wt[,22])

mean.vent.129wt_1 <- mean.vent.129wt[-4, ]# delete row from dataframe

boxplot(mean.vent.129wt_1[,21],mean.vent.129wt_1[,22])

ncol(mean.vent.129s)

boxplot(mean.vent.129wt[,22],mean.vent.129s[,24])#boxplot right ventrical

boxplot(mean.vent.129wt[,21],mean.vent.129s[,23])#boxplot left ventrical

boxplot(mean.vent.129s[,23],mean.vent.129s[,24])

summary(mean.vent.129s[,23])

summary(mean.vent.129s[,24])

# Ctgf has the largest variation btw ko and wt 
# in KO Ctgf expression is expressed almost equally in both ventricles at 80-100 units
# in WT Ctgf expression is ~200 in Left Ventricle and ~450 in the right ventricle

boxplot(mean.vent.129s[,23],mean.vent.129s[,24], main = "129 Sgcg Knockout", names=c("Left Ventricle","Right Ventricle"))


boxplot(mean.vent.129wt[,21],mean.vent.129wt[,22], main = "129 Sgcg WildType", names=c("Left Ventricle","Right Ventricle"))

boxplot(mean.vent.d2s[,22],mean.vent.d2s[,23], main = "D2 Sgcg Knockout", names=c("Left Ventricle","Right Ventricle"))
ncol(mean.vent.d2s)

ncol(mean.vent.d2wt)
boxplot(mean.vent.d2wt[,23],mean.vent.d2wt[,24], main = "D2 Sgcg WildType", names=c("Left Ventricle","Right Ventricle"))

sarcomeres<-c("Actn3","Myl1","Tnni2","Tnnt3","Mybpc1","Myh1","Myh2","Myh4","Myoc","Mylpf")

lvrv.129wt.sarc<- subset(m129WT_LvVsRV.htseq.edgeR, gene %in% sarcomeres )

Mean.left.vent.exp.129wt.sarc<-rowMeans(lvrv.129wt.sarc[,7:8])

mean.vent.129wt.sarc<- cbind(lvrv.129wt.sarc,Mean.left.vent.exp.129wt.sarc)  #how to add a list to a column


Mean.right.vent.exp.129wt.sarc<-rowMeans(lvrv.129wt.sarc[,9:10])


mean.vent.129wt.sarc<-cbind(mean.vent.129wt.sarc,Mean.right.vent.exp.129wt.sarc)
View(mean.vent.129wt.sarc)


lvrv.129s.sarc<- subset(m129KO_LvVsRV.htseq.edgeR, gene %in% sarcomeres )



Mean.left.vent.exp.129s.sarc<-rowMeans(lvrv.129s.sarc[,7:9])



mean.vent.129s.sarc<- cbind(lvrv.129s.sarc,Mean.left.vent.exp.129s.sarc)  #how to add a list to a column


Mean.right.vent.exp.129s.sarc<-rowMeans(lvrv.129s.sarc[,10:12])


mean.vent.129s.sarc<-cbind(mean.vent.129s.sarc,Mean.right.vent.exp.129s.sarc)

View(mean.vent.129s.sarc)

boxplot(mean.vent.129s.sarc[,23],mean.vent.129s.sarc[,24], main = "129 Sgcg Knockout Sarcomeres", names=c("Left Ventricle","Right Ventricle"))

boxplot(mean.vent.129wt.sarc[,21],mean.vent.129wt.sarc[,22], main = "129 Sgcg WildType Sarcomeres", names=c("Left Ventricle","Right Ventricle"))

boxplot(mean.vent.129wt.sarc[,22],mean.vent.129s.sarc[,24], main = "129 Sgcg Right Ventricle Sarcomeres", names=c("Left Ventricle","Right Ventricle"))

boxplot(mean.vent.129wt.sarc[,21],mean.vent.129s.sarc[,23], main = "129 Sgcg WildType Sarcomeres", names=c("Left Ventricle","Right Ventricle"))

#All 129 sarcomere samples

boxplot(mean.vent.129wt.sarc[,21],mean.vent.129wt.sarc[,22],mean.vent.129s.sarc[,23],mean.vent.129s.sarc[,24], main = "129 Sarcomeres", names=c("LV WT","RV WT","LV KO","RV KO"), col= c("darkred","blue"))        


#D2 sarcomeres

View(mean.vent.d2s.sarc)lvrv.d2.sarc<- subset(mD2WT_LvVsRV.htseq.edgeR, gene %in% sarcomeres)



Mean.left.vent.exp.d2wt.sarc<-rowMeans(lvrv.d2.sarc[,7:9])



mean.vent.d2wt.sarc<- cbind(lvrv.d2.sarc,Mean.left.vent.exp.d2wt.sarc)  #how to add a list to a column


Mean.right.vent.exp.d2wt.sarc<-rowMeans(lvrv.d2.sarc[,10:12])


mean.vent.d2wt.sarc<-cbind(mean.vent.d2wt.sarc,Mean.right.vent.exp.d2wt.sarc)

View(mean.vent.d2wt.sarc)



lvrv.d2s.sarc<- subset(mD2KO_LvVsRV.htseq.edgeR, gene %in% sarcomeres )


Mean.left.vent.exp.d2s.sarc<-rowMeans(lvrv.d2s.sarc[,7:8])



mean.vent.d2s.sarc<- cbind(lvrv.d2s.sarc,Mean.left.vent.exp.d2s.sarc)  #how to add a list to a column


Mean.right.vent.exp.d2s.sarc<-rowMeans(lvrv.d2s.sarc[,9:11])


mean.vent.d2s.sarc<-cbind(mean.vent.d2s.sarc,Mean.right.vent.exp.d2s.sarc)

ncol(mean.vent.d2wt.sarc)
ncol(mean.vent.d2s.sarc)
View(mean.vent.d2s.sarc)

#boxplot all D2 sarcomeres

boxplot(mean.vent.d2wt.sarc[,23],mean.vent.d2wt.sarc[,24],mean.vent.d2s.sarc[,22],mean.vent.d2s.sarc[,23], main = "D2 Sarcomeres", names=c("LV WT","RV WT","LV KO","RV KO"),col= c("darkred","blue"))



#Ion channels 129

ionchannel<-c("Ryr1", "Scn10a", "Stac3","Cacna1d","Ryr3","Cacna2d2","Sln")# ask Andy whether or not to include Sln


lvrv.129wt.ion<- subset(m129WT_LvVsRV.htseq.edgeR, gene %in% ionchannel )


Mean.left.vent.exp.129wt.ion<-rowMeans(lvrv.129wt.ion[,7:8])

mean.vent.129wt.ion<- cbind(lvrv.129wt.ion,Mean.left.vent.exp.129wt.ion)  #how to add a list to a column


Mean.right.vent.exp.129wt.ion<-rowMeans(lvrv.129wt.ion[,9:10])


mean.vent.129wt.ion<-cbind(mean.vent.129wt.ion,Mean.right.vent.exp.129wt.ion)

View(mean.vent.129wt.ion)


lvrv.129s.ion<- subset(m129KO_LvVsRV.htseq.edgeR, gene %in% ionchannel )



Mean.left.vent.exp.129s.ion<-rowMeans(lvrv.129s.ion[,7:9])



mean.vent.129s.ion<- cbind(lvrv.129s.ion,Mean.left.vent.exp.129s.ion)  #how to add a list to a column


Mean.right.vent.exp.129s.ion<-rowMeans(lvrv.129s.ion[,10:12])


mean.vent.129s.ion<-cbind(mean.vent.129s.ion,Mean.right.vent.exp.129s.ion)

View(mean.vent.129s.ion)
#Ion channel boxplot 129

boxplot(mean.vent.129wt.ion[,21],mean.vent.129wt.ion[,22],mean.vent.129s.ion[,23],mean.vent.129s.ion[,24], main = "129 Ion Channel", names=c("LV WT","RV WT","LV KO","RV KO"),col= c("darkred","blue"))



install.packages("ggplot2")
#Ion channels D2
lvrv.d2.ion<- subset(mD2WT_LvVsRV.htseq.edgeR, gene %in% ionchannel)

Mean.left.vent.exp.d2wt.ion<-rowMeans(lvrv.d2.ion[,7:9])



mean.vent.d2wt.ion<- cbind(lvrv.d2.ion,Mean.left.vent.exp.d2wt.ion)  #how to add a list to a column


Mean.right.vent.exp.d2wt.ion<-rowMeans(lvrv.d2.ion[,10:12])


mean.vent.d2wt.ion<-cbind(mean.vent.d2wt.ion,Mean.right.vent.exp.d2wt.ion)

View(mean.vent.d2wt.ion)



lvrv.d2s.ion<- subset(mD2KO_LvVsRV.htseq.edgeR, gene %in% ionchannel )


Mean.left.vent.exp.d2s.ion<-rowMeans(lvrv.d2s.ion[,7:8])



mean.vent.d2s.ion<- cbind(lvrv.d2s.ion,Mean.left.vent.exp.d2s.ion)  #how to add a list to a column


Mean.right.vent.exp.d2s.ion<-rowMeans(lvrv.d2s.ion[,9:11])


mean.vent.d2s.ion<-cbind(mean.vent.d2s.ion,Mean.right.vent.exp.d2s.ion)

View(mean.vent.d2s.ion)



# Ion channel D2 boxplot
boxplot(mean.vent.d2wt.ion[,23],mean.vent.d2wt.ion[,24],mean.vent.d2s.ion[,22],mean.vent.d2s.ion[,23], main = "D2 Ion Channel", names=c("LV WT","RV WT","LV KO","RV KO"),col= c("darkred","blue"))

#Growth Factors

growthf<-c("Ctgf","Fgf12")


lvrv.129wt.gf<- subset(m129WT_LvVsRV.htseq.edgeR, gene %in% growthf )

Mean.left.vent.exp.129wt.gf<-rowMeans(lvrv.129wt.gf[,7:8])

mean.vent.129wt.gf<- cbind(lvrv.129wt.gf,Mean.left.vent.exp.129wt.gf)  #how to add a list to a column


Mean.right.vent.exp.129wt.gf<-rowMeans(lvrv.129wt.gf[,9:10])


mean.vent.129wt.gf<-cbind(mean.vent.129wt.gf,Mean.right.vent.exp.129wt.gf)

View(mean.vent.129wt.gf)




lvrv.129s.gf<- subset(m129KO_LvVsRV.htseq.edgeR, gene %in% growthf )



Mean.left.vent.exp.129s.gf<-rowMeans(lvrv.129s.gf[,7:9])



mean.vent.129s.gf<- cbind(lvrv.129s.gf,Mean.left.vent.exp.129s.gf)  #how to add a list to a column


Mean.right.vent.exp.129s.gf<-rowMeans(lvrv.129s.gf[,10:12])


mean.vent.129s.gf<-cbind(mean.vent.129s.gf,Mean.right.vent.exp.129s.gf)

View(mean.vent.129s.gf)


#Ctgf 129
ctgf.tot<-c(mean.vent.129wt.gf[1,21],mean.vent.129wt.gf[1,22],mean.vent.129s.gf[1,23],mean.vent.129s.gf[1,24],mean.vent.d2wt.gf[1,23],mean.vent.d2wt.gf[1,24],mean.vent.d2s.gf[1,22],mean.vent.d2s.gf[1,23])
?barplot

barplot(ctgf.tot, names.arg = c("129 LV WT","129 RV WT","129 LV KO","129 RV KO","D2 LV WT","D2 RV WT","D2 LV KO","D2 RV KO"), main = "Ctgf", col= c("darkred","blue"))


#Error bars ctgf: standard error = sd/sqrt(n) 

se.129wt.lv.ctfg<- (sd(lvrv.129wt.gf[1,7:8]))/(sqrt(2))

se.129wt.rv.ctfg<- (sd(lvrv.129wt.gf[1,9:10]))/(sqrt(2))



means.ctgf.tot <- c(201.5874946,456.614560,78.53786,94.12483, 118.4080732,168.756103, 84.7232250,112.334334)
names.ctgf.tot <- c("129 LV WT","129 RV WT","129 LV KO","129 RV KO","D2 LV WT","D2 RV WT","D2 LV KO","D2 RV KO")
standardErrors.ctgf <- c(se.129wt.lv.ctfg,se.129wt.rv.ctfg,se.129s.lv.ctgf,se.129s.rv.ctgf,se.d2wt.lv.ctgf,se.d2wt.rv.ctgf,se.d2s.lv.ctgf,se.d2s.rv.ctgf)
plotTop <- max(means.ctgf.tot+standardErrors.ctgf*2)
ctgf.bar<- barplot(means.ctgf.tot, names.arg=names.ctgf.tot, main = "Ctgf", xlab = "Tissue type", ylab = "Copies per million (cpm)",col=c("darkred","blue"), las=1, ylim=c(0,plotTop))

segments(ctgf.bar, means.ctgf.tot-standardErrors.ctgf*2, ctgf.bar, means.ctgf.tot+standardErrors.ctgf*2, lwd=2)

arrows(ctgf.bar, means.ctgf.tot + standardErrors.ctgf, ctgf.bar, means.ctgf.tot-standardErrors.ctgf, lwd=2, angle=90, code=3)

se.129s.lv.ctgf<- (sd(lvrv.129s.gf[1,7:9]))/(sqrt(3))
se.129s.rv.ctgf<-(sd(lvrv.129s.gf[1,10:12]))/(sqrt(3))

se.d2wt.lv.ctgf<- (sd(lvrv.d2wt.gf[1,7:9]))/(sqrt(3))

se.d2wt.rv.ctgf<- (sd(lvrv.d2wt.gf[1,10:12]))/(sqrt(3))

se.d2s.lv.ctgf<-(sd(lvrv.d2wt.gf[1,7:8]))/(sqrt(2))

se.d2s.rv.ctgf<-(sd(lvrv.d2wt.gf[1,9:11]))/(sqrt(3))




#D2 GF
lvrv.d2wt.gf<- subset(mD2WT_LvVsRV.htseq.edgeR, gene %in% growthf)

View(lvrv.d2wt.gf)
Mean.left.vent.exp.d2wt.gf<-rowMeans(lvrv.d2wt.gf[,7:9])



mean.vent.d2wt.gf<- cbind(lvrv.d2wt.gf,Mean.left.vent.exp.d2wt.gf)  #how to add a list to a column


Mean.right.vent.exp.d2wt.gf<-rowMeans(lvrv.d2wt.gf[,10:12])


mean.vent.d2wt.gf<-cbind(mean.vent.d2wt.gf,Mean.right.vent.exp.d2wt.gf)

View(mean.vent.d2wt.gf)



lvrv.d2s.gf<- subset(mD2KO_LvVsRV.htseq.edgeR, gene %in% growthf )


Mean.left.vent.exp.d2s.gf<-rowMeans(lvrv.d2s.gf[,7:8])



mean.vent.d2s.gf<- cbind(lvrv.d2s.gf,Mean.left.vent.exp.d2s.gf)  #how to add a list to a column


Mean.right.vent.exp.d2s.gf<-rowMeans(lvrv.d2s.gf[,9:11])


mean.vent.d2s.gf<-cbind(mean.vent.d2s.gf,Mean.right.vent.exp.d2s.gf)

View(mean.vent.d2s.gf)


ncol(mean.vent.d2s.gf)

ncol(mean.vent.d2wt.gf)
# Cgtf D2 Bargraph
cgtf.D2<-c(mean.vent.d2wt.gf[1,23],mean.vent.d2wt.gf[1,24],mean.vent.d2s.gf[1,22],mean.vent.d2s.gf[1,23])

barplot(cgtf.D2, names.arg = c("Left Ventricle WT","Right Ventricle WT","Left Ventricle KO","Right Ventricle KO"), main = "D2 Ctgf", col= c("darkred","blue"))

# Total 129 and D2 CGTF expression levels
cgtf.tot<-c(mean.vent.129wt.gf[1,21],mean.vent.129wt.gf[1,22],mean.vent.129s.gf[1,23],mean.vent.129s.gf[1,24],mean.vent.d2wt.gf[1,23],mean.vent.d2wt.gf[1,24],mean.vent.d2s.gf[1,22],mean.vent.d2s.gf[1,23])

barplot(cgtf.tot, names.arg = c("129 LV WT","129 RV WT","129 LV KO","129 RV KO","D2 LV WT","D2 RV WT","D2 LV KO","D2 RV KO"), main = "Ctgf", col= c("darkred","blue"), ylim = c(0,500))


# Total 129 and D2 Fgf12 expression levels

fgf12.tot<-c(mean.vent.129wt.gf[2,21],mean.vent.129wt.gf[2,22], 0, 0,mean.vent.d2wt.gf[2,23],mean.vent.d2wt.gf[2,24],mean.vent.d2s.gf[2,22],mean.vent.d2s.gf[2,23])

barplot(fgf12.tot, names.arg = c("129 LV WT","129 RV WT","129 LV KO", "129 RV KO","D2 LV WT","D2 RV WT","D2 LV KO","D2 RV KO"), main = "Fgf12", col= c("darkred","blue"), ylim = c(0,5))


#fgf12 error bars 

# Standard errors
se.129wt.lv.fgf12<- (sd(lvrv.129wt.gf[2,7:8]))/(sqrt(2))

se.129wt.rv.fgf12<- (sd(lvrv.129wt.gf[2,9:10]))/(sqrt(2))

se.129s.lv.fgf12<- 0
se.129s.rv.fgf12<- 0

se.d2wt.lv.fgf12<- (sd(lvrv.d2wt.gf[2,7:9]))/(sqrt(3))

se.d2wt.rv.fgf12<- (sd(lvrv.d2wt.gf[2,10:12]))/(sqrt(3))

se.d2s.lv.fgf12<-(sd(lvrv.d2wt.gf[2,7:8]))/(sqrt(2))

se.d2s.rv.fgf12<-(sd(lvrv.d2wt.gf[2,9:11]))/(sqrt(3))


means.fgf12.tot <- c(0.3138265,4.229763,0,0,0.5633563,1.726794,0.9123662,1.810492)
names.fgf12.tot <- c("129 LV WT","129 RV WT","129 LV KO","129 RV KO","D2 LV WT","D2 RV WT","D2 LV KO","D2 RV KO")
standardErrors.fgf12 <- c(se.129wt.lv.fgf12,se.129wt.rv.fgf12,se.129s.lv.fgf12,se.129s.rv.fgf12,se.d2wt.lv.fgf12,se.d2wt.rv.fgf12,se.d2s.lv.fgf12,se.d2s.rv.fgf12)
plotTop <- max(means.fgf12.tot+standardErrors.fgf12*2)
fgf12.bar<- barplot(means.fgf12.tot, names.arg=names.fgf12.tot, main = "Fgf12", xlab = "Tissue type", ylab = "Copies per million (cpm)",col=c("darkred","blue"), las=1, ylim=c(0,plotTop))

segments(ctgf.bar, means.ctgf.tot-standardErrors.ctgf*2, ctgf.bar, means.ctgf.tot+standardErrors.ctgf*2, lwd=2)

arrows(fgf12.bar, means.fgf12.tot + standardErrors.fgf12, fgf12.bar, means.fgf12.tot-standardErrors.fgf12, lwd=2, angle=90, code=3)

se.129s.lv.ctgf<- (sd(lvrv.129s.gf[1,7:9]))/(sqrt(3))
se.129s.rv.ctgf<-(sd(lvrv.129s.gf[1,10:12]))/(sqrt(3))

se.d2wt.lv.ctgf<- (sd(lvrv.d2wt.gf[1,7:9]))/(sqrt(3))

se.d2wt.rv.ctgf<- (sd(lvrv.d2wt.gf[1,10:12]))/(sqrt(3))

se.d2s.lv.ctgf<-(sd(lvrv.d2wt.gf[1,7:8]))/(sqrt(2))

se.d2s.rv.ctgf<-(sd(lvrv.d2wt.gf[1,9:11]))/(sqrt(3))


View(mean.vent.129wt.gf)
View(mean.vent.129s.gf)
View(mean.vent.d2wt.gf)
View(mean.vent.d2s.gf)


#Error bars
?arrows()



#Metabolic Genes 129

meta<-c("Atp2a1","Atp2b2")


lvrv.129wt.meta<- subset(m129WT_LvVsRV.htseq.edgeR, gene %in% meta )

Mean.left.vent.exp.129wt.meta<-rowMeans(lvrv.129wt.meta[,7:8])

mean.vent.129wt.meta<- cbind(lvrv.129wt.meta,Mean.left.vent.exp.129wt.meta)  #how to add a list to a column


Mean.right.vent.exp.129wt.meta<-rowMeans(lvrv.129wt.meta[,9:10])


mean.vent.129wt.meta<-cbind(mean.vent.129wt.meta,Mean.right.vent.exp.129wt.meta)

View(mean.vent.129wt.meta)




lvrv.129s.meta<- subset(m129KO_LvVsRV.htseq.edgeR, gene %in% meta )



Mean.left.vent.exp.129s.meta<-rowMeans(lvrv.129s.meta[,7:9])



mean.vent.129s.meta<- cbind(lvrv.129s.meta,Mean.left.vent.exp.129s.meta)  #how to add a list to a column


Mean.right.vent.exp.129s.meta<-rowMeans(lvrv.129s.meta[,10:12])


mean.vent.129s.meta<-cbind(mean.vent.129s.meta,Mean.right.vent.exp.129s.meta)

View(mean.vent.129s.meta)

#D2 Metabolic genes
lvrv.d2wt.meta<- subset(mD2WT_LvVsRV.htseq.edgeR, gene %in% meta)


Mean.left.vent.exp.d2wt.meta<-rowMeans(lvrv.d2wt.meta[,7:9])



mean.vent.d2wt.meta<- cbind(lvrv.d2wt.meta,Mean.left.vent.exp.d2wt.meta)  #how to add a list to a column


Mean.right.vent.exp.d2wt.meta<-rowMeans(lvrv.d2wt.meta[,10:12])


mean.vent.d2wt.meta<-cbind(mean.vent.d2wt.meta,Mean.right.vent.exp.d2wt.meta)

View(mean.vent.d2wt.meta)



lvrv.d2s.meta<- subset(mD2KO_LvVsRV.htseq.edgeR, gene %in% meta )


Mean.left.vent.exp.d2s.meta<-rowMeans(lvrv.d2s.meta[,7:8])



mean.vent.d2s.meta<- cbind(lvrv.d2s.meta,Mean.left.vent.exp.d2s.meta)  #how to add a list to a column


Mean.right.vent.exp.d2s.meta<-rowMeans(lvrv.d2s.meta[,9:11])


mean.vent.d2s.meta<-cbind(mean.vent.d2s.meta,Mean.right.vent.exp.d2s.meta)

View(mean.vent.d2s.meta)

#Total Atp2a1 barplot

atp2a1.tot<-c(mean.vent.129wt.meta[1,21],mean.vent.129wt.meta[1,22],mean.vent.129s.meta[1,23],mean.vent.129s.meta[1,24],mean.vent.d2wt.meta[2,23],mean.vent.d2wt.meta[2,24],mean.vent.d2s.meta[2,22],mean.vent.d2s.meta[2,23])

barplot(atp2a1.tot, names.arg = c("129 LV WT","129 RV WT","129 LV KO","129 RV KO","D2 LV WT","D2 RV WT","D2 LV KO","D2 RV KO"), main = "Atp2a1", col= c("darkred","blue"), ylim = c(0,30))


#D2 Atp2b2 barplot
atp2b2.tot<-c(mean.vent.d2wt.meta[1,23],mean.vent.d2wt.meta[1,24],mean.vent.d2s.meta[1,22],mean.vent.d2s.meta[1,23])

barplot(atp2b2.tot, names.arg = c("D2 LV WT","D2 RV WT","D2 LV KO","D2 RV KO") , main = "D2 Atp2b2", col= c("darkred","blue"), ylim = c(0,16))

#Error bars for metabolism
se.129wt.lv.atpa<- (sd(lvrv.129wt.meta[1,7:8]))/(sqrt(2))

se.129wt.rv.atpa<- (sd(lvrv.129wt.meta[1,9:10]))/(sqrt(2))

se.129s.lv.atpa<- (sd(lvrv.129s.meta[1,7:9]))/(sqrt(3))
se.129s.rv.atpa<- (sd(lvrv.129s.meta[1,10:12]))/(sqrt(3))

se.d2wt.lv.atpa<- (sd(lvrv.d2wt.meta[2,7:9]))/(sqrt(3))

se.d2wt.rv.atpa<- (sd(lvrv.d2wt.meta[2,10:12]))/(sqrt(3))

se.d2s.lv.atpa<-(sd(lvrv.d2wt.meta[2,7:8]))/(sqrt(2))

se.d2s.rv.atpa<-(sd(lvrv.d2s.meta[2,9:11]))/(sqrt(3))


means.atpa.tot <- c(8.646328,18.39475,4.931922,26.69169,7.262388,24.07611,7.365246, 7.006671)
names.atpa.tot <- c("129 LV WT","129 RV WT","129 LV KO","129 RV KO","D2 LV WT","D2 RV WT","D2 LV KO","D2 RV KO")
standardErrors.atpa <- c(se.129wt.lv.atpa,se.129wt.rv.atpa,se.129s.lv.atpa,se.129s.rv.atpa,se.d2wt.lv.atpa,se.d2wt.rv.atpa,se.d2s.lv.atpa,se.d2s.rv.atpa)
plotTop <- max(means.atpa.tot+standardErrors.atpa*2)
atpa.bar<- barplot(means.atpa.tot, names.arg=names.atpa.tot, main = "Atp2a1", xlab = "Tissue type", ylab = "Copies per million (cpm)",col=c("darkred","blue"), las=1, ylim=c(0,plotTop))

segments(ctgf.bar, means.ctgf.tot-standardErrors.ctgf*2, ctgf.bar, means.ctgf.tot+standardErrors.ctgf*2, lwd=2)

arrows(atpa.bar, means.atpa.tot + standardErrors.atpa, atpa.bar, means.atpa.tot-standardErrors.atpa, lwd=2, angle=90, code=3)


# Atp2b2

se.d2wt.lv.atpb<- (sd(lvrv.d2wt.meta[1,7:9]))/(sqrt(3))

se.d2wt.rv.atpb<- (sd(lvrv.d2wt.meta[1,10:12]))/(sqrt(3))

se.d2s.lv.atpb<-(sd(lvrv.d2wt.meta[1,7:8]))/(sqrt(2))

se.d2s.rv.atpb<-(sd(lvrv.d2wt.meta[1,9:11]))/(sqrt(3))

means.atpb.tot <- c(4.835167,14.55923,3.399185,9.265303)
names.atpb.tot <- c("D2 LV WT","D2 RV WT","D2 LV KO","D2 RV KO")
standardErrors.atpb <- c(se.d2wt.lv.atpb,se.d2wt.rv.atpb,se.d2s.lv.atpb,se.d2s.rv.atpb)
plotTop <- max(means.atpb.tot+standardErrors.atpb*2)
atpb.bar<- barplot(means.atpb.tot, names.arg=names.atpb.tot, main = "Atp2b2", xlab = "Tissue type", ylab = "Copies per million (cpm)",col=c("darkred","blue"), las=1, ylim=c(0,plotTop))

arrows(atpb.bar, means.atpb.tot + standardErrors.atpb, atpb.bar, means.atpb.tot-standardErrors.atpb, lwd=2, angle=90, code=3)


#129 inflammation

infl<-"Angpt1"


lvrv.129wt.infl<- subset(m129WT_LvVsRV.htseq.edgeR, gene %in% infl )

Mean.left.vent.exp.129wt.infl<-rowMeans(lvrv.129wt.infl[,7:8])

mean.vent.129wt.infl<- cbind(lvrv.129wt.infl,Mean.left.vent.exp.129wt.infl)  #how to add a list to a column


Mean.right.vent.exp.129wt.infl<-rowMeans(lvrv.129wt.infl[,9:10])


mean.vent.129wt.infl<-cbind(mean.vent.129wt.infl,Mean.right.vent.exp.129wt.infl)

View(mean.vent.129wt.infl)




lvrv.129s.infl<- subset(m129KO_LvVsRV.htseq.edgeR, gene %in% infl )



Mean.left.vent.exp.129s.infl<-rowMeans(lvrv.129s.infl[,7:9])



mean.vent.129s.infl<- cbind(lvrv.129s.infl,Mean.left.vent.exp.129s.infl)  #how to add a list to a column


Mean.right.vent.exp.129s.infl<-rowMeans(lvrv.129s.infl[,10:12])


mean.vent.129s.infl<-cbind(mean.vent.129s.infl,Mean.right.vent.exp.129s.infl)

View(mean.vent.129s.infl)


#D2 inflammation
lvrv.d2wt.infl<- subset(mD2WT_LvVsRV.htseq.edgeR, gene %in% infl)


Mean.left.vent.exp.d2wt.infl<-rowMeans(lvrv.d2wt.infl[,7:9])



mean.vent.d2wt.infl<- cbind(lvrv.d2wt.infl,Mean.left.vent.exp.d2wt.infl)  #how to add a list to a column


Mean.right.vent.exp.d2wt.infl<-rowMeans(lvrv.d2wt.infl[,10:12])


mean.vent.d2wt.infl<-cbind(mean.vent.d2wt.infl,Mean.right.vent.exp.d2wt.infl)

View(mean.vent.d2wt.infl)



lvrv.d2s.infl<- subset(mD2KO_LvVsRV.htseq.edgeR, gene %in% infl )


Mean.left.vent.exp.d2s.infl<-rowMeans(lvrv.d2s.infl[,7:8])



mean.vent.d2s.infl<- cbind(lvrv.d2s.infl,Mean.left.vent.exp.d2s.infl)  #how to add a list to a column


Mean.right.vent.exp.d2s.infl<-rowMeans(lvrv.d2s.infl[,9:11])


mean.vent.d2s.infl<-cbind(mean.vent.d2s.infl,Mean.right.vent.exp.d2s.infl)

View(mean.vent.d2s.infl)


angpt1.tot<-c(mean.vent.129wt.infl[1,21],mean.vent.129wt.infl[1,22],mean.vent.129s.infl[1,23],mean.vent.129s.infl[1,24],mean.vent.d2wt.infl[1,23],mean.vent.d2wt.infl[1,24],mean.vent.d2s.infl[1,22],mean.vent.d2s.infl[1,23])

barplot(angpt1.tot, names.arg = c("129 LV WT","129 RV WT","129 LV KO","129 RV KO","D2 LV WT","D2 RV WT","D2 LV KO","D2 RV KO"), main = "Angpt1", col= c("darkred","blue"), ylim = c(0,60))

?arrows()

#Inflammation Barplot
se.129wt.lv.ang<- (sd(lvrv.129wt.infl[1,7:8]))/(sqrt(2))

se.129wt.rv.ang<- (sd(lvrv.129wt.infl[1,9:10]))/(sqrt(2))

se.129s.lv.ang<- (sd(lvrv.129s.infl[1,7:9]))/(sqrt(3))
se.129s.rv.ang<- (sd(lvrv.129s.infl[1,10:12]))/(sqrt(3))

se.d2wt.lv.ang<- (sd(lvrv.d2wt.infl[1,7:9]))/(sqrt(3))

se.d2wt.rv.ang<- (sd(lvrv.d2wt.infl[1,10:12]))/(sqrt(3))

se.d2s.lv.ang<-(sd(lvrv.d2wt.infl[1,7:8]))/(sqrt(2))

se.d2s.rv.ang<-(sd(lvrv.d2s.infl[1,9:11]))/(sqrt(3))


means.ang.tot <- c(36.96996, 50.01708,40.72475,53.33809,10.79351,36.18386,11.468,36.75122)
names.ang.tot <- c("129 LV WT","129 RV WT","129 LV KO","129 RV KO","D2 LV WT","D2 RV WT","D2 LV KO","D2 RV KO")
standardErrors.ang <- c(se.129wt.lv.ang,se.129wt.rv.ang,se.129s.lv.ang,se.129s.rv.ang,se.d2wt.lv.ang,se.d2wt.rv.ang,se.d2s.lv.ang,se.d2s.rv.ang)

plotTop <- max(means.ang.tot+standardErrors.ang*2)
ang.bar<- barplot(means.ang.tot, names.arg=names.ang.tot, main = "Angpt1", xlab = "Tissue type", ylab = "Copies per million (cpm)",col=c("darkred","blue"), las=1, ylim=c(0,plotTop))

segments(ctgf.bar, means.ctgf.tot-standardErrors.ctgf*2, ctgf.bar, means.ctgf.tot+standardErrors.ctgf*2, lwd=2)

arrows(ang.bar, means.ang.tot + standardErrors.ang, ang.bar, means.ang.tot-standardErrors.ang, lwd=2, angle=90, code=3)




# finding total number of unique genes, building pipeline
total.gene.129ko<-m129KO_LvVsRV.htseq.edgeR$gene

length(total.gene.129ko)

tot.gen.129<- c(total.gene.129ko,total.gene.129wt)

length(tot.gen.129)

length(unique(tot.gen.129))

tot.gen.wt<- c(total.gene.129wt,total.gene.D2wt) # found out D2wt was not factored... lead to far more unique genes then anticipated

length(tot.gen.wt)

length(unique(tot.gen.wt))

total.gene.129wt<-m129WT_LvVsRV.htseq.edgeR$gene


total.gene.D2ko<-mD2KO_LvVsRV.htseq.edgeR$gene

total.gene.D2ko

total.gen.ko<- c(total.gene.D2ko,total.gene.129ko)

length(total.gene.D2ko)

length(total.gene.129wt)

length(total.gen.ko)

length(unique(total.gen.ko))





tot.gen.wt<-c(total.gene.129wt,total.gene.D2wt)

length(tot.gen.wt)

length(unique(tot.gen.wt))

total.gene.D2wt

length(total.gene.D2wt)

total.gene.exp<-c(total.gene.129ko,total.gene.129wt,total.gene.D2ko,total.gene.D2wt)

un.d2<- c(total.gene.D2ko,total.gene.D2wt)

length(unique(un.d2))


length(unique(total.gene.exp))

?unique


# Attain number of genes with significant p value 
sig.129ko<-subset(m129KO_LvVsRV.htseq.edgeR, PValue<.05)
length(sig.129ko$gene)

sig.129wt<-subset(m129WT_LvVsRV.htseq.edgeR, PValue<.05)
length(sig.129wt$gene)
View(sig.129wt)

sig.D2ko<-subset(mD2KO_LvVsRV.htseq.edgeR, PValue<.05)
length(sig.D2ko$gene)
sig.D2wt<-subset(mD2WT_LvVsRV.htseq.edgeR, PValue<.05)
length(sig.D2wt$gene)


sig.129ko.gene<-sig.129ko$gene
View(sig.129ko.gene)
sig.129wt.gene<-sig.129wt$gene
View(sig.129wt.gene)
sig.D2ko.gene<-sig.D2ko$gene
View(sig.D2ko.gene)
sig.D2wt.gene<-sig.D2wt$gene
View(sig.D2wt.gene)


s.d2.129.t<- c(sig.D2ko.gene,sig.D2wt.gene,sig.129ko.gene,sig.129wt.gene)


#write


length(s.d2.129.t)

length(unique(s.d2.129.t))

# Fold change greater than 1

sig.129ko_1<-subset(sig.129ko, sig.129ko$logFC >= 1 | sig.129ko$logFC <= -1)
sig.129wt_1<-subset(sig.129wt, sig.129wt$logFC >= 1 | sig.129wt$logFC <= -1)
sig.D2ko_1<-subset(sig.D2ko, sig.D2ko$logFC >= 1 | sig.D2ko$logFC <= -1)
sig.D2wt_1<-subset(sig.D2wt, sig.D2wt$logFC >= 1 | sig.D2wt$logFC <= -1)



View(sig.129ko_1)
View(sig.129wt_1)
View(sig.D2wt_1)
View(sig.D2ko_1)

sig.genes.fc1<- c(sig.129ko_1,sig.129wt_1,sig.D2ko_1,sig.D2wt_1)

length(unique(sig.genes.fc1$gene))
#148 unique genes


write.table(sig.129ko_2, "~/Desktop/Output folder/sig.129ko_2.txt", append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

sig.129ko_2<-na.omit(sig.129ko_1)

View(sig.129ko_2)

#Heatmap

listgenes<-c("Actn3","Atp2a1","Atp2b2","Fgf12","Myl4","Myl7","Myl1","Mylpf","Ryr1","Scn10a","Tnni2","Tnnt3","Mybpc1","Drd2","Mstn","Myh1","Myh2","Myh4","Myoc","Stac3","Oprd1","Angpt1","F5","Slit2","Cacna1d","Ryr3","Cacna2d2","Sln","Ctgf")

lvrv.129wt.sub<- subset(m129WT_LvVsRV.htseq.edgeR, gene %in% listgenes )

View(lvrv.129wt.sub)

Mean.left.vent.exp.129wt<-rowMeans(lvrv.129wt.sub[,7:8])

mean.vent.129wt<- cbind(lvrv.129wt.sub,Mean.left.vent.exp.129wt)  #how to add a list to a column


Mean.right.vent.exp.129wt<-rowMeans(lvrv.129wt.sub[,9:10])


mean.vent.129wt<-cbind(mean.vent.129wt,Mean.right.vent.exp.129wt)

View(mean.vent.129wt)

ncol(mean.vent.129wt)
heat.129wt<-mean.vent.129wt$gene
heat.129wt<-cbind(heat.129wt, mean.vent.129wt[,21])
heat.129wt<-cbind(heat.129wt,mean.vent.129wt[,22])
write.table(heat.129wt,"~/Desktop/Output folder/heat129wt.txt", append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
