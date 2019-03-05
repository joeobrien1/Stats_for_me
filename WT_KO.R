# parse out Names, descriptions logFC, logcpm
setwd("~/Desktop/WT_KO/")
FilePath.i <- "~/Desktop/WT_KO/"  
Filename.i<- list.files(pattern="*.txt")
TotalFile.i = paste(FilePath.i,Filename.i,sep="")# paste fucntion combines strings
TotalFile.i

for (i in TotalFile.i){
  print (i)
}

# Pipeline, parse out data w FC, sig, cpm
#sig
for (i in TotalFile.i){
  for (i in TotalFile.i){
    data.1<-read.delim(i) #read in file, total file changes
    sig<-subset(data.1, PValue<.05)
    sig.1<- paste(i, "sig.txt",sep="")
    write.table(sig, sig.1 ,sep="\t",row.names=FALSE)  #write table as an output into directory file was collected from.. MUST CHANGE NAME!!!
  }
}

#FC
setwd("~/Desktop/Sig.wtko/")
FilePath.i <- "~/Desktop/Sig.wtko/"  
Filename.i<- list.files(pattern="*.txt")
TotalFile.i = paste(FilePath.i,Filename.i,sep="")# paste fucntion combines strings
TotalFile.i


for (i in TotalFile.i){
  data.1<-read.delim(i) #read in file, total file changes
  sig.fc<-subset(data.1, data.1$logFC >= 1 | data.1$logFC <= -1)
  sig.fc.i<- paste(i, "fc.txt",sep="")
  write.table(sig.fc, sig.fc.i ,sep="\t",row.names=FALSE)  #write table as an output into directory file was collected from.. MUST CHANGE NAME!!!
}
#cpm
setwd("~/Desktop/fc.wtko/")
FilePath.i <- "~/Desktop/fc.wtko/"  
Filename.i<- list.files(pattern="*.txt")
TotalFile.i = paste(FilePath.i,Filename.i,sep="")# paste fucntion combines strings
TotalFile.i

for (i in TotalFile.i){
  print (i)
}
for (i in TotalFile.i){
  data.1<-read.delim(i) #read in file, total file changes
  sig.fc.cpm<-subset(data.1, logCPM >= 1)
  sig.fc.cpm.i<- paste(i, "cpm.txt",sep="")
  write.table(sig.fc.cpm, sig.fc.cpm.i ,sep="\t",row.names=FALSE)  #write table as an output into directory file was collected from.. MUST CHANGE NAME!!!
}

# Order by FC

setwd("~/Desktop/WT_KO_final/")
FilePath.i <- "~/Desktop/WT_KO_final/"  
Filename.i<- list.files(pattern="*.txt")
TotalFile.i = paste(FilePath.i,Filename.i,sep="")# paste fucntion combines strings
TotalFile.i

for (i in TotalFile.i){
  print (i)
}
#isolate desired columns
for (i in TotalFile.i){
  print(i)
  data.1<-read.delim(i) #read in file, total file changes
  descriptions = as.matrix(data.1$description)
  genenames = as.matrix(data.1$gene)
  logfc = as.matrix(data.1$logFC)
  logcpm = as.matrix(data.1$logCPM)
  gene.descriptions<- cbind(genenames, descriptions)
  gene.descriptions<- cbind(gene.descriptions, logfc) # addition of third column
  gene.descriptions<- cbind(gene.descriptions, logcpm)
  colnames(gene.descriptions)<- c("gene","description", "logFC", "logCPM") # how to name columns 
  gene.descriptions_2<-str_replace(gene.descriptions[,2], "\\[.*\\]", "") #get rid of bracketed portions
  gene.descriptions_3<- cbind(genenames, gene.descriptions_2)
  gene.descriptions_3<- cbind(gene.descriptions_3, logfc)
  gene.descriptions_3<- cbind(gene.descriptions_3, logcpm)
  colnames(gene.descriptions_3)<- c("gene","description", "logFC", "logCPM")
  gname<- paste(i, "igene.des.txt",sep="")
  write.table(gene.descriptions_3, gname ,sep="\t",row.names=FALSE)  #write table as an output
}



setwd("~/Desktop/wtko.des/")
FilePath.i <- "~/Desktop/wtko.des/"  
Filename.i<- list.files(pattern="*.txt")
TotalFile.i = paste(FilePath.i,Filename.i,sep="")# paste fucntion combines strings
TotalFile.i
#order by fc
for (i in TotalFile.i){
  print (i)
}
for (i in TotalFile.i) {
  data.1<-read.delim(i) #read in file, total file changes
  logfc.i<- data.1[order(data.1$logFC, decreasing = TRUE),]
  logfc.name<- paste(i, "i.fc.wt.high.2.txt",sep="")
  write.table(logfc.i, logfc.name ,sep="\t",row.names=FALSE)  #write table as an output into directory file was collected from.. MUST CHANGE NAME!!!
}
#heatmap
install.packages("RColorBrewer", dependencies = TRUE)
setwd("~/Desktop/WT_KO_final/")
my_palette <- colorRampPalette(c("red", "white", "green"))(n = 1000)
FilePath <- "~/Desktop/WT_KO_final/"  
Filename<- list.files(pattern="*.txt")
TotalFile = paste(FilePath,Filename,sep="")
View(TotalFile)
i <- 1
library(readr)
library(dplyr)
library(gplots)
library(plotrix)
library(RColorBrewer)
for (i in TotalFile) {
  print (i)
}

for (i in TotalFile) {
  data.1<-read.delim(i) #read in file, total file changes
  ParsedGenes = as.data.frame(data.1[data.1$gene,])
  ParsedGenes<- na.omit(ParsedGenes)
  rownames(ParsedGenes) = ParsedGenes$gene
  ParsedGenes <- ParsedGenes[ order((ParsedGenes$logFC), decreasing = FALSE), ] # order by fold change
  if (i == "~/Desktop/WT_KO_final/M129WT_129KO.htseq.edgeR.txtsig.txtfc.txtcpm.txt"){
    testscale <- t(scale(t(ParsedGenes[,7:12]), center = TRUE, scale = TRUE))
  }
  
  #testscale = t(scale(t(ParsedGenes[,7:10]), center = TRUE, scale = TRUE))
  #View(testscale)
  # t means transpose
  if (i == "~/Desktop/WT_KO_final/D2WT_D2KO.htseq.edgeR.txtsig.txtfc.txtcpm.txt"){
    testscale <- t(scale(t(ParsedGenes[,7:12]), center = TRUE, scale = TRUE))
  }
  heatmap.2(as.matrix(testscale), Rowv = FALSE, Colv = FALSE, scale = "none", 
            col = my_palette, 
            trace = "none", dendrogram='none') #colsep = 1:ncol(testscale)) 
  #rowsep = 1:nrow(testscale), sepwidth=c(0.1,0.1), sepcolor='Black')
}

#LV
setwd("~/Desktop/WT_KO_LV/")

# parse out Names, descriptions logFC, logcpm
setwd("~/Desktop/WT_KO_LV/")
FilePath.i <- "~/Desktop/WT_KO_LV/"  
Filename.i<- list.files(pattern="*.txt")
TotalFile.i = paste(FilePath.i,Filename.i,sep="")# paste fucntion combines strings
TotalFile.i

for (i in TotalFile.i){
  print (i)
}

# Pipeline, parse out data w FC, sig, cpm
#sig

for (i in TotalFile.i){
    data.1<-read.delim(i) #read in file, total file changes
    sig<-subset(data.1, PValue<.05)
    sig.1<- paste(i, "sig.LV.txt",sep="")
    write.table(sig, sig.1 ,sep="\t",row.names=FALSE)  #write table as an output into directory file was collected from.. MUST CHANGE NAME!!!
  }

#FC
setwd("~/Desktop/WT_KO_LV_sig/")
FilePath.i <- "~/Desktop/WT_KO_LV_sig/"  
Filename.i<- list.files(pattern="*.txt")
TotalFile.i = paste(FilePath.i,Filename.i,sep="")# paste fucntion combines strings
TotalFile.i


for (i in TotalFile.i){
  data.1<-read.delim(i) #read in file, total file changes
  sig.fc<-subset(data.1, data.1$logFC >= 1 | data.1$logFC <= -1)
  sig.fc.i<- paste(i, "fc.lv.txt",sep="")
  write.table(sig.fc, sig.fc.i ,sep="\t",row.names=FALSE)  #write table as an output into directory file was collected from.. MUST CHANGE NAME!!!
}

#CPM

setwd("~/Desktop/WT_KO_LV_FC/")
FilePath.i <- "~/Desktop/WT_KO_LV_FC/"  
Filename.i<- list.files(pattern="*.txt")
TotalFile.i = paste(FilePath.i,Filename.i,sep="")# paste fucntion combines strings
TotalFile.i

for (i in TotalFile.i){
  print (i)
}
for (i in TotalFile.i){
  data.1<-read.delim(i) #read in file, total file changes
  sig.fc.cpm<-subset(data.1, logCPM >= 1)
  sig.fc.cpm.i<- paste(i, "cpm.lv.txt",sep="")
  write.table(sig.fc.cpm, sig.fc.cpm.i ,sep="\t",row.names=FALSE)  #write table as an output into directory file was collected from.. MUST CHANGE NAME!!!
}

# Order by FC

setwd("~/Desktop/WT_KO_final_LV/")
FilePath.i <- "~/Desktop/WT_KO_final_LV/"  
Filename.i<- list.files(pattern="*.txt")
TotalFile.i = paste(FilePath.i,Filename.i,sep="")# paste fucntion combines strings
TotalFile.i

for (i in TotalFile.i){
  print (i)
}
#isolate desired columns
for (i in TotalFile.i){
  print(i)
  data.1<-read.delim(i) #read in file, total file changes
  descriptions = as.matrix(data.1$description)
  genenames = as.matrix(data.1$gene)
  logfc = as.matrix(data.1$logFC)
  logcpm = as.matrix(data.1$logCPM)
  gene.descriptions<- cbind(genenames, descriptions)
  gene.descriptions<- cbind(gene.descriptions, logfc) # addition of third column
  gene.descriptions<- cbind(gene.descriptions, logcpm)
  colnames(gene.descriptions)<- c("gene","description", "logFC", "logCPM") # how to name columns 
  gene.descriptions_2<-str_replace(gene.descriptions[,2], "\\[.*\\]", "") #get rid of bracketed portions
  gene.descriptions_3<- cbind(genenames, gene.descriptions_2)
  gene.descriptions_3<- cbind(gene.descriptions_3, logfc)
  gene.descriptions_3<- cbind(gene.descriptions_3, logcpm)
  colnames(gene.descriptions_3)<- c("gene","description", "logFC", "logCPM")
  gname<- paste(i, "igene.des.lv.txt",sep="")
  write.table(gene.descriptions_3, gname ,sep="\t",row.names=FALSE)  #write table as an output
}



setwd("~/Desktop/LV_igene_des/")
FilePath.i <- "~/Desktop/LV_igene_des/"  
Filename.i<- list.files(pattern="*.txt")
TotalFile.i = paste(FilePath.i,Filename.i,sep="")# paste fucntion combines strings
TotalFile.i
#order by fc
for (i in TotalFile.i){
  print (i)
}
for (i in TotalFile.i) {
  data.1<-read.delim(i) #read in file, total file changes
  logfc.i<- data.1[order(data.1$logFC, decreasing = TRUE),]
  logfc.name<- paste(i, "i.fc.wt.high.3.txt",sep="")
  write.table(logfc.i, logfc.name ,sep="\t",row.names=FALSE)  #write table as an output into directory file was collected from.. MUST CHANGE NAME!!!
}
#heatmap

setwd("~/Desktop/WT_KO_final_LV/")
my_palette <- colorRampPalette(c("red", "white", "green"))(n = 1000)
FilePath <- "~/Desktop/WT_KO_final_LV/"  
Filename<- list.files(pattern="*.txt")
TotalFile = paste(FilePath,Filename,sep="")
View(TotalFile)
i <- 1
library(readr)
library(dplyr)
library(gplots)
library(plotrix)
library(RColorBrewer)
for (i in TotalFile) {
  print (i)
}

for (i in TotalFile) {
  data.1<-read.delim(i) #read in file, total file changes
  ParsedGenes = as.data.frame(data.1[data.1$gene,])
  ParsedGenes<- na.omit(ParsedGenes)
  rownames(ParsedGenes) = ParsedGenes$gene
  ParsedGenes <- ParsedGenes[ order((ParsedGenes$logFC), decreasing = FALSE), ] # order by fold change
  if (i == "~/Desktop/WT_KO_final_LV/M129WT_129KO.htseq.edgeR.txtsig.LV.txtfc.lv.txtcpm.lv.txt"){
    testscale <- t(scale(t(ParsedGenes[,8:12]), center = TRUE, scale = TRUE))
  }
  
  #testscale = t(scale(t(ParsedGenes[,7:10]), center = TRUE, scale = TRUE))
  #View(testscale)
  # t means transpose
  if (i == "~/Desktop/WT_KO_final_LV/D2WT_D2KO.htseq.edgeR.txtsig.LV.txtfc.lv.txtcpm.lv.txt"){
    testscale <- t(scale(t(ParsedGenes[,8:12]), center = TRUE, scale = TRUE))
  }
  heatmap.2(as.matrix(testscale), Rowv = FALSE, Colv = FALSE, scale = "none", 
            col = my_palette, 
            trace = "none", dendrogram='none') #colsep = 1:ncol(testscale)) 
  #rowsep = 1:nrow(testscale), sepwidth=c(0.1,0.1), sepcolor='Black')
}



