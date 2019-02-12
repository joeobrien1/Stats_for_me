assign(RNAexp[i], read.delim(RNAexp[i]))
WorkingName = na.omit(get(RNAexp[i]))
rna.gene = subset(WorkingName, listgenes)
sn_1 <- substr(RNAexp[i],1,nchar(RNAexp[i])-4) 
wname <- na.omit(get(RNAexp[i]))
tot.name = paste("~/Desktop/LVRV_Data/", sn_1, ".txt", sep = "_")


# set working directory
setwd("~/Desktop/LVRV_Data/")
FilePath <- "~/Desktop/LVRV_Data/"  
Filename<- list.files(pattern="*.txt")
TotalFile = paste(FilePath,Filename,sep="")# paste fucntion combines strings
TotalFile
Filename<- list.files(pattern="*.txt")
View(TotalFile)
subset(TotalFile, listgenes)
data.1<-read.delim(TotalFile) #read in file 
listgenes<-c("Actn3","Atp2a1","Atp2b2","Fgf12","Myl4","Myl7","Myl1","Mylpf","Ryr1","Scn10a","Tnni2","Tnnt3","Mybpc1","Drd2","Mstn","Myh1","Myh2","Myh4","Myoc","Stac3","Oprd1","Angpt1","F5","Slit2","Cacna1d","Ryr3","Cacna2d2","Sln","Ctgf")
TotalFile
for (i in TotalFile){
  print (i)
}

#How to loop a subset of your data into an output table format 
for (i in TotalFile){
  print(i)
  data.1<-read.delim(i) #read in file, total file changes
  i.gene<- subset(data.1, gene %in% listgenes) #subset to selected genes
  #if (i == "~/Desktop/LVRV_Data/m129KO_LvVsRV.htseq.edgeR.txt") {4
   # FilteredData = i.gene[,c(3:4,5:7)]
  #}
  #if i == "~/Desktop/LVRV_Data/m129WT_LvVsRV.htseq.edgeR.txt":
   # FilteredData = i.gene[,c(3:4,5:7)]
                  # if statements conditional, will only act if top statement is true
  gname<- paste(i, "igene.txt",sep="")
  #write.table(FilteredData, gname ,sep="\t",row.names=FALSE)#write table as an output
}

# how to take individual columns from data frame and construct a new dataframe
lvrv.129wt.sub<- subset(m129WT_LvVsRV.htseq.edgeR.txtigene, gene %in% listgenes )
View(lvrv.129wt.sub)
descriptions = as.matrix(lvrv.129wt.sub$description)
genenames = as.matrix(lvrv.129wt.sub$gene)
gene.descriptions<- cbind(genenames, descriptions)
colnames(gene.descriptions)<- c("gene","description") # how to name columns 
View(gene.descriptions)
gene.descriptions[,2]
library(stringr)
gene.descriptions_2<-str_replace(gene.descriptions[,2], "\\[.*\\]", "") #get rid of bracketed portions
gene.descriptions_3<- cbind(genenames, gene.descriptions_2)
colnames(gene.descriptions_3)<- c("gene","description")
View(gene.descriptions_3)

# Create and output tables
ed.table<- "edtable.txt"
write.table(gene.descriptions_3, ed.table,sep = "\t",row.names = FALSE)

# for loop, every interesting gene
setwd("~/Desktop/lvrv_output/")
FilePath.i <- "~/Desktop/lvrv_output/"  
Filename.i<- list.files(pattern="*.txt")
TotalFile.i = paste(FilePath.i,Filename.i,sep="")# paste fucntion combines strings
TotalFile.i
View(data.1)

# how to loop every file you are interested in 
for (i in TotalFile.i){
  print(i)
  data.1<-read.delim(i) #read in file, total file changes
  descriptions = as.matrix(data.1$description)
  genenames = as.matrix(data.1$gene)
  logfc = as.matrix(data.1$logFC)
  gene.descriptions<- cbind(genenames, descriptions)
  gene.descriptions<- cbind(gene.descriptions, logfc) # addition of third column
  colnames(gene.descriptions)<- c("gene","description", "logFC") # how to name columns 
  library(stringr)
  gene.descriptions_2<-str_replace(gene.descriptions[,2], "\\[.*\\]", "") #get rid of bracketed portions
  gene.descriptions_3<- cbind(genenames, gene.descriptions_2)
  gene.descriptions_3<- cbind(gene.descriptions_3, logfc)
  colnames(gene.descriptions_3)<- c("gene","description", "logFC")
  gname<- paste(i, "igenesub2.txt",sep="")
  write.table(gene.descriptions_3, gname ,sep="\t",row.names=FALSE)  #write table as an output
}


# compiling all unique datasets
total.data<- rbind(m129KO_LvVsRV.htseq.edgeR.txtigene.txtigenesub2,m129WT_LvVsRV.htseq.edgeR.txtigene.txtigenesub2, mD2KO_LvVsRV.htseq.edgeR.txtigene.txtigenesub2,mD2WT_LvVsRV.htseq.edgeR.txtigene.txtigenesub2)
tdg<-total.data$gene
View(unique(tdg))



