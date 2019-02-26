#Making Heatmaps

library(readr)


#Pull out desired gene and sort by length.  Extract Transcript IDs
Data <- read_delim("~/Dropbox/Ellis_RNASEQ/Analysis_Redo/Ctrl_DM1.htseq.edgeR.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Data2 <- read_delim("C:/Users/Hang/Dropbox/Ellis_RNASEQ/Analysis_Redo/DM1_DM2.htseq.edgeR.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Combined = merge(Data, Data2, by.x="X1", by.y="X1", all= TRUE)

#DM1 andDM2
#FormattedCombined = Combined[,c(1,16,5,7:10,28:29)]
#DM1 only
#FormattedCombined = Combined[,c(1,16,5,7:10)]
#DM2 only
FormattedCombined = Combined[,c(1,16,5,7:8,28:29)]

#ListGenes = c("ANXA1","FBLN1","PITPNM1","S100A6","SPOCK2","SNCA","THBS1","THBS4","CXCL12", "CYR61","ITGA1", "ITGA2","EMILIN1","CRIM1","DDR2","EPHA4","NTRK2","PDGFRB")

#CalciumDM1
#ListGenes= c("ANXA1","FBLN1","PITPNM1","S100A6","SLC24A3","SNCA","SPOCK1","SPOCK2","SYT11","THBS1","THBS4")

#AdhesionDM1
#ListGenes= c("EMILIN1","ITGA1","ITGA2")

#IntegrinDM1
#ListGenes = c("CXCL12","CYR61","EMILIN1","SFRP2","THBS1","THBS4")

#EndopeptidaseDM2
#ListGenes = c("ITIH2","SERPINA5","SERPINB9","SERPINF2")

#AdhesionDM2
#ListGenes = c("EMILIN1","ITGA1","ITGA2","ITGA11")

#MapKActivityDM2

#ListGenes = c("CRIM1","DDR2","EPHA4","NTRK2","PDGFRB")
#ParsedGenes = as.data.frame(FormattedCombined[FormattedCombined$gene.x %in% ListGenes,])
#rownames(ParsedGenes) = ParsedGenes$gene.x

#SortedGenes = ParsedGenes[ListGenes,]
library(RColorBrewer)
setwd("~/Desktop/interesting_gene/")


ParsedGenes = as.data.frame(m129WT_LvVsRV.htseq.edgeR.txtigene[m129WT_LvVsRV.htseq.edgeR.txtigene$gene,])
rownames(ParsedGenes) = ParsedGenes$gene
ParsedGenes <- ParsedGenes[ order(row.names(ParsedGenes)), ]#place in alphabetical order
testscale = t(scale(t(ParsedGenes[,7:10]), center = TRUE, scale = TRUE))
View(testscale)
# t means transpose
heatmap.2(as.matrix(testscale), Rowv = FALSE, Colv = FALSE, scale = "none", 
col = my_palette, 
trace = "none", dendrogram='none') #colsep = 1:ncol(testscale)) 
#rowsep = 1:nrow(testscale), sepwidth=c(0.1,0.1), sepcolor='Black')
#Add gradient
install.packages("RColorBrewer", dependencies = TRUE)
library(RColorBrewer)
my_palette <- colorRampPalette(c("red", "white", "green"))(n = 1000)


?heatmap.2
# set a loop for heat maps
setwd( "~/Desktop/interesting_gene/")
library(RColorBrewer)

#Add gradient
install.packages("RColorBrewer", dependencies = TRUE)
my_palette <- colorRampPalette(c("red", "white", "green"))(n = 1000)
FilePath <- "~/Desktop/interesting_gene/"  
Filename<- list.files(pattern="*.txt")
TotalFile = paste(FilePath,Filename,sep="")
i <- 1
for (i in TotalFile) {
  print (i)
}

for (i in TotalFile) {
  data.1<-read.delim(i) #read in file, total file changes
  ParsedGenes = as.data.frame(data.1[data.1$gene,])
  rownames(ParsedGenes) = ParsedGenes$gene
  ParsedGenes <- ParsedGenes[order(row.names(ParsedGenes)), ]
  #View(ParsedGenes)
  if (i == "~/Desktop/interesting_gene/m129WT_LvVsRV.htseq.edgeR.txtigene.txt") {
    testscale = t(scale(t(ParsedGenes[,7:10]), center = TRUE, scale = TRUE))}
  if (i == "~/Desktop/interesting_gene/m129KO_LvVsRV.htseq.edgeR.txtigene.txt"){
    testscale = t(scale(t(ParsedGenes[,7:12]), center = TRUE, scale = TRUE))}
  #testscale = t(scale(t(ParsedGenes[,7:10]), center = TRUE, scale = TRUE))
  #View(testscale)
  # t means transpose
  if (i == "~/Desktop/interesting_gene/mD2WT_LvVsRV.htseq.edgeR.txtigene.txt"){
    testscale = t(scale(t(ParsedGenes[,7:11]), center = TRUE, scale = TRUE))
  }
  if (i == "~/Desktop/interesting_gene/mD2KO_LvVsRV.htseq.edgeR.txtigene.txt"){
    testscale = t(scale(t(ParsedGenes[,7:11]), center = TRUE, scale = TRUE))
  }
  heatmap.2(as.matrix(testscale), Rowv = FALSE, Colv = FALSE, scale = "none", 
            col = my_palette, 
            trace = "none", dendrogram='none') #colsep = 1:ncol(testscale)) 
  #rowsep = 1:nrow(testscale), sepwidth=c(0.1,0.1), sepcolor='Black')
  }


#View(testscale)
?scale

library(dplyr)
library(gplots)
library(plotrix)

phenotype_raw<- read.csv("~/QTL/D2129-2_Final/QTLRel/Original QTLRel Input Files/D2129-2a_pdat.csv", na.strings="-", stringsAsFactors=FALSE)
genotype_raw <- read.csv("~/QTL/D2129-2_Final/QTLRel/Original QTLRel Input Files/D2129-2a_gdat.csv", na.strings="-", stringsAsFactors=TRUE)
gmap = read.csv("~/QTL/D2129-2_Final/QTLRel/Original QTLRel Input Files/D2129-2a_gmap.csv", na.strings="-", stringsAsFactors=FALSE)

#MalePhenos = subset(phenotype_raw, Sex =="M")
#phenotype_raw = MalePhenos

#FemalePhenos = subset(phenotype_raw, Sex =="F")

#phenotype_raw = FemalePhenos


#for (i in 3:length(colnames(phenotype_raw))) { 
#  lo
#  
#  TestDF = TestDF[,-2]
  #rownames(TestDF) = TestDF$Row.names
#  TestDF[is.na(TestDF)] <- 4
  
  
  MarkerList=c()
  for (marker in 1:nrow(TestDF)) {
    if (max(TestDF[marker,]) == 2 && max(TestDF[marker,]) == 2) {
      #NewTestDF = TestDF[-marker, ]
      MarkerList = c(MarkerList, marker)
    }
  }
  
  for (i in c(3,4,5,6,7,10,12,13,14)) { 
    PhenoMasses = phenotype_raw[,c(3,4,5,6,7,10,12,13,14)]
    #10 = QUAD
    #3 = BODY
    PhenoName = colnames(phenotype_raw[i])
    BodyMasses = data.frame(phenotype_raw[,i, drop = FALSE])
    BodyMasses = subset(BodyMasses, !is.na(BodyMasses[,1]))
    gmapchr1 = subset(gmap, gmap$chr == "14")
    #gmapchr1chopped = gmapchr1
    gmapchr1chopped = gmapchr1[28:56,]
    #CHr10
    #gmapchr1chopped = gmapchr1[43:64,]
    TransGmap = t(gmapchr1chopped)
    orderedsnps = TransGmap[1,]
    newgenotype_raw = genotype_raw[,(orderedsnps)]
    MergedBodyMass = merge(BodyMasses, newgenotype_raw, by.x = 0, by.y =0)
    library(dplyr)
    library(gplots)
    Weight = MergedBodyMass[,2]
    TestDF = arrange(MergedBodyMass, desc(Weight))
    rownames(TestDF) = paste(TestDF$Row.names, "-", round(TestDF[,2], digits = 2), sep = "")
    TestDF = TestDF[,-2]
    #rownames(TestDF) = TestDF$Row.names
    TestDF[is.na(TestDF)] <- 4
    TestDF = TestDF[,-1]
    colnames(TestDF) = TransGmap[3,]
    
    #Now need to plot for each phenotype... Chromosome & p values of markers on that chromosome
    somePDFPath <- paste("/Volumes/General_Use/Synpo2L_Project/QTL_Heatmap/", PhenoName, ".pdf", sep = "")
    pdf(file=somePDFPath,height=24,width=24) 
    heatmap.2(as.matrix(TestDF), Rowv = FALSE, Colv = FALSE, scale = "none", col = c("Red", "Yellow", "Pink", "Green", "White"), trace = "none", dendrogram='none', main = PhenoName, ylab = "Weight Ranking (Top = Heaviest)")
    dev.off()
  }
  
  #In Geno Files... 1 = D2, 2 = HET, 3 = 129, 4 = MISSING  ## 32nd snp is main one 28-56 is QTL region 11-33cM, 1-66 is ANOVA region from R Script
  #heatmap.2(as.matrix(TestDF), Rowv = FALSE, Colv = FALSE, scale = "none", col = c("Red", "Yellow", "Pink", "Green", "White"), trace = "none", dendrogram='none', colsep = 1:ncol(TestDF), rowsep = 1:nrow(TestDF), sepwidth=c(0.5,0.5), sepcolor='Black')
  
  #test=summary(aov(RV ~ UNC140105795))
  #pval = test[[1]][, 5][1] 
  #setwd("/Volumes/General_Use/Synpo2L_Project/QTL_Heatmap/")
  #for (i in 2:ncol(TestDF)) {
    #meanaggdata <-aggregate(TestDF$RV, by=list(TestDF[,i]), FUN=mean, na.rm=TRUE)
    #sderroraggdata <-aggregate(TestDF$RV, by=list(TestDF[,i]), FUN=std.error, na.rm=TRUE)
    #FinalAggData = cbind(meanaggdata)
   # FinalAggData = cbind(meanaggdata,sderroraggdata)
  #  write.table(aggdata, paste0(colnames(TestDF)[i],".csv"), sep=",", col.names=FALSE, row.names=FALSE)
 # }
  
  
  #p_value_group = c()
  #MarkerNames = c()
  #whatever = data.frame()
  #for (i in 3:ncol(MergedBodyMass)) {
    #x<-as.factor(MergedBodyMass[,i])
    #fit<- lm(MergedBodyMass$RV~x)
   # p_value<-anova(fit)$"Pr(>F)"[1]
  #  MarkerNames = c(MarkerNames, (colnames(MergedBodyMass)[i]))
   # p_value_group<-c(p_value_group,p_value)
  #}
  #test = cbind(MarkerNames)
  #test = cbind(test,p_value_group )
#  write.table(test, paste0("/Volumes/General_Use/Synpo2L_Project/QTL_Heatmap/",PhenoName,".csv"), sep=",", col.names=FALSE, row.names=FALSE)
   #plot(p_value_group)
  #axis(3, at=1:18, labels=MarkerNames)
  #abline(1/.05)
  
  #x<-as.factor(MergedBodyMass[,i])
  #fit<- lm(as.numeric(MergedBodyMass$RV~x)
           #p_value<- summary(fit)$coefficients[,4][2]
   #        p_value<-anova(fit)$"Pr(>F)"[1]
  ##         p_value_group<-c(p_value_group,p_value)
  setwd("~/Desktop/Full_i_gene/")
  install.packages("RColorBrewer", dependencies = TRUE)
  
  my_palette <- colorRampPalette(c("red", "white", "green"))(n = 1000)
  FilePath <- "~/Desktop/Full_i_gene/"  
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

  View(ParsedGenes)
  for (i in TotalFile) {
    data.1<-read.delim(i) #read in file, total file changes
    ParsedGenes = as.data.frame(data.1[data.1$gene,])
    ParsedGenes<- na.omit(ParsedGenes)
    rownames(ParsedGenes) = ParsedGenes$gene
    ParsedGenes <- ParsedGenes[ order(row.names(ParsedGenes)), ]
    if (i == "~/Desktop/Full_i_gene/m129WT_LvVsRV.htseq.edgeR.txtirvlv.txt"){
      testscale <- t(scale(t(ParsedGenes[,7:10]), center = TRUE, scale = TRUE))
      }
    if (i == "~/Desktop/Full_i_gene/m129KO_LvVsRV.htseq.edgeR.txtirvlv.txt"){
      testscale <- t(scale(t(ParsedGenes[,7:12]), center = TRUE, scale = TRUE))
      }
    #testscale = t(scale(t(ParsedGenes[,7:10]), center = TRUE, scale = TRUE))
    #View(testscale)
    # t means transpose
    if (i == "~/Desktop/Full_i_gene/mD2WT_LvVsRV.htseq.edgeR.txtirvlv.txt"){
      testscale <- t(scale(t(ParsedGenes[,7:11]), center = TRUE, scale = TRUE))
    }
    if (i == "~/Desktop/Full_i_gene/mD2KO_LvVsRV.htseq.edgeR.txtirvlv.txt"){
      testscale <- t(scale(t(ParsedGenes[,7:11]), center = TRUE, scale = TRUE))
    }
    heatmap.2(as.matrix(testscale), Rowv = FALSE, Colv = FALSE, scale = "none", 
              col = my_palette, 
              trace = "none", dendrogram='none') #colsep = 1:ncol(testscale)) 
    #rowsep = 1:nrow(testscale), sepwidth=c(0.1,0.1), sepcolor='Black')
  }
  
  # order interesting genes by fold change

  fc.129wt<- 
  
  for (i in TotalFile) {
    data.1<-read.delim(i) #read in file, total file changes
    ParsedGenes = as.data.frame(data.1[data.1$gene,])
    ParsedGenes<- na.omit(ParsedGenes)
    rownames(ParsedGenes) = ParsedGenes$gene
    ParsedGenes <- ParsedGenes[ order(row.names(ParsedGenes)), ]
    byfc.i <- ParsedGenes[order(logFC),]
  }
 
  
           