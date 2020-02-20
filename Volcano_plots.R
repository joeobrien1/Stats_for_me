wt129lvrv<- m129WT_LvVsRV.htseq.edgeR
colnames(wt129lvrv)
dim(wt129lvrv)
wt129lvrv <- wt129lvrv[c("gene", "PValue", "logFC")]
head(wt129lvrv)
wt129lvrv["group"]<-"Not Significant"
wt129lvrv[which(wt129lvrv["PValue"] < 0.05 & abs(wt129lvrv["logFC"]) < 1 ),"group"]<- "Significant"
wt129lvrv[which(wt129lvrv["PValue"] > 0.05 & abs(wt129lvrv["logFC"]) > 1 ),"group"] <- "Fold Change"
wt129lvrv[which(wt129lvrv["PValue"] < 0.05 & abs(wt129lvrv["logFC"]) > 1 ),"group"] <- "Significant and Fold Change"#sig p value and fold change
top_peaks <- wt129lvrv[with(wt129lvrv, order(logFC, PValue)),][1:5,]
top_peaks <- rbind(top_peaks, wt129lvrv[with(wt129lvrv, order(-logFC, PValue)),][1:5,])
# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly
a <- list()
for (i in seq_len(nrow(top_peaks))) {
  m <- top_peaks[i, ]
  a[[i]] <- list(
    x = m[["logFC"]],
    y = -log10(m[["PValue"]]),
    text = m[["gene"]],
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 0.5,
    ax = 20,
    ay = -40
  )
}
library("lazyeval")
library(ggplot2)
install.packages("dplyr") 
install.packages("plotly")
install.packages("magrittr")
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
library(plotly)

p <- plot_ly(data = wt129lvrv, x = wt129lvrv$logFC, y = -log10(wt129lvrv$PValue), mode = "markers", 
  text= wt129lvrv$gene, color = wt129lvrv$group) %>% 
  layout(title ="129 WT") %>%
  layout(annotations = a)


p

#Volcano Plot 129 

s129lvrv<- m129KO_LvVsRV.htseq.edgeR
colnames(s129lvrv)
dim(s129lvrv)
s129lvrv <- s129lvrv[c("gene", "PValue", "logFC")]
head(s129lvrv)
s129lvrv["group"]<-"Not Significant"
s129lvrv[which(s129lvrv["PValue"] < 0.05 & abs(s129lvrv["logFC"]) < 1 ),"group"]<- "Significant"
s129lvrv[which(s129lvrv["PValue"] > 0.05 & abs(s129lvrv["logFC"]) > 1 ),"group"] <- "Fold Change"
s129lvrv[which(s129lvrv["PValue"] < 0.05 & abs(s129lvrv["logFC"]) > 1 ),"group"] <- "Significant and Fold Change"#sig p value and fold change
top_peaks1 <- s129lvrv[with(s129lvrv, order(logFC, PValue)),][1:5,]
top_peaks1 <- rbind(top_peaks1, s129lvrv[with(s129lvrv, order(-logFC, PValue)),][1:5,])
# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly
b <- list()
for (i in seq_len(nrow(top_peaks1))) {
  n <- top_peaks1[i, ]
  b[[i]] <- list(
    x = m[["logFC"]],
    y = -log10(m[["PValue"]]),
    text = m[["gene"]],
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 0.5,
    ax = 20,
    ay = -40
  )
}


v <- plot_ly(data = s129lvrv, x = s129lvrv$logFC, y = -log10(s129lvrv$PValue), mode = "markers", 
             text= s129lvrv$gene, color = s129lvrv$group) %>% 
  layout(title ="129 KO") %>%
  layout(annotations = b)


v
