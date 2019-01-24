abs.129wt.129ko<- subset(`129WT_129KO_Abs.rsem.edgeR`, gene=="Plcb1" )



View(abs.129wt.129ko)


write.table(abs.129wt.129ko, "~/Desktop/Output folder/abs.129wt.129ko.txt", append = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)


lv.129wt.129ko<-subset(`129WT_129KO_LV.rsem.edgeR`,gene=="Plcb1")

View (lv.129wt.129ko)


write.table(lv.129wt.129ko, "~/Desktop/Output folder/lv129wtko.txt", append = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)


rv.129wt.129ko<-subset(`129WT_129KO_RV.rsem.edgeR`,gene=="Plcb1")



View(rv.129wt.129ko)

write.table(rv.129wt.129ko, "~/Desktop/Output folder/rv.129wt.129ko.txt", append = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)




quad.129wt.129ko<-subset(`129WT_129KO_QUAD.rsem.edgeR`,gene=="Plcb1")


View(quad.129wt.129ko)

write.table(rv.129wt.129ko, "~/Desktop/Output folder/quad.129wt.129ko.txt", append = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)


