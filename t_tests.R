#code snippet for calculating mean & SEM, running t-tests, and Bonferroni correction

#df_coding and df_intron are data frames containing position (Pos) and icSHAPE enrichment score (Enrich) 
#for spliced junctions and retained introns, respectively

sem <- function(datavec){
  datavec_sem <- sd(datavec, na.rm=TRUE) /
    sqrt(length(datavec[!is.na(datavec)]))
  return(datavec_sem)
}

poscount <- length(unique(df_coding$Pos))

for (i in unique(df_coding$Pos) ) {
  subdata_coding <- subset(df_coding,Pos == i)$Enrich  
  subdata_intron <- subset(df_intron,Pos == i)$Enrich
  t_out <- t.test(subdata_coding,subdata_intron)
  x <- paste(i,t_out$estimate[[1]],t_out$estimate[[2]],sem(subdata_coding),sem(subdata_intron),t_out$p.value*poscount,sep="\t")
  write.table(x,"output_file.tsv",append = TRUE,quote = FALSE,sep="\t",row.names = FALSE,col.names = FALSE)
}

