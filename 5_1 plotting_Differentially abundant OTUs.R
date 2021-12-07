

## let's plot MA plot

get_resSig <- function(phyloseq1){
  zero.clean.filt <- phyloseq::filter_taxa(phyloseq1, function(x) sum(x) != 0, TRUE)
  sum(taxa_sums(zero.clean.filt) == 0)
  
  ## CODE for CSS normalization using preloaded data
  sort(sample_sums(zero.clean.filt))
  
  obj <- phyloseq_to_metagenomeSeq(zero.clean.filt)
  obj 
  
  ## fitZig sample
  obj <-  cumNorm(obj, p = cumNormStatFast(obj))
  normFactor <-  normFactors(obj)
  normFactor <-  log2(normFactor/median(normFactor) + 1)
  settings <-  zigControl(maxit = 30, verbose = TRUE)
  
  Type  <-  pData(obj)$Diagnosis
  
  mod  <-  model.matrix(~Type)
  colnames(mod)  <-  levels(Type)
  colnames(mod)
  
  res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)
  
  
  zigFit = res$fit
  finalMod= res$fit$design
  
  
  contrast.matrix.1 =makeContrasts(control - RA, levels = finalMod)
  
  
  fit2_diagnosis = contrasts.fit(zigFit, contrasts=contrast.matrix.1) 
  
  fit3 = eBayes(fit2_diagnosis)
  
  topTable(fit3, coef="control - RA")
  
  res <- topTable(fit3,coef=1,adjust="fdr",n=Inf)
  head(res)

  log2AverageAbundance <- psmelt(phyloseq1) %>% group_by(OTU) %>% summarise(log2AverageAbundance=log2(mean(Abundance)))
  log2AverageAbundance
  Ta <- psmelt(phyloseq1) %>% group_by(OTU) %>% select(Phylum,Class,Order,Family, Genus)
  Ta <- unique(Ta)
  Ta
  
  resSig = res[!is.na(res$adj.P.Val), ]
  resSig = data.frame(resSig)
  head(resSig)
  resSig <- tibble::rownames_to_column(resSig, 'OTU')
  resSig <- left_join(resSig, log2AverageAbundance,by= c('OTU','OTU'))
  resSig <- left_join(resSig,Ta,by=c('OTU','OTU'))
  return(resSig)
}

resSig <- get_resSig(bac.clean.ss)

head(resSig)

resSig$significant = (resSig$adj.P.Val < .05)
resSig$significant = as.factor(resSig$significant)
resSig$significant[is.na(resSig$significant)] = F

resSig$significant2 <- 0
resSig$significant2[which(resSig$logFC > 0 & resSig$adj.P.Val < .05)] <- "Control"
resSig$significant2[which(resSig$logFC < 0 & resSig$adj.P.Val < .05)] <- "RA"
resSig$significant2[which(resSig$significant2 == "0")] <- "Non-significant"

resSig$Genus <- as.character(resSig$Genus)
resSig$Genus[is.na(resSig$Genus)] <- "unidentified"



ggplot(resSig, aes(x=log10(AveExpr), y=logFC, color=significant2)) +
  geom_point(size = 2.5) + theme(aspect.ratio=1)+
  geom_hline(yintercept = 0, linetype='dashed', color='black') +
  # stat_smooth(se = FALSE, method = "loess", color = "red3") +
  scale_color_manual(values=c("Non-significant"="light grey","Control" = "#6699CC", "RA" = "#CC9900"))+
  xlab('\n log10(Average abundance)')+
  ylab("log2(Fold Change)\n") +theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


resSig2 <-resSig
resSig2$Family2 <-0
resSig2$Family2[which(resSig2$significant2 == "Non-significant")] <- "Non-significant"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Lachnospiraceae")]<- "Lachnospiraceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Ruminococcaceae")]<- "Ruminococcaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Bifidobacteriaceae")]<- "Bifidobacteriaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Streptococcaceae")]<- "Streptococcaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Peptostreptococcaceae")]<- "Peptostreptococcaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Erysipelotrichaceae")]<- "Erysipelotrichaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Clostridiaceae 1")]<- "Clostridiaceae 1"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Lactobacillaceae")]<- "Lactobacillaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Enterococcaceae")]<- "Enterococcaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Coriobacteriaceae")]<- "Coriobacteriaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Enterobacteriaceae")]<- "Enterobacteriaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Bacteroidaceae")]<- "Bacteroidaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Leuconostocaceae")]<- "Leuconostocaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Veillonellaceae")]<- "Veillonellaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Eggerthellaceae")]<- "Eggerthellaceae"
resSig2$Family2[which(resSig2$significant2 != "Non-significant" & resSig2$Family == "Prevotellaceae")]<- "Prevotellaceae"
resSig2$Family2[which(resSig2$Family2== "0")] <- "Low abundance"


resSig2$Family2 <- factor(resSig2$Family2, levels = c("Non-significant",vec.reorder.fam)) 

ggplot(resSig2, aes(x=log10(AveExpr), y=logFC, color=Family2)) +
  geom_point(size = 2.5) + theme(aspect.ratio=1)+
  geom_hline(yintercept = 0, linetype='dashed', color='black') +
  # stat_smooth(se = FALSE, method = "loess", color = "red3") +
  scale_color_manual(values=c("Non-significant"= "light grey","Lachnospiraceae"= "#4ECDC4","Ruminococcaceae"= "#1A535C",
                              "Bifidobacteriaceae"= "#0E103D","Streptococcaceae"= "#69306D","Peptostreptococcaceae"= "#A5668B",
                              "Erysipelotrichaceae"= "#D3BCC0","Clostridiaceae 1"= "#F2D7EE", "Lactobacillaceae"= "#04471C",
                              "Enterococcaceae"= "#E55812","Coriobacteriaceae"= "#95C623","Enterobacteriaceae"= "#5F7FC7",
                              "Bacteroidaceae"= "#ED7B84","Leuconostocaceae"= "#F92A82","Veillonellaceae"= "#CCA43B","Eggerthellaceae"= "#6D9F71",
                              "Prevotellaceae"= "#242F40","Low abundance" = "#999999"))+
  xlab('\n log10(Average abundance)')+
  ylab("log2(Fold Change)\n") +theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


length(resSig2$Genus[which(resSig2$significant2 == "Control" & resSig2$Family == "Lachnospiraceae")]) #115
length(resSig2$Genus[which(resSig2$significant2 == "RA" & resSig2$Family == "Lachnospiraceae")]) #70

intersect(resSig2$Genus[which(resSig2$significant2 == "Control" & resSig2$Family == "Lachnospiraceae")], resSig2$Genus[which(resSig2$significant2 == "RA" & resSig2$Family == "Lachnospiraceae")])


length(resSig2$Genus[which(resSig2$significant2 == "Control" & resSig2$Family == "Ruminococcaceae")]) #108
length(resSig2$Genus[which(resSig2$significant2 == "RA" & resSig2$Family == "Ruminococcaceae")]) #68

intersect(resSig2$Genus[which(resSig2$significant2 == "Control" & resSig2$Family == "Ruminococcaceae")], resSig2$Genus[which(resSig2$significant2 == "RA" & resSig2$Family == "Ruminococcaceae")])


length(resSig2$Genus[which(resSig2$significant2 == "Control" & resSig2$Family == "Bifidobacteriaceae")]) #8
length(resSig2$Genus[which(resSig2$significant2 == "RA" & resSig2$Family == "Bifidobacteriaceae")]) #0

intersect(resSig2$Genus[which(resSig2$significant2 == "Control" & resSig2$Family == "Bifidobacteriaceae")], resSig2$Genus[which(resSig2$significant2 == "RA" & resSig2$Family == "Bifidobacteriaceae")])


length(resSig2$Genus[which(resSig2$significant2 == "Control" & resSig2$Family == "Streptococcaceae")]) #16
length(resSig2$Genus[which(resSig2$significant2 == "RA" & resSig2$Family == "Streptococcaceae")]) #6

intersect(resSig2$Genus[which(resSig2$significant2 == "Control" & resSig2$Family == "Streptococcaceae")], resSig2$Genus[which(resSig2$significant2 == "RA" & resSig2$Family == "Streptococcaceae")])


length(resSig2$Genus[which(resSig2$significant2 == "Control" & resSig2$Family == "Peptostreptococcaceae")]) #16
length(resSig2$Genus[which(resSig2$significant2 == "RA" & resSig2$Family == "Peptostreptococcaceae")]) #2

intersect(resSig2$Genus[which(resSig2$significant2 == "Control" & resSig2$Family == "Peptostreptococcaceae")], resSig2$Genus[which(resSig2$significant2 == "RA" & resSig2$Family == "Peptostreptococcaceae")])


length(resSig2$OTU[which(resSig2$significant2 == "Control")]) #430
length(resSig2$OTU[which(resSig2$significant2 == "RA")]) #307

intersect(resSig2$Genus[which(resSig2$significant2 == "Control" & resSig2$Family == "Peptostreptococcaceae")], resSig2$Genus[which(resSig2$significant2 == "RA" & resSig2$Family == "Peptostreptococcaceae")])



####Fungal community
resSig.f <- get_resSig(fun.clean.ss)

head(resSig.f)

resSig.f$significant = (resSig.f$adj.P.Val < .05)
resSig.f$significant = as.factor(resSig.f$significant)
resSig.f$significant[is.na(resSig.f$significant)] = F

resSig.f$significant2 <- 0
resSig.f$significant2[which(resSig.f$logFC > 0 & resSig.f$adj.P.Val < .05)] <- "Control"
resSig.f$significant2[which(resSig.f$logFC < 0 & resSig.f$adj.P.Val < .05)] <- "RA"
resSig.f$significant2[which(resSig.f$significant2 == "0")] <- "Non-significant"

resSig.f$Genus <- as.character(resSig.f$Genus)
resSig.f$Genus[is.na(resSig.f$Genus)] <- "unidentified"



ggplot(resSig.f, aes(x=log10(AveExpr), y=logFC, color=significant2)) +
  geom_point(size = 2.5) + theme(aspect.ratio=1)+
  geom_hline(yintercept = 0, linetype='dashed', color='black') +
  # stat_smooth(se = FALSE, method = "loess", color = "red3") +
  scale_color_manual(values=c("Non-significant"="light grey","Control" = "#6699CC", "RA" = "#CC9900"))+
  xlab('\n log10(Average abundance)')+
  ylab("log2(Fold Change)\n") +theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

length(resSig.f$significant2[which(resSig.f$significant2 == "Control")]) #345
length(resSig.f$significant2[which(resSig.f$significant2 == "RA")]) #238

