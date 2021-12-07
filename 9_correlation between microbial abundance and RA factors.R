## Only RA samples
bac.clean.ss.RA
bac.clean.ss

##Correlation abundance and RA factors
##Bacteria
##Family
## Calculate the relative abundance
good_family_all <-tax_glom(bac.clean.ss, taxrank = "Family", NArm = FALSE)
##  Melt it into a dataframe 
goodsamps_phy_melt_all <- psmelt(good_family_all)  

goodsamps_phy_melt_all.m<- tidyr::unite(goodsamps_phy_melt_all, "Class", Phylum, Class, sep = "_")
goodsamps_phy_melt_all.m<- tidyr::unite(goodsamps_phy_melt_all.m, "Order", Class, Order, sep = "_")
goodsamps_phy_melt_all.m<- tidyr::unite(goodsamps_phy_melt_all.m, "Family", Order, Family, sep = "_")
head(goodsamps_phy_melt_all.m)

sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt_all.m, select = c("Sample","Family","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
head(sub_goodsamps_phy_melt)
### Merge family sums and the total numbers -> so we can calculate the Relative Abundance

sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
head(sub_phy_melt_totals)

SUBTOTALS <-plyr::ddply(sub_phy_melt_totals, c("Sample", "Family"), summarise, mean_RA = mean(RelAbundance))
head(SUBTOTALS)

library(reshape2)
SUBTOTALS.trans<- dcast(SUBTOTALS, Sample ~ Family, value.var = "mean_RA")
names(SUBTOTALS.trans)[1] <- c("SampleID")

ps1.meta.all <- meta(bac.clean.ss)
ps1.meta.all <- data.frame(ps1.meta.all)
ps1.meta.numeric <- ps1.meta.all %>% select(SampleID, Age, BMI, RA_factor, anti_CCP)


ps1.meta.all.family <- merge(ps1.meta.numeric, SUBTOTALS.trans, by = "SampleID")
row.names(ps1.meta.all.family) <- ps1.meta.all.family$SampleID
ps1.meta.all.family<-ps1.meta.all.family[-c(1)]

cor_5.func <- Hmisc::rcorr(as.matrix(ps1.meta.all.family), type="spearman")
M.func <- cor_5.func$r
p_mat.func <- cor_5.func$P


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

cor.bac.family <- flattenCorrMatrix(M.func,p_mat.func)
cor.bac.family$sig <- stars.pval(cor.bac.family$p)

head(cor.bac.family)

RA_variables<-colnames(ps1.meta.all.family[c(1:4)])
familynames.bac<-colnames(ps1.meta.all.family[-c(1:4)])

cor.bac.family.trim <- subset(cor.bac.family, row%in%RA_variables)
cor.bac.family.trim <- subset(cor.bac.family.trim, column%in%familynames.bac)

cor.bac.family.trim$cor<- round(cor.bac.family.trim$cor, 5)

cor.bac.family.trim.sig<- tidyr::unite(cor.bac.family.trim, "Rho", cor, sig, sep = "")

cor.bac.family.trim.sig.trans <- dcast(cor.bac.family.trim.sig, row ~ column, value.var = "Rho")
cor.bac.family.trim.sig.trans.p <- dcast(cor.bac.family.trim.sig, row ~ column, value.var = "p")

write.xlsx(cor.bac.family.trim.sig.trans, "Relative abundance_Spearman correlation between soil and bacterial family.xlsx")
write.xlsx(cor.bac.family.trim.sig.trans.p, "Relative abundance_Spearman p-values for spearman correlation between soil and bacterial family.xlsx")



##Fungi
## Correlation based on relative abundance
## Calculate the relative abundance
good_family_all.f <-tax_glom(fun.clean.ss, taxrank = "Family", NArm = FALSE)
##  Melt it into a dataframe 
goodsamps_phy_melt_all.f <- psmelt(good_family_all.f)  

goodsamps_phy_melt_all.f.m<- tidyr::unite(goodsamps_phy_melt_all.f, "Class", Phylum, Class, sep = "_")
goodsamps_phy_melt_all.f.m<- tidyr::unite(goodsamps_phy_melt_all.f.m, "Order", Class, Order, sep = "_")
goodsamps_phy_melt_all.f.m<- tidyr::unite(goodsamps_phy_melt_all.f.m, "Family", Order, Family, sep = "_")

sub_goodsamps_phy_melt.f <- subset(goodsamps_phy_melt_all.f.m, select = c("Sample","Family","Abundance"))
TOTALS.f <- plyr::ddply(sub_goodsamps_phy_melt.f, c("Sample"), summarise, total = sum(Abundance))   
head(sub_goodsamps_phy_melt.f)
### Merge family sums and the total numbers -> so we can calculate the Relative Abundance

sub_phy_melt_totals.f <- merge(sub_goodsamps_phy_melt.f, TOTALS, by = "Sample")  
sub_phy_melt_totals.f$RelAbundance <- sub_phy_melt_totals.f$Abundance/sub_phy_melt_totals.f$total  
head(sub_phy_melt_totals.f)

SUBTOTALS.f <-plyr::ddply(sub_phy_melt_totals.f, c("Sample", "Family"), summarise, mean_RA = mean(RelAbundance))
head(SUBTOTALS.f)

SUBTOTALS.trans.f<- dcast(SUBTOTALS.f, Sample ~ Family, value.var = "mean_RA")
names(SUBTOTALS.trans.f)[1] <- c("SampleID")
ps2.meta.all <- meta(fun.clean.ss)
ps2.meta.all.family <- merge(ps2.meta.all, SUBTOTALS.trans.f, by = "SampleID")

ps2.meta.all.family <- ps2.meta.all.family[-c(1:6)]
head(ps2.meta.all.family)

sapply(ps2.meta.all.family, class) ## soil properties are class 'factor'
indx.div <- sapply(ps2.meta.all.family, is.factor)
ps2.meta.all.family[indx.div] <- lapply(ps2.meta.all.family[indx.div], function(x) as.numeric(as.character(x)))
head(ps2.meta.all.family)
sapply(ps2.meta.all.family, class)

cor_5.f <- Hmisc::rcorr(as.matrix(ps2.meta.all.family), type="spearman")
M.f <- cor_5.f$r
p_mat.f <- cor_5.f$P


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

cor.fun.family <- flattenCorrMatrix(M.f,p_mat.f)
cor.fun.family$sig <- stars.pval(cor.fun.family$p)

head(cor.fun.family)

physicochemicals<-colnames(ps2.meta.all.family[c(1:12)])
familynames.fun<-colnames(ps2.meta.all.family[-c(1:12)])

cor.fun.family.trim <- subset(cor.fun.family, row%in%physicochemicals)
cor.fun.family.trim <- subset(cor.fun.family.trim, column%in%familynames.fun)

cor.fun.family.trim$cor<- round(cor.fun.family.trim$cor, 4)

cor.fun.family.trim.sig<- tidyr::unite(cor.fun.family.trim, "Rho", cor, sig, sep = "")

cor.fun.family.trim.sig.trans <- dcast(cor.fun.family.trim.sig, row ~ column, value.var = "Rho")
cor.fun.family.trim.sig.trans.p <- dcast(cor.fun.family.trim.sig, row ~ column, value.var = "p")

write.xlsx(cor.fun.family.trim.sig.trans, "Relative abundance_Spearman correlation between soil and fungal family.xlsx")
write.xlsx(cor.fun.family.trim.sig.trans.p, "Relative abundance_Spearman p-values for spearman correlation between soil and fungal family.xlsx")


##Genus
##Bacteria
good_genus_all <-tax_glom(bac.clean.ss, taxrank = "Genus", NArm = FALSE)
##  Melt it into a dataframe 
goodsamps_phy_melt_all <- psmelt(good_genus_all)  

goodsamps_phy_melt_all$Phylum<-as.character(goodsamps_phy_melt_all$Phylum)
goodsamps_phy_melt_all$Class<-as.character(goodsamps_phy_melt_all$Class)
goodsamps_phy_melt_all$Order<-as.character(goodsamps_phy_melt_all$Order)
goodsamps_phy_melt_all$Family<-as.character(goodsamps_phy_melt_all$Family)
goodsamps_phy_melt_all$Genus<-as.character(goodsamps_phy_melt_all$Genus)

goodsamps_phy_melt_all$Phylum[is.na(goodsamps_phy_melt_all$Phylum)] <- "unidentified"
goodsamps_phy_melt_all$Class[is.na(goodsamps_phy_melt_all$Class)] <- "unidentified"
goodsamps_phy_melt_all$Order[is.na(goodsamps_phy_melt_all$Order)] <- "unidentified"
goodsamps_phy_melt_all$Family[is.na(goodsamps_phy_melt_all$Family)] <- "unidentified"
goodsamps_phy_melt_all$Genus[is.na(goodsamps_phy_melt_all$Genus)] <- "unidentified"

goodsamps_phy_melt_all$Genus2 <- ifelse(goodsamps_phy_melt_all$Genus == "unidentified",ifelse(goodsamps_phy_melt_all$Family=="unidentified",paste0(goodsamps_phy_melt_all$Order,'_',"Unidentified genus"),paste0(goodsamps_phy_melt_all$Family,'_',"Unidentified genus")),paste0(goodsamps_phy_melt_all$Genus))
goodsamps_phy_melt_all$Genus2[grep("Prevotella ", goodsamps_phy_melt_all$Genus2)] <- "Prevotella"
goodsamps_phy_melt_all$Genus2[grep("Ruminococcus ", goodsamps_phy_melt_all$Genus2)] <- "Ruminococcus"
goodsamps_phy_melt_all$Genus2[grep("Coprococcus ", goodsamps_phy_melt_all$Genus2)] <- "Coprococcus"
goodsamps_phy_melt_all$Genus2[grep("Ruminiclostridium ", goodsamps_phy_melt_all$Genus2)] <- "Ruminiclostridium"
goodsamps_phy_melt_all$Genus2[grep("Tyzzerella ", goodsamps_phy_melt_all$Genus2)] <- "Tyzzerella"

# goodsamps_phy_melt_all.m<- tidyr::unite(goodsamps_phy_melt_all, "Class", Phylum, Class, sep = "_")
# goodsamps_phy_melt_all.m<- tidyr::unite(goodsamps_phy_melt_all.m, "Order", Class, Order, sep = "_")
# goodsamps_phy_melt_all.m<- tidyr::unite(goodsamps_phy_melt_all.m, "Family", Order, Family, sep = "_")
# goodsamps_phy_melt_all.m<- tidyr::unite(goodsamps_phy_melt_all.m, "Genus", Family, Genus, sep = "_")
# 
head(goodsamps_phy_melt_all)

sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt_all, select = c("Sample","Genus","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
head(sub_goodsamps_phy_melt)
### Merge genus sums and the total numbers -> so we can calculate the Relative Abundance

sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
head(sub_phy_melt_totals)

SUBTOTALS <-plyr::ddply(sub_phy_melt_totals, c("Sample", "Genus"), summarise, mean_RA = mean(RelAbundance))
head(SUBTOTALS)

SUBTOTALS.trans<- dcast(SUBTOTALS, Sample ~ Genus, value.var = "mean_RA")
names(SUBTOTALS.trans)[1] <- c("SampleID")

ps1.meta.all <- b.meta.ra.edit
ps1.meta.all <- data.frame(ps1.meta.all)
ps1.meta.numeric <- ps1.meta.all %>% select(SampleID, Age, BMI, RA_factor, anti_CCP, CRP, ESR, WBC, Hb, Plt, BUN, Cr, MTX_dose,Total_cholesterol, Triglyceride, HDL)


ps1.meta.all.genus <- merge(ps1.meta.numeric, SUBTOTALS.trans, by = "SampleID")
row.names(ps1.meta.all.genus) <- ps1.meta.all.genus$SampleID
ps1.meta.all.genus<-ps1.meta.all.genus[-c(1)]

cor_5.func <- Hmisc::rcorr(as.matrix(ps1.meta.all.genus), type="spearman")
M.func <- cor_5.func$r
p_mat.func <- cor_5.func$P




cor.bac.genus <- flattenCorrMatrix(M.func,p_mat.func)
cor.bac.genus$sig <- gtools::stars.pval(cor.bac.genus$p)

head(cor.bac.genus)

RA_variables<-colnames(ps1.meta.all.genus[c(1:15)])
genusnames.bac<-colnames(ps1.meta.all.genus[-c(1:15)])

cor.bac.genus.trim <- subset(cor.bac.genus, row%in%RA_variables)
cor.bac.genus.trim <- subset(cor.bac.genus.trim, column%in%genusnames.bac)

cor.bac.genus.trim$cor<- round(cor.bac.genus.trim$cor, 5)
head(cor.bac.genus.trim)


##Only significant
sig_cor.bac.genus.trim <- subset(cor.bac.genus.trim, p < 0.05)
sig_cor.bac.genus.trim.trans <- dcast(sig_cor.bac.genus.trim , row ~ column, value.var = "cor")
sig_cor.bac.genus.trim.trans.p <- dcast(sig_cor.bac.genus.trim , row ~ column, value.var = "p")
write.xlsx(t(sig_cor.bac.genus.trim.trans), "siginificant correlation between bacterial genus and quantitative factors.xlsx")


## all
cor.bac.genus.trim.sig<- tidyr::unite(cor.bac.genus.trim, "Rho", cor, sig, sep = "")

cor.bac.genus.trim.sig.trans <- dcast(cor.bac.genus.trim.sig, row ~ column, value.var = "Rho")
cor.bac.genus.trim.sig.trans.p <- dcast(cor.bac.genus.trim.sig, row ~ column, value.var = "p")
head(t(cor.bac.genus.trim.sig.trans))
t(cor.bac.genus.trim.sig.trans)
write.xlsx(t(cor.bac.genus.trim.sig.trans), "correlation between bacterial genus and quantitative factors.xlsx")

write.xlsx(cor.bac.genus.trim.sig.trans, "Relative abundance_Spearman correlation between RA and bacterial genus.xlsx")
write.xlsx(cor.bac.genus.trim.sig.trans.p, "Relative abundance_Spearman p-values for spearman correlation between RA and bacterial genus.xlsx")


##Linear regression between RA-related factors and microbial RA


##RA factor
ggplot(ps1.meta.all.genus, aes(x=MTX_dose, y=`Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Lachnospiraceae UCG-001`)) +
  xlab('\n MTX dose')+
  ylab("Relative abundance \n") +
  geom_point(size=3, alpha=0.7) +
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  #scale_x_continuous(breaks=seq(0,1,0.2))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred")

#Anti CCP
ggplot(ps1.meta.all.genus, aes(x=anti_CCP, y=`Firmicutes_Erysipelotrichia_Erysipelotrichales_Erysipelotrichaceae_Turicibacter`)) +
  xlab('\n anti CCP')+
  ylab("Relative abundance \n") +
  geom_point(size=3, alpha=0.7) +
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  #scale_x_continuous(breaks=seq(0,1,0.2))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred")


##Fungi
## Correlation based on relative abundance
## Calculate the relative abundance
good_genus_all.f <-tax_glom(fun.clean.ss, taxrank = "Genus", NArm = FALSE)
##  Melt it into a dataframe 
goodsamps_phy_melt_all.f <- psmelt(good_genus_all.f)  

goodsamps_phy_melt_all.f$Phylum<-as.character(goodsamps_phy_melt_all.f$Phylum)
goodsamps_phy_melt_all.f$Class<-as.character(goodsamps_phy_melt_all.f$Class)
goodsamps_phy_melt_all.f$Order<-as.character(goodsamps_phy_melt_all.f$Order)
goodsamps_phy_melt_all.f$Family<-as.character(goodsamps_phy_melt_all.f$Family)
goodsamps_phy_melt_all.f$Genus<-as.character(goodsamps_phy_melt_all.f$Genus)

goodsamps_phy_melt_all.f$Phylum[is.na(goodsamps_phy_melt_all.f$Phylum)] <- "unidentified"
goodsamps_phy_melt_all.f$Class[is.na(goodsamps_phy_melt_all.f$Class)] <- "unidentified"
goodsamps_phy_melt_all.f$Order[is.na(goodsamps_phy_melt_all.f$Order)] <- "unidentified"
goodsamps_phy_melt_all.f$Family[is.na(goodsamps_phy_melt_all.f$Family)] <- "unidentified"
goodsamps_phy_melt_all.f$Genus[is.na(goodsamps_phy_melt_all.f$Genus)] <- "unidentified"

goodsamps_phy_melt_all.f$Genus2 <- ifelse(goodsamps_phy_melt_all.f$Genus == "unidentified",ifelse(goodsamps_phy_melt_all.f$Family=="unidentified",paste0(goodsamps_phy_melt_all.f$Order,'_',"Unidentified genus"),paste0(goodsamps_phy_melt_all.f$Family,'_',"Unidentified genus")),paste0(goodsamps_phy_melt_all.f$Genus))
head(goodsamps_phy_melt_all.f)

sub_goodsamps_phy_melt.f <- subset(goodsamps_phy_melt_all.f, select = c("Sample","Genus","Abundance"))
TOTALS.f <- plyr::ddply(sub_goodsamps_phy_melt.f, c("Sample"), summarise, total = sum(Abundance))   
head(sub_goodsamps_phy_melt.f)
### Merge genus sums and the total numbers -> so we can calculate the Relative Abundance

sub_phy_melt_totals.f <- merge(sub_goodsamps_phy_melt.f, TOTALS, by = "Sample")  
sub_phy_melt_totals.f$RelAbundance <- sub_phy_melt_totals.f$Abundance/sub_phy_melt_totals.f$total  
head(sub_phy_melt_totals.f)

SUBTOTALS.f <-plyr::ddply(sub_phy_melt_totals.f, c("Sample", "Genus"), summarise, mean_RA = mean(RelAbundance))
head(SUBTOTALS.f)

SUBTOTALS.trans.f<- dcast(SUBTOTALS.f, Sample ~ Genus, value.var = "mean_RA")
names(SUBTOTALS.trans.f)[1] <- c("SampleID")
ps2.meta.all <- meta(fun.clean.ss)
ps2.meta.numeric <- ps2.meta.all %>% select(SampleID, Age, BMI, RA_factor, anti_CCP, CRP, ESR, WBC, Hb, Plt, BUN, Cr, MTX_dose,Total_cholesterol, Triglyceride, HDL)


ps2.meta.all.genus <- merge(ps2.meta.numeric, SUBTOTALS.trans.f, by = "SampleID")
row.names(ps2.meta.all.genus) <- ps2.meta.all.genus$SampleID
ps2.meta.all.genus<-ps2.meta.all.genus[-c(1)]

cor_5.f <- Hmisc::rcorr(as.matrix(ps2.meta.all.genus), type="spearman")
M.f <- cor_5.f$r
p_mat.f <- cor_5.f$P


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

cor.fun.genus <- flattenCorrMatrix(M.f,p_mat.f)
cor.fun.genus$sig <- gtools::stars.pval(cor.fun.genus$p)

head(cor.fun.genus)


RA_variables<-colnames(ps2.meta.all.genus[c(1:15)])
genusnames.fun<-colnames(ps2.meta.all.genus[-c(1:15)])

cor.fun.genus.trim <- subset(cor.fun.genus, row%in%RA_variables)
cor.fun.genus.trim <- subset(cor.fun.genus.trim, column%in%genusnames.fun)

cor.fun.genus.trim$cor<- round(cor.fun.genus.trim$cor, 5)
head(cor.fun.genus.trim)

write.csv(cor.fun.genus.trim, "Raw data table for correlation heatmap.csv")

##Only significant
sig_cor.fun.genus.trim <- subset(cor.fun.genus.trim, p < 0.05)
sig_cor.fun.genus.trim.trans <- dcast(sig_cor.fun.genus.trim , row ~ column, value.var = "cor")
sig_cor.fun.genus.trim.trans.p <- dcast(sig_cor.fun.genus.trim , row ~ column, value.var = "p")
write.xlsx(t(sig_cor.fun.genus.trim.trans), "siginificant correlation between fungal genus and quantitative factors.xlsx")


## all
cor.fun.genus.trim.sig<- tidyr::unite(cor.fun.genus.trim, "Rho", cor, sig, sep = "")

cor.fun.genus.trim.sig.trans <- dcast(cor.fun.genus.trim.sig, row ~ column, value.var = "Rho")
cor.fun.genus.trim.sig.trans.p <- dcast(cor.fun.genus.trim.sig, row ~ column, value.var = "p")
head(t(cor.fun.genus.trim.sig.trans))
t(cor.fun.genus.trim.sig.trans)
write.xlsx(t(cor.fun.genus.trim.sig.trans), "correlation between fungal genus and quantitative factors.xlsx")

write.xlsx(cor.fun.genus.trim.sig.trans, "Relative abundance_Spearman correlation between RA and funterial genus.xlsx")
write.xlsx(cor.fun.genus.trim.sig.trans.p, "Relative abundance_Spearman p-values for spearman correlation between RA and funterial genus.xlsx")



## Correlation with normalized abundance
good_genus_all <-tax_glom(bac.clean.nolog, taxrank = "Genus", NArm = FALSE)
##  Melt it into a dataframe 
goodsamps_phy_melt_all <- psmelt(good_genus_all)  

goodsamps_phy_melt_all.m<- tidyr::unite(goodsamps_phy_melt_all, "Class", Phylum, Class, sep = "_")
goodsamps_phy_melt_all.m<- tidyr::unite(goodsamps_phy_melt_all.m, "Order", Class, Order, sep = "_")
goodsamps_phy_melt_all.m<- tidyr::unite(goodsamps_phy_melt_all.m, "Family", Order, Family, sep = "_")
goodsamps_phy_melt_all.m<- tidyr::unite(goodsamps_phy_melt_all.m, "Genus", Family, Genus, sep = "_")
head(goodsamps_phy_melt_all.m)

sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt_all.m, select = c("Sample","Genus","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
head(sub_goodsamps_phy_melt)
### Merge genus sums and the total numbers -> so we can calculate the Relative Abundance

sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
head(sub_phy_melt_totals)

SUBTOTALS <-plyr::ddply(sub_phy_melt_totals, c("Sample", "Genus"), summarise, mean_A = mean(Abundance))
head(SUBTOTALS)

SUBTOTALS.trans<- dcast(SUBTOTALS, Sample ~ Genus, value.var = "mean_A")
names(SUBTOTALS.trans)[1] <- c("SampleID")

ps1.meta.all <- b.meta.ra.edit
ps1.meta.all <- data.frame(ps1.meta.all)
ps1.meta.numeric <- ps1.meta.all %>% select(SampleID, Age, BMI, RA_factor, anti_CCP, CRP, ESR, WBC, Hb, Plt, BUN, Cr, MTX_dose,Total_Cholesterol, Triglyceride, HDL)


ps1.meta.all.genus <- merge(ps1.meta.numeric, SUBTOTALS.trans, by = "SampleID")
row.names(ps1.meta.all.genus) <- ps1.meta.all.genus$SampleID
ps1.meta.all.genus<-ps1.meta.all.genus[-c(1)]

cor_5.func <- Hmisc::rcorr(as.matrix(ps1.meta.all.genus), type="spearman")
M.func <- cor_5.func$r
p_mat.func <- cor_5.func$P


cor.bac.genus <- flattenCorrMatrix(M.func,p_mat.func)
cor.bac.genus$sig <- gtools::stars.pval(cor.bac.genus$p)

head(cor.bac.genus)

RA_variables<-colnames(ps1.meta.all.genus[c(1:15)])
genusnames.bac<-colnames(ps1.meta.all.genus[-c(1:15)])

cor.bac.genus.trim <- subset(cor.bac.genus, row%in%RA_variables)
cor.bac.genus.trim <- subset(cor.bac.genus.trim, column%in%genusnames.bac)

cor.bac.genus.trim$cor<- round(cor.bac.genus.trim$cor, 5)
head(cor.bac.genus.trim)


##Only significant
sig_cor.bac.genus.trim <- subset(cor.bac.genus.trim, p < 0.05)
sig_cor.bac.genus.trim.trans <- dcast(sig_cor.bac.genus.trim , row ~ column, value.var = "cor")
sig_cor.bac.genus.trim.trans.p <- dcast(sig_cor.bac.genus.trim , row ~ column, value.var = "p")
write.xlsx(t(sig_cor.bac.genus.trim.trans), "normalized abundance_siginificant correlation between bacterial genus and quantitative factors.xlsx")


## all
cor.bac.genus.trim.sig<- tidyr::unite(cor.bac.genus.trim, "Rho", cor, sig, sep = "")

cor.bac.genus.trim.sig.trans <- dcast(cor.bac.genus.trim.sig, row ~ column, value.var = "Rho")
cor.bac.genus.trim.sig.trans.p <- dcast(cor.bac.genus.trim.sig, row ~ column, value.var = "p")
head(t(cor.bac.genus.trim.sig.trans))
t(cor.bac.genus.trim.sig.trans)
write.xlsx(t(cor.bac.genus.trim.sig.trans), "correlation between bacterial genus and quantitative factors.xlsx")
