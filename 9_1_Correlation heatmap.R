### Correlation heatmap
###Bacteria
## Relative abundance
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
unique(goodsamps_phy_melt_all$Genus2)
sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt_all, select = c("Sample","Genus2","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
head(sub_goodsamps_phy_melt)
### Merge genus sums and the total numbers -> so we can calculate the Relative Abundance

sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
head(sub_phy_melt_totals)

SUBTOTALS <-plyr::ddply(sub_phy_melt_totals, c("Sample", "Genus2"), summarise, mean_RA = mean(RelAbundance))
head(SUBTOTALS)

SUBTOTALS.trans<- dcast(SUBTOTALS, Sample ~ Genus2, value.var = "mean_RA")
names(SUBTOTALS.trans)[1] <- c("SampleID")

ps1.meta.all <- b.meta.ra.edit
ps1.meta.all <- data.frame(ps1.meta.all)
ps1.meta.numeric <- ps1.meta.all %>% select(SampleID, Age, BMI, RA_factor, anti_CCP, CRP, ESR, WBC, Hb, Plt, BUN, Cr, MTX_dose,Total_cholesterol, Triglyceride, HDL)


ps1.meta.all.genus <- merge(ps1.meta.numeric, SUBTOTALS.trans, by = "SampleID")
row.names(ps1.meta.all.genus) <- ps1.meta.all.genus$SampleID
ps1.meta.all.genus<-ps1.meta.all.genus[-c(1)]

subject.factor<- c("Age", "BMI", "RA_factor", "anti_CCP", "CRP", "ESR", "WBC", "Hb", "Plt", "BUN", "Cr", "MTX_dose","Total_cholesterol", "Triglyceride", "HDL")

###Correlation
remain.cols<- colnames(ps1.meta.all.genus)[-which(colSums(ps1.meta.all.genus)==0)]
length(remain.cols)
genus.cols <- remain.cols[-which(remain.cols %in% subject.factor)]

sub1<-subset(t(ps1.meta.all.genus), rownames(t(ps1.meta.all.genus))%in%remain.cols)

ps1.meta.all.genus.t <-t(sub1)

cor.2 <- Hmisc::rcorr(as.matrix(ps1.meta.all.genus.t), type="spearman")
cor.2.r <- cor.2$r
cor.2.p <- cor.2$P
reorder_cormat <- function(cormat2){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat2)/2)
  hc <- hclust(dd)
  cormat2 <-cormat2[hc$order, hc$order]
}

get_upper_tri <- function(cormat2){
  cormat2[lower.tri(cormat2)]<- NA
  return(cormat2)
}

cormat <- reorder_cormat(cor.2.r)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
head(melted_cormat$Var1)
melted_cormat<-subset(melted_cormat, Var2 %in% subject.factor)
melted_cormat<-subset(melted_cormat, Var1 %in% genus.cols)


##p-value table
cormat.p <- reorder_cormat(cor.2.p)
upper_tri.p <- get_upper_tri(cormat.p)
# Melt the correlation matrix
melted_cormat.p <- melt(upper_tri.p, na.rm = TRUE)

melted_cormat.p<-subset(melted_cormat.p, Var2 %in% subject.factor)
melted_cormat.p<-subset(melted_cormat.p, Var1 %in% genus.cols)
head(melted_cormat.p)
head(melted_cormat)

names(melted_cormat.p)[3] <- "P_value"
melted_cormat.r.p <-merge(melted_cormat, melted_cormat.p, by = c("Var1" = "Var1", "Var2" = "Var2"))
melted_cormat.r.p
melted_cormat.r.p$sig <- gtools::stars.pval(melted_cormat.r.p$P_value)

###Ordering subject factors
melted_cormat.r.p$Var2 <- factor(melted_cormat.r.p$Var2, levels = subject.factor)



ggplot(melted_cormat.r.p, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spaearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(aspect.ratio = 3) + 
  geom_text(aes(Var2, Var1, label = sig), color = "black", size = 4) +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = "top",
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


### Normalized abundance
good_genus_all <-tax_glom(bac.clean.nolog, taxrank = "Genus", NArm = FALSE)
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

sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt_all, select = c("Sample","Genus2","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
head(sub_goodsamps_phy_melt)

sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
head(sub_phy_melt_totals)

SUBTOTALS <-plyr::ddply(sub_phy_melt_totals, c("Sample", "Genus2"), summarise, mean_A = mean(Abundance))
head(SUBTOTALS)

SUBTOTALS.trans<- dcast(SUBTOTALS, Sample ~ Genus2, value.var = "mean_A")
names(SUBTOTALS.trans)[1] <- c("SampleID")

ps1.meta.all <- b.meta.ra.edit
ps1.meta.all <- data.frame(ps1.meta.all)
ps1.meta.numeric <- ps1.meta.all %>% select(SampleID, Age, BMI, RA_factor, anti_CCP, CRP, ESR, WBC, Hb, Plt, BUN, Cr, MTX_dose,Total_cholesterol, Triglyceride, HDL)


ps1.meta.all.genus <- merge(ps1.meta.numeric, SUBTOTALS.trans, by = "SampleID")
row.names(ps1.meta.all.genus) <- ps1.meta.all.genus$SampleID
ps1.meta.all.genus<-ps1.meta.all.genus[-c(1)]

subject.factor<- c("Age", "BMI", "RA_factor", "anti_CCP", "CRP", "ESR", "WBC", "Hb", "Plt", "BUN", "Cr", "MTX_dose","Total_cholesterol", "Triglyceride", "HDL")

###Corelation
remain.cols<- colnames(ps1.meta.all.genus)[-which(colSums(ps1.meta.all.genus)==0)]
length(remain.cols)
genus.cols <- remain.cols[-which(remain.cols %in% subject.factor)]

sub1<-subset(t(ps1.meta.all.genus), rownames(t(ps1.meta.all.genus))%in%remain.cols)

ps1.meta.all.genus.t <-t(sub1)

cor.2 <- Hmisc::rcorr(as.matrix(ps1.meta.all.genus.t), type="spearman")
cor.2.r <- cor.2$r
cor.2.p <- cor.2$P


cormat <- reorder_cormat(cor.2.r)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
head(melted_cormat$Var1)
melted_cormat<-subset(melted_cormat, Var2 %in% subject.factor)
melted_cormat<-subset(melted_cormat, Var1 %in% genus.cols)


##p-value table
cormat.p <- reorder_cormat(cor.2.p)
upper_tri.p <- get_upper_tri(cormat.p)
# Melt the correlation matrix
melted_cormat.p <- melt(upper_tri.p, na.rm = TRUE)

melted_cormat.p<-subset(melted_cormat.p, Var2 %in% subject.factor)
melted_cormat.p<-subset(melted_cormat.p, Var1 %in% genus.cols)
head(melted_cormat.p)
head(melted_cormat)

names(melted_cormat.p)[3] <- "P_value"
melted_cormat.r.p <-merge(melted_cormat, melted_cormat.p, by = c("Var1" = "Var1", "Var2" = "Var2"))
melted_cormat.r.p
melted_cormat.r.p$sig <- gtools::stars.pval(melted_cormat.r.p$P_value)

###Ordering subject factors
melted_cormat.r.p$Var2 <- factor(melted_cormat.r.p$Var2, levels = subject.factor)



ggplot(melted_cormat.r.p, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spaearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(aspect.ratio = 3) + 
  geom_text(aes(Var2, Var1, label = sig), color = "black", size = 4) +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = "top",
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))














###Fungi
## Relative abundance
good_genus_all <-tax_glom(fun.its1.clean.ss, taxrank = "Genus", NArm = FALSE)
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

head(goodsamps_phy_melt_all)
unique(goodsamps_phy_melt_all$Genus2)
sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt_all, select = c("Sample","Genus2","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
head(sub_goodsamps_phy_melt)
### Merge genus sums and the total numbers -> so we can calculate the Relative Abundance

sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
head(sub_phy_melt_totals)

SUBTOTALS <-plyr::ddply(sub_phy_melt_totals, c("Sample", "Genus2"), summarise, mean_RA = mean(RelAbundance))
head(SUBTOTALS)

SUBTOTALS.trans<- dcast(SUBTOTALS, Sample ~ Genus2, value.var = "mean_RA")
names(SUBTOTALS.trans)[1] <- c("SampleID")

ps2.meta.all <- meta(fun.its1.clean.ss)
ps2.meta.all <- data.frame(ps2.meta.all)
ps2.meta.numeric <- ps2.meta.all %>% dplyr::select(SampleID, Age, BMI, RA_factor, anti_CCP, CRP, ESR, WBC, Hb, Plt, BUN, Cr, MTX_dose,Total_cholesterol, Triglyceride, HDL)

ps2.meta.all.genus <- merge(ps2.meta.numeric, SUBTOTALS.trans, by = "SampleID")
row.names(ps2.meta.all.genus) <- ps2.meta.all.genus$SampleID
ps2.meta.all.genus<-ps2.meta.all.genus[-c(1)]

subject.factor<- c("Age", "BMI", "RA_factor", "anti_CCP", "CRP", "ESR", "WBC", "Hb", "Plt", "BUN", "Cr", "MTX_dose","Total_cholesterol", "Triglyceride", "HDL")

###Corelation
remain.cols<- colnames(ps2.meta.all.genus)

genus.cols <- remain.cols[-c(1:15)]
ps2.meta.all.genus.t <- t(ps2.meta.all.genus)
cor.2 <- Hmisc::rcorr(as.matrix(ps2.meta.all.genus), type="spearman")
cor.2.r <- cor.2$r
cor.2.p <- cor.2$P

reorder_cormat <- function(cormat2){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat2)/2)
  hc <- stats::hclust(dd)
  cormat2 <-cormat2[hc$order, hc$order]
}

get_upper_tri <- function(cormat2){
  cormat2[lower.tri(cormat2)]<- NA
  return(cormat2)
}

cor.2.r[is.na(cor.2.r)] <- 0
cormat <- reorder_cormat(cor.2.r)
dim(cor.2.r)
dim(cormat)
upper_tri <- get_upper_tri(cormat)
dim(melted_cormat)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = T)
unique(melted_cormat$Var2)

melted_cormat_anti_BMI <- subset(melted_cormat, Var2 %in% c("anti_CCP", "BMI"))
melted_cormat_others <- subset(melted_cormat, Var1 %in% c("Age","RA_factor", "CRP", "ESR", "WBC", "Hb", "Plt", "BUN", "Cr", "MTX_dose","Total_cholesterol", "Triglyceride", "HDL"))

melted_cormat_others<-melted_cormat_others[,c(2,1,3)]
names(melted_cormat_others)[1] <- "Var1"
names(melted_cormat_others)[2] <- "Var2"

melted_cormat <- rbind(melted_cormat_others, melted_cormat_anti_BMI)


melted_cormat<-subset(melted_cormat, Var2 %in% subject.factor)
melted_cormat<-subset(melted_cormat, Var1 %in% genus.cols)


##p-value table
cor.2.p[is.na(cor.2.p)] <- 1
cormat.p <- reorder_cormat(cor.2.p)
upper_tri.p <- get_upper_tri(cormat.p)
# Melt the correlation matrix
melted_cormat.p <- melt(upper_tri.p, na.rm = TRUE)

melted_cormat.p<-subset(melted_cormat.p, Var2 %in% subject.factor)
melted_cormat.p<-subset(melted_cormat.p, Var1 %in% genus.cols)
head(melted_cormat.p)
head(melted_cormat)

names(melted_cormat.p)[3] <- "P_value"
melted_cormat.r.p <-merge(melted_cormat, melted_cormat.p, by = c("Var1" = "Var1", "Var2" = "Var2"))
melted_cormat.r.p
melted_cormat.r.p$sig <- gtools::stars.pval(melted_cormat.r.p$P_value)

###Ordering subject factors
melted_cormat.r.p$Var2 <- factor(melted_cormat.r.p$Var2, levels = subject.factor)

melted_cormat.r.p.sig <- subset(melted_cormat.r.p, P_value < 0.05)

ggplot(melted_cormat.r.p.sig, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spaearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  #theme(aspect.ratio = 3) + 
  geom_text(aes(Var2, Var1, label = sig), color = "black", size = 4) +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = "top",
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


### Normalized abundance
good_genus_all <-tax_glom(fun.clean.nolog, taxrank = "Genus", NArm = FALSE)
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

sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt_all, select = c("Sample","Genus2","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
head(sub_goodsamps_phy_melt)

sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
head(sub_phy_melt_totals)

SUBTOTALS <-plyr::ddply(sub_phy_melt_totals, c("Sample", "Genus2"), summarise, mean_A = mean(Abundance))
head(SUBTOTALS)

SUBTOTALS.trans<- dcast(SUBTOTALS, Sample ~ Genus2, value.var = "mean_A")
names(SUBTOTALS.trans)[1] <- c("SampleID")

ps2.meta.all <- meta(fun.clean.nolog)
ps2.meta.all <- data.frame(ps2.meta.all)
ps2.meta.numeric <- ps2.meta.all %>% select(SampleID, Age, BMI, RA_factor, anti_CCP, CRP, ESR, WBC, Hb, Plt, BUN, Cr, MTX_dose,Total_cholesterol, Triglyceride, HDL)


ps2.meta.all.genus <- merge(ps2.meta.numeric, SUBTOTALS.trans, by = "SampleID")
row.names(ps2.meta.all.genus) <- ps2.meta.all.genus$SampleID
ps2.meta.all.genus<-ps2.meta.all.genus[-c(1)]

subject.factor<- c("Age", "BMI", "RA_factor", "anti_CCP", "CRP", "ESR", "WBC", "Hb", "Plt", "BUN", "Cr", "MTX_dose","Total_cholesterol", "Triglyceride", "HDL")

###Corelation
remain.cols<- colnames(ps2.meta.all.genus)
genus.cols <- remain.cols[-which(remain.cols %in% subject.factor)]

sub1<-subset(t(ps2.meta.all.genus), rownames(t(ps2.meta.all.genus))%in%remain.cols)

ps2.meta.all.genus.t <-t(sub1)

cor.2 <- Hmisc::rcorr(as.matrix(ps2.meta.all.genus.t), type="spearman")
cor.2.r <- cor.2$r
cor.2.p <- cor.2$P

cor.2.r[is.na(cor.2.r)]<-0
cormat <- reorder_cormat(cor.2.r)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
head(melted_cormat$Var1)
melted_cormat_anti_BMI <- subset(melted_cormat, Var2 %in% c("anti_CCP", "BMI"))
melted_cormat_others <- subset(melted_cormat, Var1 %in% c("Age","RA_factor", "CRP", "ESR", "WBC", "Hb", "Plt", "BUN", "Cr", "MTX_dose","Total_cholesterol", "Triglyceride", "HDL"))

melted_cormat_others<-melted_cormat_others[,c(2,1,3)]
names(melted_cormat_others)[1] <- "Var1"
names(melted_cormat_others)[2] <- "Var2"

melted_cormat <- rbind(melted_cormat_others, melted_cormat_anti_BMI)


melted_cormat<-subset(melted_cormat, Var2 %in% subject.factor)
melted_cormat<-subset(melted_cormat, Var1 %in% genus.cols)


##p-value table
cor.2.p[is.na(cor.2.p)] <-1
cormat.p <- reorder_cormat(cor.2.p)
upper_tri.p <- get_upper_tri(cormat.p)
# Melt the correlation matrix
melted_cormat.p <- melt(upper_tri.p, na.rm = TRUE)

melted_cormat.p<-subset(melted_cormat.p, Var2 %in% subject.factor)
melted_cormat.p<-subset(melted_cormat.p, Var1 %in% genus.cols)
head(melted_cormat.p)
head(melted_cormat)

names(melted_cormat.p)[3] <- "P_value"
melted_cormat.r.p <-merge(melted_cormat, melted_cormat.p, by = c("Var1" = "Var1", "Var2" = "Var2"))
melted_cormat.r.p
melted_cormat.r.p$sig <- gtools::stars.pval(melted_cormat.r.p$P_value)

###Ordering subject factors
melted_cormat.r.p$Var2 <- factor(melted_cormat.r.p$Var2, levels = subject.factor)

melted_cormat.r.p.sig <- subset(melted_cormat.r.p, P_value < 0.05)


ggplot(melted_cormat.r.p.sig, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spaearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(aspect.ratio = 3) + 
  geom_text(aes(Var2, Var1, label = sig), color = "black", size = 4) +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = "top",
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


write.csv(melted_cormat.r.p.sig,"Fig. 5D_spearman correlation.csv")
