## Set color scheme

grey_g<-colorRampPalette(colors=c("#CCCCCC",'#333333'))(3)
blue_g1<-colorRampPalette(colors=c("#00CCFF",'#0066FF'))(3)
blue_g2 <- colorRampPalette(colors=c("#3399CC",'#003366'))(3)
pink_g<-colorRampPalette(colors=c("#FFCCCC",'#FF6699'))(2)
green_g<-colorRampPalette(colors=c("#66CC99",'#006600'))(3)
purple_g<-colorRampPalette(colors=c("#CC99CC",'#663366'))(2)


# TI   orange
# NI   brown
# Gla #CC9900
# AA : 6   blue
# BB : 1  red
# CC : 2   pink
# CCDD : 3   green
# EE : 1   
# FF : 1
# GG : 1
# HHJJ : 2   purple

my_beta_gradient <- c(
  grey_g,
  blue_g1,
  blue_g2,
  'yellow3',
  pink_g,
  green_g,
  '#990000',  '#CCCC99',"#FF9900",
  purple_g)


##Making order right
bac.clean.ss

##for full data

# https://rpubs.com/marschmi/133626
library(plyr)
###Phylum level
box_phylum_cultural <- function(phyloseq2, pl){
  
  good_phylum <-tax_glom(phyloseq2, taxrank = "Phylum")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Diagnosis","Phylum","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
    sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  sub_phy_melt_totals$RelAbundance
   ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 100 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Phylum == pl)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Diagnosis", "MaxAbund")
  
  x <- subset(ps1.meta, Diagnosis=='control')$simp
  y <- subset(ps1.meta, Diagnosis=='RA')$simp
  wilcox.test(x, y, conf.int = TRUE)
  
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    xlab(paste('Mann-Whitney, P-value = ', kw$p.value))+
    ylab("Relative abundance (%) \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

##Use the function
box_phylum(phy.clean.ss, 'Chloroflexi')

box_phylum(phy.clean.ss.2017.f.edit, 'Proteobacteria')

box_phylum_cultural(phy.clean.ss.2017.f, 'Acidobacteria')

box_phylum_cultural(phy.clean.ss.2017.f.edit, 'Armatimonadetes')


#### Class level
box_class <- function(phy.firstpart.5, cla){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Class")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Soil_texture","Class","Abundance"))
  TOTALS <- ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Class == cla)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Soil_texture), max)
  
  colnames(test) <- c("Soil_texture", "MaxAbund")
  
  model <- aov(PercentAbund~Soil_texture,data=sub_phy_melt_totals.proteo)
  TukeyHSD(model)
  hsd <- HSD.test(model, "Soil_texture", group=T)
  hsd$groups
  hsd1 <- tibble::rownames_to_column(hsd$groups, var = 'Soil_texture')
  hsd1 <- data.frame(hsd1)
  hsd1 <- merge(hsd1,test, by = 'Soil_texture')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Soil_texture, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Soil_texture,y=MaxAbund, label=groups), vjust=-1) +
    xlab('')+
    ylab("Relative Abundance ?? \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

## Use the function
box_class(phy.clean.ss.2017.f, 'Nitrospira')
###Class cultural practice
box_phylum_cultural <- function(phy.firstpart.5, phyl){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Phylum")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Phylum","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Phylum == phyl)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Cultural_practice", "MaxAbund")
  
  model <- aov(PercentAbund~Cultural_practice,data=sub_phy_melt_totals.proteo)
  TukeyHSD(model)
  hsd <- HSD.test(model, "Cultural_practice", group=T)
  hsd$groups
  hsd1 <- tibble::rownames_to_column(hsd$groups, var = 'Cultural_practice')
  hsd1 <- data.frame(hsd1)
  hsd1 <- merge(hsd1,test, by = 'Cultural_practice')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Cultural_practice,y=MaxAbund, label=groups), vjust=-1) +
    xlab('')+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}
box_phylum_cultural(bac.clean.ss, "Acidobacteria")

###Class cultural practice
box_class_cultural <- function(phy.firstpart.5, cla){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Class")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Class","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Class == cla)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Cultural_practice", "MaxAbund")
  
  model <- aov(PercentAbund~Cultural_practice,data=sub_phy_melt_totals.proteo)
  TukeyHSD(model)
  hsd <- HSD.test(model, "Cultural_practice", group=T)
  hsd$groups
  hsd1 <- tibble::rownames_to_column(hsd$groups, var = 'Cultural_practice')
  hsd1 <- data.frame(hsd1)
  hsd1 <- merge(hsd1,test, by = 'Cultural_practice')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Cultural_practice,y=MaxAbund, label=groups), vjust=-1) +
    xlab('')+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

## Use the function
box_class_cultural (bac.clean.ss, 'Gammaproteobacteria')
box_class_cultural (fun.clean.ss.2017.f.edit, 'Agaricomycetes')


## Genus ANOVA
###Class cultural practice
box_genus_cultural <- function(phy.firstpart.5, genu){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Genus")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Genus","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Genus == genu)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Cultural_practice", "MaxAbund")
  
  model <- aov(PercentAbund~Cultural_practice,data=sub_phy_melt_totals.proteo)
  TukeyHSD(model)
  hsd <- HSD.test(model, "Cultural_practice", group=T)
  hsd$groups
  hsd1 <- tibble::rownames_to_column(hsd$groups, var = 'Cultural_practice')
  hsd1 <- data.frame(hsd1)
  hsd1 <- merge(hsd1,test, by = 'Cultural_practice')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Cultural_practice,y=MaxAbund, label=groups), vjust=-1) +
    xlab('')+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

box_genus_cultural(arch.clean.ss, "Methanosaeta")
#### Family level
box_family <- function(phy.firstpart.5, fam){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Family", NArm = FALSE)
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Label","Family","Abundance"))
  TOTALS <- ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Family == fam)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Label), max)
  colnames(test) <- c("Label", "MaxAbund")
  
  model <- aov(PercentAbund~Label,data=sub_phy_melt_totals.proteo)
  TukeyHSD(model)
  hsd <- HSD.test(model, "Label", group=T)
  hsd$groups
  hsd1 <- tibble::rownames_to_column(hsd$groups, var = 'Label')
  hsd1 <- data.frame(hsd1)
  hsd1 <- merge(hsd1,test, by = 'Label')
  
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Label, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Label,y=MaxAbund, label=groups), vjust=-1) +
    xlab('')+
    ylab("Relative Abundance ??? \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

## Use the function
box_family(phy.clean.ss.5, 'Cytophagaceae')

colSums(otu_table(phy.thirdpart.5))
tax_table(phy.thirdpart.5)['20dd5fe2a2a31b7d5a1a43a84c16cbf2']
(otu_table(phy.thirdpart.5))[,'F24.1']
df <- goodsamps_phy_melt
df[df$SampleID=='F24.1',]$Abundance

otu_table(good_phylum)[,'F24.1']

#### Genus level
box_genus <- function(phy.firstpart.5, ge){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Genus", NArm = FALSE)
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Label","Genus","Abundance"))
  TOTALS <- ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Genus == ge)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Label), max)
  colnames(test) <- c("Label", "MaxAbund")
  
  model <- aov(PercentAbund~Label,data=sub_phy_melt_totals.proteo)
  TukeyHSD(model)
  hsd <- HSD.test(model, "Label", group=T)
  hsd$groups
  hsd1 <- tibble::rownames_to_column(hsd$groups, var = 'Label')
  hsd1 <- data.frame(hsd1)
  hsd1 <- merge(hsd1,test, by = 'Label')
  
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Label, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Label,y=MaxAbund, label=groups), vjust=-1) +
    xlab('')+
    ylab("Relative Abundance ?? \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

## Use the function
# Bacteria
box_genus(phy.clean.ss.5, 'Pseudonocardia')

box_genus(phy.thirdpart.5, 'Pantoea')

box_genus(phy.firstpart.5, 'Xanthomonas')

# Fungi
box_genus(phy.clean.ss.5, 'Coprinopsis')

box_genus(phy.thirdpart.5, 'Sarocladium')

box_genus(phy.firstpart.5, 'Sarocladium')

###### Box plot at the OTU level
phy.rel <- transform_sample_counts(phy.firstpart.5, function(x) (1000*x)/sum(x))
phy.rel
otu_table(phy.rel)
boxplot_abundance(phy.rel, x = "Label", y = "SH182780.07FU_AB693782_reps", line = NULL, violin = FALSE, na.rm = FALSE,
                  show.points = TRUE) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none") +theme(axis.line = element_line(colour = "black"),
                                                                                                                                                   panel.grid.major = element_blank(),
                                                                                                                                                   panel.grid.minor = element_blank(),
                                                                                                                                                   panel.background = element_blank())

phy.rel.total <- transform_sample_counts(phy.clean.ss.5, function(x) (1000*x)/sum(x))
phy.rel
otu_table(phy.rel)
boxplot_abundance(phy.rel.total, x = "Label", y = "SH182780.07FU_AB693782_reps", line = NULL, violin = FALSE, na.rm = FALSE,
                  show.points = TRUE) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none") +theme(axis.line = element_line(colour = "black"),
                                                                                                                                                    panel.grid.major = element_blank(),
                                                                                                                                                    panel.grid.minor = element_blank(),
                                                                                                                                                   panel.background = element_blank())

##Phylum level
box_phylum_cultural <- function(phy.firstpart.5, phyl){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Phylum")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Phylum","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Phylum == phyl)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Group", "MaxAbund")
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(PercentAbund ~ Cultural_practice, data = sub_phy_melt_totals.proteo)
  kw$p.value
  kw$p.value<- round(kw$p.value, 4)
  
  #library(FSA)
  DT = dunnTest(PercentAbund ~ Cultural_practice,
                data=sub_phy_melt_totals.proteo,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  hsd1 <- merge(dunn,test, by = 'Group')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
    #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

box_phylum_cultural(bac.clean.ss, "Acidobacteria")
box_phylum_cultural(bac.clean.ss, "Proteobacteria")
box_phylum_cultural(bac.clean.ss, "Actinobacteria")
box_phylum_cultural(bac.clean.ss, "Chloroflexi")
box_phylum_cultural(bac.clean.ss, "Nitrospirae")
box_phylum_cultural(bac.clean.ss, "Verrucomicrobia")
box_phylum_cultural(bac.clean.ss, "Planctomycetes")

box_phylum_cultural(arch.clean.ss, "Euryarchaeota")
box_phylum_cultural(arch.clean.ss, "Crenarchaeota")


##Class level
box_class_cultural <- function(phy.firstpart.5, cla){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Class")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Class","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Class == cla)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Group", "MaxAbund")
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(PercentAbund ~ Cultural_practice, data = sub_phy_melt_totals.proteo)
  kw$p.value
  kw$p.value<- round(kw$p.value, 4)
  
  #library(FSA)
  DT = dunnTest(PercentAbund ~ Cultural_practice,
                data=sub_phy_melt_totals.proteo,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  hsd1 <- merge(dunn,test, by = 'Group')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
    xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}


box_class_cultural(bac.clean.ss,"Gammaproteobacteria")
box_class_cultural(fun.clean.ss,"Agaricomycetes")
box_class_cultural(fun.clean.ss,"Tremellomycetes")
box_class_cultural(fun.clean.ss,"Sordariomycetes")
box_class_cultural(fun.clean.ss,"Mortierellomycetes")

### Order level
box_order_cultural <- function(phy.firstpart.5, ord){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Order")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Order","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Order == ord)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Group", "MaxAbund")
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(PercentAbund ~ Cultural_practice, data = sub_phy_melt_totals.proteo)
  kw$p.value
  kw$p.value<- round(kw$p.value, 4)
  
  #library(FSA)
  DT = dunnTest(PercentAbund ~ Cultural_practice,
                data=sub_phy_melt_totals.proteo,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  hsd1 <- merge(dunn,test, by = 'Group')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
    xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

box_order_cultural(bac.clean.ss, "Rhizobiales")


### Family level
box_family_cultural <- function(phy.firstpart.5, fam){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Family")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Family","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Family == fam)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Group", "MaxAbund")
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(PercentAbund ~ Cultural_practice, data = sub_phy_melt_totals.proteo)
  kw$p.value
  kw$p.value<- round(kw$p.value, 4)
  
  #library(FSA)
  DT = dunnTest(PercentAbund ~ Cultural_practice,
                data=sub_phy_melt_totals.proteo,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  hsd1 <- merge(dunn,test, by = 'Group')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
    xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

box_family_cultural(phy.clean.ss, "Rhizobiales")
box_family_cultural(arch.clean.ss, "Rice Cluster I")

library(plyr)

### Genus level
box_genus_diagnosis <- function(phy.firstpart.5, gen){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Genus")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Diagnosis","Genus","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 100 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Genus == gen)
  sub_phy_melt_totals.proteo
  
  x <- subset(sub_phy_melt_totals.proteo, Diagnosis=='control')$PercentAbund
  y <- subset(sub_phy_melt_totals.proteo, Diagnosis=='RA')$PercentAbund
  wilcox.test(x, y, conf.int = TRUE)
  
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Diagnosis, y=PercentAbund, fill = Diagnosis)) + geom_boxplot(width = 0.8) +
    theme_bw() + theme(aspect.ratio=2)+ scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
    geom_point(position='jitter',shape=1, alpha=.5)+ggtitle(gen) +
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Relative abundance (%) \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  
  return(p)
}
vec.reorder

box_genus_diagnosis(bac.clean.ss, "Streptococcus")

box_genus_diagnosis(bac.clean.ss, "Blautia")

box_genus_diagnosis(bac.clean.ss, "Bifidobacterium")

box_genus_diagnosis(bac.clean.ss, "Clostridium sensu stricto 1")

box_genus_diagnosis(bac.clean.ss, "Anaerostipes")

box_genus_diagnosis(bac.clean.ss, "Dorea")

box_genus_diagnosis(bac.clean.ss, "Weissella")

box_genus_diagnosis(bac.clean.ss, "Prevotella")

box_genus_diagnosis(bac.clean.ss, "Bacteroides")




##MTX

good_phylum <-tax_glom(bac.clean.ss, taxrank = "Genus")
good_phylum

##  Melt it into a dataframe 
goodsamps_phy_melt <- psmelt(good_phylum)  
goodsamps_phy_melt
sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "MTX","Genus","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
sub_goodsamps_phy_melt
### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
sub_phy_melt_totals
## Calculate the relative abundance
sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
##  Calculate the Percent Abundance
sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 100 
sub_phy_melt_totals


sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Genus == "Subdoligranulum")
sub_phy_melt_totals.proteo

sub_phy_melt_totals.proteo$MTX <- ifelse(sub_phy_melt_totals.proteo$MTX == "0", "untreated", ifelse(sub_phy_melt_totals.proteo$MTX == "1", "treated", "Control"))
sub_phy_melt_totals.proteo$MTX <-as.factor(sub_phy_melt_totals.proteo$MTX)


max.diversity <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$MTX), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

kw<-kruskal.test(PercentAbund ~ MTX, data = sub_phy_melt_totals.proteo)
kw$p.value

#library(FSA)
DT = dunnTest(PercentAbund ~ MTX,
              data=sub_phy_melt_totals.proteo,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')


p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=MTX, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+ theme(aspect.ratio=1.5)+
  xlab(paste('P = ', kw$p.value))+geom_text(data=hsd1,aes(x=Group,y=MaxDiversity,label=Letter),vjust=-1)+
  ylab("Relative abundance (%) \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")
p




###### Fungal community ######
### Genus level
box_genus_diagnosis <- function(phy.firstpart.5, gen){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Genus")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Diagnosis","Genus","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 100 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Genus == gen)
  sub_phy_melt_totals.proteo
  
  x <- subset(sub_phy_melt_totals.proteo, Diagnosis=='control')$PercentAbund
  y <- subset(sub_phy_melt_totals.proteo, Diagnosis=='RA')$PercentAbund
  wilcox.test(x, y, conf.int = TRUE)
  
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Diagnosis, y=PercentAbund, fill = Diagnosis)) + geom_boxplot(width = 0.8) +
    theme_bw() + theme(aspect.ratio=2)+ scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
    geom_point(position='jitter',shape=1, alpha=.5)+ggtitle(gen) +
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Relative abundance (%) \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  
  return(p)
}
vec.reorder

box_genus_diagnosis(fun.clean.ss, "Candida")

box_genus_diagnosis(fun.clean.ss, "Saccharomyces")

box_genus_diagnosis(fun.clean.ss, "Aspergillus")

box_genus_diagnosis(fun.clean.ss, "Kazachstania")

box_genus_diagnosis(fun.clean.ss, "Penicillium")

box_genus_diagnosis(fun.clean.ss, "Issatchenkia")

box_genus_diagnosis(fun.clean.ss, "Mucor")

box_genus_diagnosis(fun.clean.ss, "Rhodotorula")

box_genus_diagnosis(fun.clean.ss, "Meyerozyma")

box_genus_diagnosis(fun.clean.ss, "Trichosporon")



### Grouped box plot
  good_phylum <-tax_glom(bac.clean.ss, taxrank = "Genus")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  
  goodsamps_phy_melt$Phylum<-as.character(goodsamps_phy_melt$Phylum)
  goodsamps_phy_melt$Class<-as.character(goodsamps_phy_melt$Class)
  goodsamps_phy_melt$Order<-as.character(goodsamps_phy_melt$Order)
  goodsamps_phy_melt$Family<-as.character(goodsamps_phy_melt$Family)
  goodsamps_phy_melt$Genus<-as.character(goodsamps_phy_melt$Genus)
  
  goodsamps_phy_melt$Phylum[is.na(goodsamps_phy_melt$Phylum)] <- "unidentified"
  goodsamps_phy_melt$Class[is.na(goodsamps_phy_melt$Class)] <- "unidentified"
  goodsamps_phy_melt$Order[is.na(goodsamps_phy_melt$Order)] <- "unidentified"
  goodsamps_phy_melt$Family[is.na(goodsamps_phy_melt$Family)] <- "unidentified"
  goodsamps_phy_melt$Genus[is.na(goodsamps_phy_melt$Genus)] <- "unidentified"
  
  goodsamps_phy_melt$Genus2 <- goodsamps_phy_melt$Genus
  goodsamps_phy_melt$Genus2[grep("Prevotella ", goodsamps_phy_melt$Genus2)] <- "Prevotella"
  goodsamps_phy_melt$Genus2[grep("Ruminococcus ", goodsamps_phy_melt$Genus2)] <- "Ruminococcus"
  goodsamps_phy_melt$Genus2[grep("Coprococcus ", goodsamps_phy_melt$Genus2)] <- "Coprococcus"
  goodsamps_phy_melt$Genus2[grep("Ruminiclostridium ", goodsamps_phy_melt$Genus2)] <- "Ruminiclostridium"
  goodsamps_phy_melt$Genus2[grep("Tyzzerella ", goodsamps_phy_melt$Genus2)] <- "Tyzzerella"
  
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Diagnosis","Genus2","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  
  ## filtering genus based on relative abundance
  
  sub_phy_melt_totals.2 <- sub_phy_melt_totals %>% group_by(Genus2) %>% summarise(meanRelAbundance = mean(RelAbundance))
  sub_phy_melt_totals.2 <- subset(sub_phy_melt_totals.2, sub_phy_melt_totals.2$meanRelAbundance >0.01)
  sub_phy_melt_totals <- subset(sub_phy_melt_totals, sub_phy_melt_totals$Genus2 %in% sub_phy_melt_totals.2$Genus2)
  ord <- sub_phy_melt_totals %>% group_by(Genus2) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
  vec <- ord$Genus2
  
  sub_phy_melt_totals$Genus2 <- factor(sub_phy_melt_totals$Genus2, levels = rev(vec))
  
  library(ggpubr)
  p <- ggboxplot(data = sub_phy_melt_totals, x="Genus2", y="RelAbundance", fill = "Diagnosis") +
    theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
    #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
    ylab("Relative abundance\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")+
   stat_compare_means(aes(group = Diagnosis), label = "p.signif", method = "wilcox.test")
  #stat_compare_means(,aes(label = ..p.signif..), method = "wilcox.test")
p

bac.genus.group.p<-compare_means(RelAbundance ~ Diagnosis, sub_phy_melt_totals, method = "wilcox.test",group.by="Genus2")
write.xlsx(bac.genus.group.p, "Bacterial genus comparison_above 1% relative abundance.xlsx")
dev.off()


##fungi
good_phylum <-tax_glom(fun.clean.ss, taxrank = "Genus")
good_phylum

##  Melt it into a dataframe 
goodsamps_phy_melt <- psmelt(good_phylum)  
goodsamps_phy_melt

goodsamps_phy_melt$Phylum<-as.character(goodsamps_phy_melt$Phylum)
goodsamps_phy_melt$Class<-as.character(goodsamps_phy_melt$Class)
goodsamps_phy_melt$Order<-as.character(goodsamps_phy_melt$Order)
goodsamps_phy_melt$Family<-as.character(goodsamps_phy_melt$Family)
goodsamps_phy_melt$Genus<-as.character(goodsamps_phy_melt$Genus)

goodsamps_phy_melt$Phylum[is.na(goodsamps_phy_melt$Phylum)] <- "unidentified"
goodsamps_phy_melt$Class[is.na(goodsamps_phy_melt$Class)] <- "unidentified"
goodsamps_phy_melt$Order[is.na(goodsamps_phy_melt$Order)] <- "unidentified"
goodsamps_phy_melt$Family[is.na(goodsamps_phy_melt$Family)] <- "unidentified"
goodsamps_phy_melt$Genus[is.na(goodsamps_phy_melt$Genus)] <- "unidentified"

sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Diagnosis","Genus","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
sub_goodsamps_phy_melt
### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
sub_phy_melt_totals
## Calculate the relative abundance
sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  

## filtering genus based on relative abundance

sub_phy_melt_totals.2 <- sub_phy_melt_totals %>% group_by(Genus) %>% summarise(meanRelAbundance = mean(RelAbundance))
sub_phy_melt_totals.2 <- subset(sub_phy_melt_totals.2, sub_phy_melt_totals.2$meanRelAbundance >0.004)
sub_phy_melt_totals <- subset(sub_phy_melt_totals, sub_phy_melt_totals$Genus %in% sub_phy_melt_totals.2$Genus)
ord <- sub_phy_melt_totals %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus

sub_phy_melt_totals$Genus <- factor(sub_phy_melt_totals$Genus, levels = rev(vec))

p <- ggboxplot(data = sub_phy_melt_totals, x="Genus", y="RelAbundance", fill = "Diagnosis") +
  theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  ylab("Relative abundance\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+
  stat_compare_means(aes(group = Diagnosis), label = "p.signif", method = "wilcox.test")

p

fun.genus.group.p<-compare_means(RelAbundance ~ Diagnosis, sub_phy_melt_totals, method = "wilcox.test",group.by="Genus")
write.xlsx(fun.genus.group.p, "fungal genus comparison_above 1% relative abundance.xlsx")
dev.off()

###its1
##fungi
good_phylum <-tax_glom(fun.its1.clean.ss, taxrank = "Genus")
good_phylum

##  Melt it into a dataframe 
goodsamps_phy_melt <- psmelt(good_phylum)  
goodsamps_phy_melt

goodsamps_phy_melt$Phylum<-as.character(goodsamps_phy_melt$Phylum)
goodsamps_phy_melt$Class<-as.character(goodsamps_phy_melt$Class)
goodsamps_phy_melt$Order<-as.character(goodsamps_phy_melt$Order)
goodsamps_phy_melt$Family<-as.character(goodsamps_phy_melt$Family)
goodsamps_phy_melt$Genus<-as.character(goodsamps_phy_melt$Genus)

goodsamps_phy_melt$Phylum[is.na(goodsamps_phy_melt$Phylum)] <- "unidentified"
goodsamps_phy_melt$Class[is.na(goodsamps_phy_melt$Class)] <- "unidentified"
goodsamps_phy_melt$Order[is.na(goodsamps_phy_melt$Order)] <- "unidentified"
goodsamps_phy_melt$Family[is.na(goodsamps_phy_melt$Family)] <- "unidentified"
goodsamps_phy_melt$Genus[is.na(goodsamps_phy_melt$Genus)] <- "unidentified"

sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Diagnosis","Genus","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
sub_goodsamps_phy_melt
### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
sub_phy_melt_totals
## Calculate the relative abundance
sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  

## filtering genus based on relative abundance
sub_phy_melt_totals.2 <- sub_phy_melt_totals %>% group_by(Genus) %>% summarise(meanRelAbundance = mean(RelAbundance))
sub_phy_melt_totals.2 <- subset(sub_phy_melt_totals.2, sub_phy_melt_totals.2$meanRelAbundance >0.006)
sub_phy_melt_totals <- subset(sub_phy_melt_totals, sub_phy_melt_totals$Genus %in% sub_phy_melt_totals.2$Genus)
ord <- sub_phy_melt_totals %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus

sub_phy_melt_totals$Genus <- factor(sub_phy_melt_totals$Genus, levels = rev(vec))

p <- ggboxplot(data = sub_phy_melt_totals, x="Genus", y="RelAbundance", fill = "Diagnosis") +
  theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  ylab("Relative abundance\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+
  stat_compare_means(aes(group = Diagnosis), label = "p.signif", method = "wilcox.test")

p

fun.genus.group.p<-compare_means(RelAbundance ~ Diagnosis, sub_phy_melt_totals, method = "wilcox.test",group.by="Genus")
write.xlsx(fun.genus.group.p, "fungal ITS1 genus comparison_top 21 relative abundance.xlsx")
dev.off()







##MTX

good_phylum <-tax_glom(fun.clean.ss, taxrank = "Genus")
good_phylum

##  Melt it into a dataframe 
goodsamps_phy_melt <- psmelt(good_phylum)  
goodsamps_phy_melt
sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "MTX","Genus","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
sub_goodsamps_phy_melt
### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
sub_phy_melt_totals
## Calculate the relative abundance
sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
##  Calculate the Percent Abundance
sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 100 
sub_phy_melt_totals


sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Genus == "Subdoligranulum")
sub_phy_melt_totals.proteo

sub_phy_melt_totals.proteo$MTX <- ifelse(sub_phy_melt_totals.proteo$MTX == "0", "untreated", ifelse(sub_phy_melt_totals.proteo$MTX == "1", "treated", "Control"))
sub_phy_melt_totals.proteo$MTX <-as.factor(sub_phy_melt_totals.proteo$MTX)


max.diversity <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$MTX), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

kw<-kruskal.test(PercentAbund ~ MTX, data = sub_phy_melt_totals.proteo)
kw$p.value

#library(FSA)
DT = dunnTest(PercentAbund ~ MTX,
              data=sub_phy_melt_totals.proteo,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')


p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=MTX, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+ theme(aspect.ratio=1.5)+
  xlab(paste('P = ', kw$p.value))+geom_text(data=hsd1,aes(x=Group,y=MaxDiversity,label=Letter),vjust=-1)+
  ylab("Relative abundance (%) \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")
p



##treatment2

good_phylum <-tax_glom(fun.clean.ss, taxrank = "Genus")
good_phylum

##  Melt it into a dataframe 
goodsamps_phy_melt <- psmelt(good_phylum)  
goodsamps_phy_melt
sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "treatment2","Genus","Abundance"))
TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
sub_goodsamps_phy_melt
### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
sub_phy_melt_totals
## Calculate the relative abundance
sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
##  Calculate the Percent Abundance
sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 100 
sub_phy_melt_totals


sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Genus == "Candida")
sub_phy_melt_totals.proteo

max.diversity <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$treatment2), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

sub_phy_melt_totals.proteo$treatment2 <- as.factor(sub_phy_melt_totals.proteo$treatment2)

kw<-kruskal.test(PercentAbund ~ treatment2, data = sub_phy_melt_totals.proteo)
kw$p.value

#library(FSA)
DT = dunnTest(PercentAbund ~ treatment2,
              data=sub_phy_melt_totals.proteo,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')


p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=treatment2, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+ theme(aspect.ratio=1.5)+
  xlab(paste('P = ', kw$p.value))+geom_text(data=hsd1,aes(x=Group,y=MaxDiversity,label=Letter),vjust=-1)+
  ylab("Relative abundance (%) \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")
p
