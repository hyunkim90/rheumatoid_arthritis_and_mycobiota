
## Ridge plot


## testing


## get imp.tax
get_imp_tax <- function(phy.clean.ss.5, classifier){
  imp <- data.frame(importance(classifier))
  imp <- tibble::rownames_to_column(imp, 'OTU')
  imp.desc <- imp %>% arrange(desc(MeanDecreaseGini))
  tax <- tax_table(phy.clean.ss.5)
  tax <- as.data.frame(tax)
  tax
  tax <- tibble::rownames_to_column(tax,'OTU')
  
  df.clean.ss.5 <- psmelt(phy.clean.ss.5)
  head(df.clean.ss.5)
  abun <- df.clean.ss.5 %>% group_by(OTU) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
  imp.tax <- left_join(imp.desc, tax, by=c('OTU', 'OTU')) %>% left_join(abun, by=c('OTU','OTU'))
  return(imp.tax)
}

#Subset of Conventional and no fertilization 
imp_tax_bac <- get_imp_tax(bac.clean.ss, classifier_diagnosis)

imp_tax_fun <- get_imp_tax(fun.clean.ss, classifier_diagnosis.f)

imp_tax_fun <- get_imp_tax(fun.its1.clean.ss, classifier_diagnosis.f)


#head(mydata)
imp_tax_bac.edit <- subset(imp_tax_bac, OTU != "result1")
head(imp_tax_bac.edit)
imp.arrange.b <- imp_tax_bac.edit %>% arrange(desc(MeanDecreaseGini))
imp.arrange.b$rank <- paste0('RF',1:dim(imp_tax_bac.edit)[1])
another.b <- imp.arrange.b[c('OTU','rank')]


imp_tax_fun.edit <- subset(imp_tax_fun, OTU != "result2")

imp.arrange.f <- imp_tax_fun.edit %>% arrange(desc(MeanDecreaseGini))
imp.arrange.f$rank <- paste0('RF',1:dim(imp_tax_fun.edit)[1])
another.f <- imp.arrange.f[c('OTU','rank')]
head(another.f)




##Construct the plot object

imp.total.arrange <- imp_tax_bac %>% arrange(desc(total))
imp.total.arrange$number <- 1:dim(imp_tax_bac)[1]
imp.total.arrange
num <- imp.total.arrange[c('OTU','number')]
num 

##Construct the plot object
#bac
Ab <- psmelt(bac.clean.ss) %>% group_by(OTU) %>% summarise(Abun=log(sum(Abundance)))
Ta <- psmelt(bac.clean.ss) %>% group_by(OTU) %>% select(OTU, Phylum,Class,Order,Family, Genus,Species)
Ta <- unique(Ta)


imp.total.arrange <- imp_tax_fun %>% arrange(desc(total))
imp.total.arrange$number <- 1:dim(imp_tax_fun)[1]
imp.total.arrange
num <- imp.total.arrange[c('OTU','number')]
num 

#fun
Ab <- psmelt(fun.its1.clean.ss) %>% group_by(OTU) %>% summarise(Abun=log(sum(Abundance)))
Ta <- psmelt(fun.its1.clean.ss) %>% group_by(OTU) %>% select(OTU, Phylum,Class,Order,Family, Genus,Species)
Ta <- unique(Ta)



## indexing the enrichment classification
##Bacteria
resSig <- read.table("lefse_bacterial_otu_effect_size_2.tsv", sep = '\t', header =T)
imp_tax_bac.id<-merge(imp_tax_bac, bac.list.OTU_id, by = "OTU")

mydata.lefse <- merge(imp_tax_bac.id,resSig, by = "OTU" )
mydata.lefse <- mydata.lefse %>% mutate(threshold = ifelse(Significant =='control', 'Control', ifelse(Significant =='RA',"RA", "Non-differential")))
mydata.lefse$threshold[is.na(mydata.lefse$threshold)] <- "Non-differential"

sub.mydata.lefse <- mydata.lefse %>% arrange(desc(MeanDecreaseGini)) %>% head(70) #b:head(70), f:head(80)
unique(sub.mydata.lefse$threshold)
# }
write.xlsx(mydata.lefse, "RF prediction_bacteria.xlsx")

###Fungi
resSig.fun <- read.table("lefse_fungal_otu_effect_size.tsv", sep = '\t', header =T)
imp_tax_fun.id<-merge(imp_tax_fun, fun.list.OTU_id, by = "OTU")

mydata.lefse.fun <- merge(imp_tax_fun.id,resSig.fun, by = "OTU_id" )
mydata.lefse.fun <- mydata.lefse.fun %>% mutate(threshold = ifelse(Significant =='control', 'Control', ifelse(Significant =='RA',"RA", "Non-differential")))
mydata.lefse.fun$threshold[is.na(mydata.lefse.fun$threshold)] <- "Non-differential"

sub.mydata.lefse.fun <- mydata.lefse.fun %>% arrange(desc(MeanDecreaseGini)) %>% head(70) #b:head(70), f:head(80)
unique(sub.mydata.lefse.fun$threshold)
write.xlsx(mydata.lefse.fun, "RF prediction_ITS2.xlsx")
##ITS1
resSig.fun <- read.table("lefse_fungal_its1_otu_effect_size.tsv", sep = '\t', header =T)
imp_tax_fun.id<-merge(imp_tax_fun, fun.list.its1.OTU_id, by = "OTU")

mydata.lefse.fun <- merge(imp_tax_fun.id,resSig.fun, by = "OTU_id" )
mydata.lefse.fun <- mydata.lefse.fun %>% mutate(threshold = ifelse(Significant =='control', 'Control', ifelse(Significant =='RA',"RA", "Non-differential")))
mydata.lefse.fun$threshold[is.na(mydata.lefse.fun$threshold)] <- "Non-differential"

sub.mydata.lefse.fun <- mydata.lefse.fun %>% arrange(desc(MeanDecreaseGini)) %>% head(40)
unique(sub.mydata.lefse.fun$threshold)


write.xlsx(mydata.lefse.fun, "RF prediction_ITS1.xlsx")

## get df.ridge
get_df_ridge <- function(phy.clean.ss.5, imp.tax, sub.mydata,n){
    df.otu <- phy.clean.ss.5 %>% psmelt()
  head(df.otu)
  # we need to group by samples
  df.otu.rel <- df.otu %>%  
    group_by(Sample) %>%                         # Filter out at absolute read of 20       
    mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance
  
  df.otu.rel
  imp.top20 <- imp.tax %>% arrange(desc(MeanDecreaseGini)) %>% head(n) #b: 70, f: 30
  
  imp20 <- imp.top20$OTU
  imp20
  
  imp.top20
  
  df.selected.rel <- df.otu.rel %>% filter(OTU %in% imp20)
  df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
  str(df.ridge)
  b.thresh <- sub.mydata[c('OTU','threshold')]
  b.thresh
  df.ridge_2 <- left_join(df.ridge, b.thresh, by=c('OTU','OTU'))
  
  unique(df.ridge$OTU)
  unique(b.thresh$OTU)
  
  df.ridge_2
  str(df.ridge_2)
  dim(df.ridge_2)
  43*20
  
  return(df.ridge_2)
}
df.ridge_2.bac.lefse<- get_df_ridge(bac.clean.ss, imp_tax_bac, sub.mydata.lefse, 70)

df.ridge_2.fun.lefse<- get_df_ridge(fun.clean.ss, imp_tax_fun, sub.mydata.lefse.fun,70)

df.ridge_2.fun.its1.lefse<- get_df_ridge(fun.its1.clean.ss, imp_tax_fun, sub.mydata.lefse.fun,40)

####Assigning OTU_id
###bacteria
df.ridge_2.bac.lefse <- merge(df.ridge_2.bac.lefse, bac.list.OTU_id, by  = "OTU")

b.manipulate <- sub.mydata.lefse[c('threshold','MeanDecreaseGini','OTU_id')]
ord.mani <- b.manipulate %>% group_by(OTU_id) %>%arrange(MeanDecreaseGini)

b.id <- ord.mani$OTU_id

df.ridge_2.bac.lefse$OTU_id <- factor(df.ridge_2.bac.lefse$OTU_id, levels=b.id)
unique(df.ridge_2.bac.lefse$threshold)

###Fungi
df.ridge_2.fun.lefse <- merge(df.ridge_2.fun.lefse, fun.list.OTU_id, by  = "OTU")

f.manipulate <- sub.mydata.lefse.fun[c('threshold','MeanDecreaseGini','OTU_id')]
f.id <- rev(f.manipulate$OTU_id)
df.ridge_2.fun.lefse$OTU_id <- factor(df.ridge_2.fun.lefse$OTU_id, levels=f.id)

##ITS1
###Fungi
df.ridge_2.fun.its1.lefse <- merge(df.ridge_2.fun.its1.lefse, fun.list.its1.OTU_id, by  = "OTU")

f.manipulate <- sub.mydata.lefse.fun[c('threshold','MeanDecreaseGini','OTU_id')]
f.id <- rev(f.manipulate$OTU_id)
df.ridge_2.fun.its1.lefse$OTU_id <- factor(df.ridge_2.fun.its1.lefse$OTU_id, levels=f.id)

library(ggplot2)
library(ggridges)

plot_ridge <- function(df_ridge){
  ggplot(df_ridge, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=threshold)) + 
    # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
    geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5)+
    xlab('\n Diagnosis')+
    ylab("Relative abundance (%) \n") +
    theme(legend.text=element_text(size=12)) + 
    # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8)))+
    guides(size=FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
    theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
    theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
    theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
    scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
    
    geom_vline(xintercept=30.5, color="slategray3", linetype='dashed',size=1)
  
}

df.ridge_2.fun.lefse$Sample <- factor(df.ridge_2.fun.lefse$Sample, levels = c("F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12",
                                                          "F13","F14","F15","F16","F17","F18","F19","F20","F21","F22","F23",
                                                          "F24","F25","F26","F27","F28","F29","F30","F36","F37","F38","F39",
                                                          "F40","F41","F42","F43", "F44","F45","F46","F47","F48","F49","F50",
                                                          "F51","F52","F53","F54","F55","F56","F57","F58","F59","F60","F61",
                                                          "F62","F63","F64","F65","F71","F72","F73","F74","F75","F76","F77",
                                                          "F78","F79","F81","F82","F83","F84","F87","F88","F89","F90","F91",
                                                          "F93","F94","F96","F97","F98","F99","F100","F101","F102","F103","F104",
                                                          "F105","F106","F107","F108","F109","F110","F111","F112","F113",
                                                          "F114","F115","F116","F117","F118","F119","F120","F121","F122","F123",
                                                          "F124","F125","F126","F128","F129","F130","F131","F132","F133","F134",
                                                          "F135","F136","F137","F138","F139","F140","F141","F142","F143",
                                                          "F144","F145"))
plot_ridge(df.ridge_2.bac.lefse)

plot_ridge(df.ridge_2.fun.lefse)

df.ridge_2.fun.its1.lefse$Sample <- factor(df.ridge_2.fun.its1.lefse$Sample, levels = c("F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12",
                                                                              "F13","F14","F15","F16","F17","F18","F19","F20","F21","F22","F23",
                                                                              "F24","F25","F26","F27","F28","F29","F30","F36","F37","F38","F39",
                                                                              "F40","F41","F42","F43", "F44","F45","F46","F47","F48","F49","F50",
                                                                              "F51","F52","F53","F54","F55","F56","F57","F58","F59","F60","F61",
                                                                              "F62","F63","F64","F65","F71","F72","F73","F74","F75","F76","F77",
                                                                              "F78","F79","F81","F82","F83","F84","F87","F88","F89","F90","F91",
                                                                              "F93","F94","F96","F97","F98","F99","F100","F101","F102","F103","F104",
                                                                              "F105","F106","F107","F108","F109","F110","F111","F112","F113",
                                                                              "F114","F115","F116","F117","F118","F119","F120","F121","F122","F123",
                                                                              "F124","F125","F126","F128","F129","F130","F131","F132","F133","F134",
                                                                              "F135","F136","F137","F138","F139","F140","F141","F142","F143",
                                                                              "F144","F145"))

plot_ridge(df.ridge_2.fun.its1.lefse)


daOTU_bac <- read.xlsx("Bac_daOTUs_and_RF.xlsx",1)

daOTU_bac.id <- merge(daOTU_bac, bac.list.OTU_id, by = "OTU")
write.xlsx(daOTU_bac.id, "Bac_daOTUs_and_RF_with id.xlsx")



daOTU_fun <- read.xlsx("Fun_daOTUs_with RF.xlsx",1)

daOTU_fun.id <- merge(daOTU_fun, fun.list.OTU_id, by = "OTU")
write.xlsx(daOTU_fun.id, "Fun_daOTUs_and_RF_with id.xlsx")


write.csv(df.ridge_2.bac.lefse,"Raw data for ridge plot_random forest with lefse.csv")
write.csv(df.ridge_2.fun.lefse,"Raw data for ridge plot_random forest with lefse_fungi.csv")





#### Ridge plot for commonly detected OTUs (20 bacteria and 25 fungi)
### For bacteria
df.ridge_2.bac.lefse.da<-subset(df.ridge_2.bac.lefse, threshold %in% c("Control","RA"))
b.id <- rev(bac.list$OTU_id)
df.ridge_2.bac.lefse.da$OTU_id <- factor(df.ridge_2.bac.lefse.da$OTU_id, levels=b.id)

plot_ridge(df.ridge_2.bac.lefse.da)+theme(aspect.ratio = 2)



df.ridge_2.fun.lefse.da<-subset(df.ridge_2.fun.lefse, threshold %in% c("Control","RA"))
f.id <- rev(fun.list$OTU_id)
df.ridge_2.fun.lefse.da$OTU_id <- factor(df.ridge_2.fun.lefse.da$OTU_id, levels=f.id)

plot_ridge(df.ridge_2.fun.lefse.da)+theme(aspect.ratio = 2)


write.csv(df.ridge_2.fun.lefse.da,"Fig. 3B.csv")
