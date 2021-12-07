### Control vs nontreated RA vs treated RA
control.subject
non.treated.RA

b.meta$treatment <- 0
b.meta$treatment[which(b.meta$SampleID %in% control.subject)] <- "Control"
b.meta$treatment[which(b.meta$SampleID %in% non.treated.RA)] <- "Not_treated"
b.meta$treatment[which(b.meta$treatment == 0)] <- "Treated"


f.meta <- meta(fun.clean.ss)
control.subject.f<-f.meta$SampleID[which(f.meta$Diagnosis == "control")]
non.treated.RA.f<-f.meta$SampleID[which(f.meta$MTX == "0" & f.meta$TNF_inhibitor == "0")]

f.meta$treatment <- 0
f.meta$treatment[which(f.meta$SampleID %in% control.subject.f)] <- "Control"
f.meta$treatment[which(f.meta$SampleID %in% non.treated.RA.f)] <- "Not_treated"
f.meta$treatment[which(f.meta$treatment == 0)] <- "Treated"


sample_data(fun.clean.ss) <- sample_data(f.meta)
sample_data(bac.clean.ss) <- sample_data(b.meta)


###Treatment methods
b.meta$treatment2 <- 0
b.meta$treatment2[which(b.meta$SampleID %in% control.subject)] <- "Control"
b.meta$treatment2[which(b.meta$SampleID %in% non.treated.RA)] <- "Not_treated"
b.meta$treatment2[which(b.meta$SampleID %in% MTX_subject)] <- "MTX_treated"
b.meta$treatment2[which(b.meta$SampleID %in% TNF_subject)] <- "TNF_treated"
b.meta$treatment2[which(b.meta$treatment2 == 0)] <- "Both"


f.meta <- meta(fun.clean.ss)
control.subject.f<-f.meta$SampleID[which(f.meta$Diagnosis == "control")]
MTX_subject.f<-f.meta$SampleID[which(f.meta$MTX == 1 & f.meta$TNF_inhibitor == 0)]
TNF_subject.f<-f.meta$SampleID[which(f.meta$MTX == 0 & f.meta$TNF_inhibitor !=0 )]

f.meta$treatment2 <- 0
f.meta$treatment2[which(f.meta$SampleID %in% control.subject.f)] <- "Control"
f.meta$treatment2[which(f.meta$SampleID %in% non.treated.RA.f)] <- "Not_treated"
f.meta$treatment2[which(f.meta$SampleID %in% MTX_subject.f)] <- "MTX_treated"
f.meta$treatment2[which(f.meta$SampleID %in% TNF_subject.f)] <- "TNF_treated"
f.meta$treatment2[which(f.meta$treatment2 == 0)] <- "Both"


sample_data(fun.clean.ss) <- sample_data(f.meta)
sample_data(bac.clean.ss) <- sample_data(b.meta)



##Relative abundance
df.genus <- fun.clean.ss %>%
  tax_glom(taxrank = "Genus", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.genus %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))


levels(df.genus$Genus)
levels(df.genus$Genus) = c(levels(df.genus$Genus), 'Low abundance')

# we need to group by samples
df.genus.rel <- df.genus %>%  
  group_by(treatment2) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.genus.rel[df.genus.rel$RelAbundance < 0.5,]$Genus <- 'Low abundance'
unique(df.genus$Genus)

ord <- df.genus.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.genus <- vec.charac[-which(vec.charac %in% c("Low abundance", "unidentified"))]
vec.Low <- c("Low abundance","unidentified")
vec.reorder <- append(vec.Low, vec.genus)


df.genus.rel$Genus <- factor(df.genus.rel$Genus, levels = vec.reorder) 
df.genus.rel$treatment2 <- factor(df.genus.rel$treatment2, levels = c('Control', 'Not_treated', 'MTX_treated', 'TNF_treated', 'Both')) 
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.genus.rel.p1 <- ggplot(df.genus.rel, aes(x=treatment2, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_gen) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=1)

df.genus.rel.p1


write.csv(df.genus.rel, "Fig. 5A.csv")

#### Statistical analysis - Bacteria
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

goodsamps_phy_melt$Genus2 <- ifelse(goodsamps_phy_melt$Genus == "unidentified",ifelse(goodsamps_phy_melt$Family=="unidentified",paste0(goodsamps_phy_melt$Order,'_',"Unidentified genus"),paste0(goodsamps_phy_melt$Family,'_',"Unidentified genus")),paste0(goodsamps_phy_melt$Genus))
goodsamps_phy_melt$Genus2[grep("Prevotella ", goodsamps_phy_melt$Genus2)] <- "Prevotella"
goodsamps_phy_melt$Genus2[grep("Ruminococcus ", goodsamps_phy_melt$Genus2)] <- "Ruminococcus"
goodsamps_phy_melt$Genus2[grep("Coprococcus ", goodsamps_phy_melt$Genus2)] <- "Coprococcus"
goodsamps_phy_melt$Genus2[grep("Ruminiclostridium ", goodsamps_phy_melt$Genus2)] <- "Ruminiclostridium"
goodsamps_phy_melt$Genus2[grep("Tyzzerella ", goodsamps_phy_melt$Genus2)] <- "Tyzzerella"


sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "treatment2","Genus2","Abundance"))
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


sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Genus2 == "Subdoligranulum")
sub_phy_melt_totals.proteo

max.diversity <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$treatment2), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

sub_phy_melt_totals.proteo$treatment2 <- as.factor(sub_phy_melt_totals.proteo$treatment2)
sub_phy_melt_totals.proteo$treatment2 <- factor(sub_phy_melt_totals.proteo$treatment2, levels = c('Control', 'Not_treated',"MTX_treated","TNF_treated", "Both"))


kw<-kruskal.test(PercentAbund ~ treatment2, data = sub_phy_melt_totals.proteo)
kw$p.value

#library(FSA)
DT = dunnTest(PercentAbund ~ treatment2,
              data=sub_phy_melt_totals.proteo,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.unadj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')


p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=treatment2, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+ theme(aspect.ratio=0.6)+
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") +
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




#### Statistical analysis - Fungi
good_phylum <-tax_glom(fun.its1.clean.ss, taxrank = "Genus")
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
sub_phy_melt_totals.proteo$treatment2 <- factor(sub_phy_melt_totals.proteo$treatment2, levels = c('Control', 'Not_treated',"MTX_treated","TNF_treated", "Both"))


kw<-kruskal.test(PercentAbund ~ treatment2, data = sub_phy_melt_totals.proteo)
kw$p.value

#library(FSA)
DT = dunnTest(PercentAbund ~ treatment2,
              data=sub_phy_melt_totals.proteo,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.unadj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')


p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=treatment2, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+ theme(aspect.ratio=0.6)+
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") +
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


write.csv(sub_phy_melt_totals.proteo, "Fig. 5B.csv")

### distribution of duration
both.duration<-b.meta$Ds_Duration[which(b.meta$treatment2 == "Both")]
both.duration <- as.numeric(as.character(both.duration))
mean(both.duration) #12.65


MTX.duration<-b.meta$Ds_Duration[which(b.meta$treatment2 == "MTX_treated")]
MTX.duration <- as.numeric(as.character(MTX.duration))
mean(MTX.duration) #8.051282

TNF.duration<-b.meta$Ds_Duration[which(b.meta$treatment2 == "TNF_treated")]
TNF.duration <- as.numeric(as.character(TNF.duration))
mean(TNF.duration) #8.95


not.duration<-b.meta$Ds_Duration[which(b.meta$treatment2 == "Not_treated")]
not.duration <- as.numeric(as.character(not.duration))
mean(not.duration) # 6.55


##Distribution of active disease
both.active<-b.meta$active_disease[which(b.meta$treatment2 == "Both")]
both.active <- as.numeric(as.character(both.active))
mean(both.active) #0.4


MTX.active<-b.meta$active_disease[which(b.meta$treatment2 == "MTX_treated")]
MTX.active <- as.numeric(as.character(MTX.active))
mean(MTX.active) #1.4

TNF.active<-b.meta$active_disease[which(b.meta$treatment2 == "TNF_treated")]
TNF.active <- as.numeric(as.character(TNF.active))
mean(TNF.active) #0.3


not.active<-b.meta$active_disease[which(b.meta$treatment2 == "Not_treated")]
not.active <- as.numeric(as.character(not.active))
mean(not.active) # 1.3




### Duration
df.genus <- fun.clean.ss %>%
  tax_glom(taxrank = "Genus", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.genus %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))


levels(df.genus$Genus)
levels(df.genus$Genus) = c(levels(df.genus$Genus), 'Low abundance')

# we need to group by samples
df.genus.rel <- df.genus %>%  
  group_by(Duration) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.genus.rel[df.genus.rel$RelAbundance < 0.05,]$Genus <- 'Low abundance'
unique(df.genus$Genus)

ord <- df.genus.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.genus <- vec.charac[-which(vec.charac %in% c("Low abundance", "unidentified"))]
vec.Low <- c("Low abundance","unidentified")
vec.reorder <- append(vec.Low, vec.genus)


df.genus.rel$Genus <- factor(df.genus.rel$Genus, levels = vec.reorder) 
df.genus.rel$Duration <- factor(df.genus.rel$Duration, levels = c('control', 'earlyRA', 'establishedRA')) 
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.genus.rel.p1 <- ggplot(df.genus.rel, aes(x=Duration, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_gen) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 15,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=1)

df.genus.rel.p1

