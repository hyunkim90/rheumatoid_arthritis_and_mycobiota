########### Alpha diversity ###################
library(microbiome)
library(knitr)
##rarefy
bac.rarefy <- rarefy_even_depth(bac.clean.ss, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

bac.clean.ss.sub <- subset_samples(bac.clean.ss, MTX == "Control"|MTX=="0"& TNF_inhibitor == "0")
bac.clean.ss.sub <- phyloseq::filter_taxa(bac.clean.ss.sub, function(x) sum(x) != 0, TRUE)
bac.rarefy.sub <-rarefy_even_depth(bac.clean.ss.sub, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

fun.rarefy <- rarefy_even_depth(fun.clean.ss, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

fun.clean.ss.sub <- subset_samples(fun.clean.ss, MTX == "Control"|MTX=="0"& TNF_inhibitor == "0")
fun.clean.ss.sub <- phyloseq::filter_taxa(fun.clean.ss.sub, function(x) sum(x) != 0, TRUE)
fun.rarefy.sub <-rarefy_even_depth(fun.clean.ss.sub, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

sample_data(bac.rarefy) <- sample_data(b.meta)

sample_data(fun.rarefy) <- sample_data(f.meta)
########## Baterial community ###############
tab_all <- microbiome::alpha(bac.rarefy, index = "all")
write.table(tab_all, "Alpha diversity_bac.txt", sep = "\t", row.names = TRUE,  quote = TRUE, na = "NA")

ps1.meta <- b.meta

## using Observed OTUs
ps1.meta$Observed <- tab_all$observed 
ps1.meta
ps1.meta$simp <- tab_all$evenness_simpson 
ps1.meta$simp
ps1.meta$Shannon <- tab_all$diversity_shannon 
ps1.meta$Shannon

########Diagnosis

##Richness (Observed OTUs)
max.diversity <- aggregate(ps1.meta$Observed, by = list(ps1.meta$Diagnosis), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

histogram(subset(ps1.meta, Diagnosis=='control')$Observed)
histogram(subset(ps1.meta, Diagnosis=='RA')$Observed)

# wilcoxon test
x <- subset(ps1.meta, Diagnosis=='control')$Observed
y <- subset(ps1.meta, Diagnosis=='RA')$Observed
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.9068

p <- ggplot(data = ps1.meta, aes(x=Diagnosis, y=Observed, fill = Diagnosis)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Observed OTU \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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



##Shannon



max.diversity <- aggregate(ps1.meta$Shannon, by = list(ps1.meta$Diagnosis), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

histogram(subset(ps1.meta, Diagnosis=='control')$Shannon)
histogram(subset(ps1.meta, Diagnosis=='RA')$Shannon)

# wilcoxon test
x <- subset(ps1.meta, Diagnosis=='control')$Shannon
y <- subset(ps1.meta, Diagnosis=='RA')$Shannon
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.8826

p <- ggplot(data = ps1.meta, aes(x=Diagnosis, y=Shannon, fill = Diagnosis)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Shannon \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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



##Simpson evenness


max.diversity <- aggregate(ps1.meta$simp, by = list(ps1.meta$Diagnosis), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

histogram(subset(ps1.meta, Diagnosis=='control')$simp)
histogram(subset(ps1.meta, Diagnosis=='RA')$simp)

# wilcoxon test
x <- subset(ps1.meta, Diagnosis=='control')$simp
y <- subset(ps1.meta, Diagnosis=='RA')$simp
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.8826


p <- ggplot(data = ps1.meta, aes(x=Diagnosis, y=simp, fill = Diagnosis)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Simpson \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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


#########Active disease
max.diversity <- aggregate(ps1.meta$Observed, by = list(ps1.meta$active_disease), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Observed ~ active_disease, data = ps1.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Observed ~ active_disease,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')

p <- ggplot(data = ps1.meta, aes(x=active_disease, y=Observed, fill = active_disease)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Observed OTU\n") + theme(aspect.ratio=0.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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

## using Shannon

##ANOVA
max.diversity <- aggregate(ps1.meta$Shannon, by = list(ps1.meta$active_disease), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Shannon ~ active_disease, data = ps1.meta)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(Shannon ~ active_disease,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=active_disease, y=Shannon, fill = active_disease)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Shannon\n") + theme(aspect.ratio=0.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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


## using evenness Simpson
ps1.meta$simp <- tab_all$evenness_simpson 
ps1.meta$simp

##ANOVA
max.diversity <- aggregate(ps1.meta$simp, by = list(ps1.meta$active_disease), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(simp ~ active_disease, data = ps1.meta)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(simp ~ active_disease,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=active_disease, y=simp, fill = active_disease)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Simpson\n") + theme(aspect.ratio=0.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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



##Duration
ps1.meta$Duration <-as.factor(ps1.meta$Duration)

max.diversity <- aggregate(ps1.meta$Observed, by = list(ps1.meta$Duration), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Observed ~ Duration, data = ps1.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Observed ~ Duration,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')

p <- ggplot(data = ps1.meta, aes(x=Duration, y=Observed, fill = Duration)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Observed OTU\n") + theme(aspect.ratio=1.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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

## using Shannon
max.diversity <- aggregate(ps1.meta$Shannon, by = list(ps1.meta$Duration), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Shannon ~ Duration, data = ps1.meta)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(Shannon ~ Duration,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=Duration, y=Shannon, fill = Duration)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Shannon\n") + theme(aspect.ratio=1.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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


## using evenness Simpson
max.diversity <- aggregate(ps1.meta$simp, by = list(ps1.meta$Duration), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(simp ~ Duration, data = ps1.meta)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(simp ~ Duration,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=Duration, y=simp, fill = Duration)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Simpson\n") + theme(aspect.ratio=1.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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

###### treatment2

ps1.meta$treatment2 <-as.factor(ps1.meta$treatment2)
ps1.meta$treatment2 <-factor(ps1.meta$treatment2, levels = c('Control','Not_treated',"MTX_treated", 'TNF_treated', 'Both'))


max.diversity <- aggregate(ps1.meta$Observed, by = list(ps1.meta$treatment2), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Observed ~ treatment2, data = ps1.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Observed ~ treatment2,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')

p <- ggplot(data = ps1.meta, aes(x=treatment2, y=Observed, fill = treatment2)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Observed OTU\n") + theme(aspect.ratio=0.6)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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

## using Shannon
max.diversity <- aggregate(ps1.meta$Shannon, by = list(ps1.meta$treatment2), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Shannon ~ treatment2, data = ps1.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Shannon ~ treatment2,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=treatment2, y=Shannon, fill = treatment2)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Shannon\n") + theme(aspect.ratio=0.6)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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


## using evenness Simpson
max.diversity <- aggregate(ps1.meta$simp, by = list(ps1.meta$treatment2), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(simp ~ treatment2, data = ps1.meta)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(simp ~ treatment2,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=treatment2, y=simp, fill = treatment2)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Simpson\n") + theme(aspect.ratio=0.6)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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




########## fungal community ###############
tab_all.f <- microbiome::alpha(fun.rarefy, index = "all")
write.table(tab_all.f, "Alpha diversity_fun.txt", sep = "\t", row.names = TRUE,  quote = TRUE, na = "NA")

ps2.meta <- f.meta

## using Observed OTUs
ps2.meta$Observed <- tab_all.f$observed 
ps2.meta


##Diagnosis

##Richness (Observed OTUs)
max.diversity <- aggregate(ps2.meta$Observed, by = list(ps2.meta$Diagnosis), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

histogram(subset(ps2.meta, Diagnosis=='control')$Observed)
histogram(subset(ps2.meta, Diagnosis=='RA')$Observed)

# wilcoxon test
x <- subset(ps2.meta, Diagnosis=='control')$Observed
y <- subset(ps2.meta, Diagnosis=='RA')$Observed
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.7316

p <- ggplot(data = ps2.meta, aes(x=Diagnosis, y=Observed, fill = Diagnosis)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Observed OTU \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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



##Shannon
ps2.meta$Shannon <- tab_all.f$diversity_shannon 
ps2.meta$Shannon


max.diversity <- aggregate(ps2.meta$Shannon, by = list(ps2.meta$Diagnosis), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

histogram(subset(ps2.meta, Diagnosis=='control')$Shannon)
histogram(subset(ps2.meta, Diagnosis=='RA')$Shannon)

# wilcoxon test
x <- subset(ps2.meta, Diagnosis=='control')$Shannon
y <- subset(ps2.meta, Diagnosis=='RA')$Shannon
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.8475

p <- ggplot(data = ps2.meta, aes(x=Diagnosis, y=Shannon, fill = Diagnosis)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Shannon \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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



##Simpson evenness
ps2.meta$simp <- tab_all.f$evenness_simpson 
ps2.meta$simp

max.diversity <- aggregate(ps2.meta$simp, by = list(ps2.meta$Diagnosis), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

histogram(subset(ps2.meta, Diagnosis=='control')$simp)
histogram(subset(ps2.meta, Diagnosis=='RA')$simp)

# wilcoxon test
x <- subset(ps2.meta, Diagnosis=='control')$simp
y <- subset(ps2.meta, Diagnosis=='RA')$simp
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.5944


p <- ggplot(data = ps2.meta, aes(x=Diagnosis, y=simp, fill = Diagnosis)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Simpson \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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


##Active disease
max.diversity <- aggregate(ps2.meta$Observed, by = list(ps2.meta$active_disease), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Observed ~ active_disease, data = ps2.meta)
kw$p.value #0.639939

#library(FSA)
DT = dunnTest(Observed ~ active_disease,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')

p <- ggplot(data = ps2.meta, aes(x=active_disease, y=Observed, fill = active_disease)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Observed OTU\n") + theme(aspect.ratio=0.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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

## using Shannon

max.diversity <- aggregate(ps2.meta$Shannon, by = list(ps2.meta$active_disease), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Shannon ~ active_disease, data = ps2.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Shannon ~ active_disease,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=active_disease, y=Shannon, fill = active_disease)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Shannon\n") + theme(aspect.ratio=0.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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


## using evenness Simpson
max.diversity <- aggregate(ps2.meta$simp, by = list(ps2.meta$active_disease), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(simp ~ active_disease, data = ps2.meta)
kw$p.value # 0.1181907
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(simp ~ active_disease,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=active_disease, y=simp, fill = active_disease)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Simpson\n") + theme(aspect.ratio=0.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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

##Active disease 2
ps2.meta$active_disease_2 <-as.factor(ps2.meta$active_disease_2)

max.diversity <- aggregate(ps2.meta$Observed, by = list(ps2.meta$active_disease_2), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Observed ~ active_disease_2, data = ps2.meta)
kw$p.value #0.7278378

#library(FSA)
DT = dunnTest(Observed ~ active_disease_2,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')

p <- ggplot(data = ps2.meta, aes(x=active_disease_2, y=Observed, fill = active_disease_2)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Observed OTU\n") + theme(aspect.ratio=1.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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

## using Shannon

max.diversity <- aggregate(ps2.meta$Shannon, by = list(ps2.meta$active_disease_2), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Shannon ~ active_disease_2, data = ps2.meta)
kw$p.value #0.9553377

#library(FSA)
DT = dunnTest(Shannon ~ active_disease_2,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=active_disease_2, y=Shannon, fill = active_disease_2)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Shannon\n") + theme(aspect.ratio=1.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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


## using evenness Simpson
max.diversity <- aggregate(ps2.meta$simp, by = list(ps2.meta$active_disease_2), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(simp ~ active_disease_2, data = ps2.meta)
kw$p.value # 0.8586115
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(simp ~ active_disease_2,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=active_disease_2, y=simp, fill = active_disease_2)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Simpson\n") + theme(aspect.ratio=1.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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



##Duration
ps2.meta$Duration <-as.factor(ps2.meta$Duration)

max.diversity <- aggregate(ps2.meta$Observed, by = list(ps2.meta$Duration), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Observed ~ Duration, data = ps2.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Observed ~ Duration,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')

p <- ggplot(data = ps2.meta, aes(x=Duration, y=Observed, fill = Duration)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Observed OTU\n") + theme(aspect.ratio=1.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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

## using Shannon

##ANOVA
max.diversity <- aggregate(ps2.meta$Shannon, by = list(ps2.meta$Duration), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Shannon ~ Duration, data = ps2.meta)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(Shannon ~ Duration,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=Duration, y=Shannon, fill = Duration)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Shannon\n") + theme(aspect.ratio=1.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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


## using evenness Simpson
max.diversity <- aggregate(ps2.meta$simp, by = list(ps2.meta$Duration), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(simp ~ Duration, data = ps2.meta)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(simp ~ Duration,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=Duration, y=simp, fill = Duration)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Simpson\n") + theme(aspect.ratio=1.5)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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


##MTX

##Richness (Observed OTUs)
max.diversity <- aggregate(ps2.meta$Observed, by = list(ps2.meta$MTX), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

ps2.meta$MTX <-as.factor(ps2.meta$MTX)


kw<-kruskal.test(Observed ~ MTX, data = ps2.meta)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(simp ~ MTX,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')

p <- ggplot(data = ps2.meta, aes(x=MTX, y=Observed, fill = MTX)) + geom_boxplot(width = 0.8) +
  theme_bw() + #scale_fill_manual(values=c("0" = "#6699CC", "1" = "#CC9900"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Observed OTU \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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



##Shannon
max.diversity <- aggregate(ps2.meta$Shannon, by = list(ps2.meta$MTX), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

kw<-kruskal.test(Shannon ~ MTX, data = ps2.meta)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(Shannon ~ MTX,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')


p <- ggplot(data = ps2.meta, aes(x=MTX, y=Shannon, fill = MTX)) + geom_boxplot(width = 0.8) +
  theme_bw() + #scale_fill_manual(values=c("0" = "#6699CC", "1" = "#CC9900"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Shannon \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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



##Simpson evenness
max.diversity <- aggregate(ps2.meta$simp, by = list(ps2.meta$MTX), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

histogram(subset(ps2.meta, MTX=='0')$simp)
histogram(subset(ps2.meta, MTX=='1')$simp)

# wilcoxon test
x <- subset(ps2.meta, MTX=='0')$simp
y <- subset(ps2.meta, MTX=='1')$simp
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.8674


p <- ggplot(data = ps2.meta, aes(x=MTX, y=simp, fill = MTX)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("0" = "#6699CC", "1" = "#CC9900"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Simpson \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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




### Regression analysis 
lmHeight = lm(Age~Shannon, data = ps2.meta)
summary(lmHeight) #Adjusted R-squared:  -0.002192  p-value: 0.3977

lmHeight = lm(Age~Observed, data = ps2.meta) 
summary(lmHeight) #Adjusted R-squared:  -0.007807 p-value: 0.9272

lmHeight = lm(Age~simp, data = ps2.meta) 
summary(lmHeight) #Adjusted R-squared:  0.01545 p-value: 0.08525


lmHeight = lm(BMI~Shannon, data = ps2.meta)
summary(lmHeight) #Adjusted R-squared:  -0.006594  p-value: 0.6884

lmHeight = lm(BMI~Observed, data = ps2.meta) 
summary(lmHeight) #Adjusted R-squared:  -0.005868 p-value: 0.6157

lmHeight = lm(BMI~simp, data = ps2.meta) 
summary(lmHeight) #Adjusted R-squared:  0.006212 p-value: 0.1821



lmHeight = lm(RA_factor~Shannon, data = ps2.meta)
summary(lmHeight) #Adjusted R-squared:  0.0006358  p-value: 0.3004

lmHeight = lm(RA_factor~Observed, data = ps2.meta) 
summary(lmHeight) #Adjusted R-squared:  -0.005645 p-value: 0.5967

lmHeight = lm(RA_factor~simp, data = ps2.meta) 
summary(lmHeight) #Adjusted R-squared:  -0.005635  p-value: 0.5958


lmHeight = lm(anti_CCP~Shannon, data = ps2.meta)
summary(lmHeight) #Adjusted R-squared:  0.03332   p-value: 0.02158

lmHeight = lm(anti_CCP~Observed, data = ps2.meta) 
summary(lmHeight) #Adjusted R-squared:  -0.004951 p-value: 0.5444

lmHeight = lm(anti_CCP~simp, data = ps2.meta) 
summary(lmHeight) #Adjusted R-squared:  0.008995   p-value: 0.1439


ggplot(ps2.meta, aes(x=Age, y=Shannon)) +
  xlab('\n Age')+
  ylab("Shannon \n") +
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

ggplot(ps2.meta, aes(x=Age, y=Observed)) +
  xlab('\n Age')+
  ylab("Shannon \n") +
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

ggplot(ps2.meta, aes(x=Age, y=simp)) +
  xlab('\n Age')+
  ylab("Shannon \n") +
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

ggplot(ps2.meta, aes(x=anti_CCP, y=Shannon)) +
  xlab('\n Age')+
  ylab("Shannon \n") +
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




###### treatment2

ps2.meta$treatment2 <-as.factor(ps2.meta$treatment2)
ps2.meta$treatment2 <-factor(ps2.meta$treatment2, levels = c('Control','Not_treated',"MTX_treated", 'TNF_treated', 'Both'))

max.diversity <- aggregate(ps2.meta$Observed, by = list(ps2.meta$treatment2), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Observed ~ treatment2, data = ps2.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Observed ~ treatment2,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')

p <- ggplot(data = ps2.meta, aes(x=treatment2, y=Observed, fill = treatment2)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Observed OTU\n") + theme(aspect.ratio=0.6)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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

## using Shannon
max.diversity <- aggregate(ps2.meta$Shannon, by = list(ps2.meta$treatment2), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(Shannon ~ treatment2, data = ps2.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Shannon ~ treatment2,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=treatment2, y=Shannon, fill = treatment2)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Shannon\n") + theme(aspect.ratio=0.6)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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


## using evenness Simpson
max.diversity <- aggregate(ps2.meta$simp, by = list(ps2.meta$treatment2), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(simp ~ treatment2, data = ps2.meta)
kw$p.value

#library(FSA)
DT = dunnTest(simp ~ treatment2,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=treatment2, y=simp, fill = treatment2)) + geom_boxplot(width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Simpson\n") + theme(aspect.ratio=0.6)+
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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


