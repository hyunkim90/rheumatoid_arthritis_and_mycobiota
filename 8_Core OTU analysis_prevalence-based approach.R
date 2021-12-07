# check the data 
##Relative abundance


############
##Read OTU tables and convert to phyloseq object
bac.clean.ss


##phyloseq files
bac.clean.ss

fun.clean.ss


# Calculate compositional version of the data
# (relative abundances)
bac.clean.ss.rel <- microbiome::transform(bac.clean.ss, "compositional")
t(otu_table(bac.clean.ss.rel))

fun.clean.ss.rel <- microbiome::transform(fun.clean.ss, "compositional")
t(otu_table(fun.clean.ss.rel))

# Filter the data to include only healthy subjects
#pseq.1 <- subset_samples(phy.clean.l, ibd_subtype == "HC" & timepoint == "1") 
#print(pseq.1)
# keep only taxa with positive sums
#pseq.2 <- prune_taxa(taxa_sums(pseq.1) > 0, pseq.1)

#print(pseq.2)

#Relative population frequencies; at 1% compositional abundance threshold:
  # 
  # head(prevalence(bac.clean.ss.rel, detection = 0.01, sort = TRUE))
  # 
  # head(prevalence(fun.clean.ss.rel, detection = 0.01, sort = TRUE))
  # 
  
#We can see that only OTU ids are listed with no taxonomic information. Absolute population frequencies (sample count):
  

#Core microbiota analysis
#If you only need the names of the core taxa, do as follows. This returns the taxa that exceed the given prevalence and detection thresholds.

  
core.bac.50 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 50/100)
core.bac.50

core.bac.60 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 60/100)
core.bac.60  

core.bac.70 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 70/100)
core.bac.70
  
core.bac.80 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 80/100)
core.bac.80

core.bac.85 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 85/100)
core.bac.85

core.bac.90 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 90/100)
core.bac.90

core.bac.95 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 95/100)
core.bac.95

tax_core_bac <- subset(bac.list, OTU%in% core.bac.85)

write.csv(tax_core_bac,"Core OTUs-85% prevalence.csv")

tax_core_bac <- subset(bac.list, OTU%in% core.bac.70)

core.fun.95 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 95/100)
core.fun.95

core.fun.90 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 90/100)
core.fun.90

core.fun.85 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 85/100)
core.fun.85

core.fun.80 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 80/100)
core.fun.80

core.fun.70 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 70/100)
core.fun.70


tax_core_fun <- subset(fun.list, OTU%in% core.fun.70)

tax_core_fun$OTU_id
###fungal ITS1 core
fun.its1.clean.ss.rel <- microbiome::transform(fun.its1.clean.ss, "compositional")
core.fun.its1.95 <- core_members(fun.its1.clean.ss.rel, detection = 0, prevalence = 95/100)
core.fun.its1.95

core.fun.its1.90 <- core_members(fun.its1.clean.ss.rel, detection = 0, prevalence = 90/100)
core.fun.its1.90

core.fun.its1.85 <- core_members(fun.its1.clean.ss.rel, detection = 0, prevalence = 85/100)
core.fun.its1.85

core.fun.its1.80 <- core_members(fun.its1.clean.ss.rel, detection = 0, prevalence = 80/100)
core.fun.its1.80

core.fun.its1.70 <- core_members(fun.its1.clean.ss.rel, detection = 0, prevalence = 70/100)
core.fun.its1.70

tax_core_fun.its1 <- subset(fun.list.its1, OTU%in% core.fun.its1.70)

##RA
bac.clean.ss.RA <- subset_samples(bac.clean.ss, Diagnosis == "RA") 
bac.clean.ss.RA <- phyloseq::filter_taxa(bac.clean.ss.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs

bac.clean.ss.RA.rel <- microbiome::transform(bac.clean.ss.RA, "compositional")
t(otu_table(bac.clean.ss.RA.rel))

core.bac.40.ra <- core_members(bac.clean.ss.RA.rel, detection = 0, prevalence = 40/100)
core.bac.40.ra

core.bac.50.ra <- core_members(bac.clean.ss.RA.rel, detection = 0, prevalence = 50/100)
core.bac.50.ra

core.bac.60.ra <- core_members(bac.clean.ss.RA.rel, detection = 0, prevalence = 60/100)
core.bac.60.ra

core.bac.70.ra <- core_members(bac.clean.ss.RA.rel, detection = 0, prevalence = 70/100)
core.bac.70.ra

core.bac.80.ra <- core_members(bac.clean.ss.RA.rel, detection = 0, prevalence = 80/100)
core.bac.80.ra

core.bac.85.ra <- core_members(bac.clean.ss.RA.rel, detection = 0, prevalence = 85/100)
core.bac.85.ra

core.bac.90.ra <- core_members(bac.clean.ss.RA.rel, detection = 0, prevalence = 90/100)
core.bac.90.ra

core.bac.95.ra <- core_members(bac.clean.ss.RA.rel, detection = 0, prevalence = 95/100)
core.bac.95.ra

tax_core_ra_bac <- subset(bac.list, OTU%in% core.bac.60.ra)

tax_bac.core.ra1 <- subset(tax_bac, rownames(tax_bac) %in% core.bac.50.ra)
tax_bac.core.ra2 <- subset(tax_bac, rownames(tax_bac) %in% core.bac.70.ra)
tax_bac.core.ra3 <- subset(tax_bac, rownames(tax_bac) %in% core.bac.40.ra)
write.table(tax_bac.core.ra1, "50% core-RA.tsv", sep='\t', quote=F)

intersect(core.bac.50.ra, core.bac.50)

otu_bac.rel<-otu_table(bac.clean.ss.rel)
otu_bac.rel.core <- subset(otu_bac.rel, rownames(otu_bac.rel)%in%core.bac.85.ra)
write.table(otu_bac.rel.core, "abundance table of 85% core-RA.tsv", sep='\t', quote=F)

otu_bac.normal<-otu_table(bac.clean.nolog)
otu_bac.normal.core <- subset(otu_bac.normal, rownames(otu_bac.normal)%in%core.bac.85.ra)
write.table(otu_bac.normal.core, "normalized abundance table of 85% core-RA.tsv", sep='\t', quote=F)


###Fungi
fun.clean.ss.RA

fun.clean.ss.RA.rel <- microbiome::transform(fun.clean.ss.RA, "compositional")
t(otu_table(fun.clean.ss.RA.rel))

core.fun.40.ra <- core_members(fun.clean.ss.RA.rel, detection = 0, prevalence = 40/100)
core.fun.40.ra

core.fun.50.ra <- core_members(fun.clean.ss.RA.rel, detection = 0, prevalence = 50/100)
core.fun.50.ra

core.fun.60.ra <- core_members(fun.clean.ss.RA.rel, detection = 0, prevalence = 60/100)
core.fun.60.ra

core.fun.70.ra <- core_members(fun.clean.ss.RA.rel, detection = 0, prevalence = 70/100)
core.fun.70.ra # "SH181823.07FU_FR819718_refs" "SH216041.07FU_AF138287_refs" "SH211137.07FU_GU721261_reps" "SH194551.07FU_GQ280331_reps"

core.fun.80.ra <- core_members(fun.clean.ss.RA.rel, detection = 0, prevalence = 80/100)
core.fun.80.ra #"SH194551.07FU_GQ280331_reps"

core.fun.85.ra <- core_members(fun.clean.ss.RA.rel, detection = 0, prevalence = 85/100)
core.fun.85.ra 

core.fun.90.ra <- core_members(fun.clean.ss.RA.rel, detection = 0, prevalence = 90/100)
core.fun.90.ra

core.fun.95.ra <- core_members(fun.clean.ss.RA.rel, detection = 0, prevalence = 95/100)
core.fun.95.ra

### fungal ITS1 -RA
fun.its1.clean.ss.RA

fun.its1.clean.ss.RA.rel <- microbiome::transform(fun.its1.clean.ss.RA, "compositional")
t(otu_table(fun.its1.clean.ss.RA.rel))

core.fun.its1.40.ra <- core_members(fun.its1.clean.ss.RA.rel, detection = 0, prevalence = 40/100)
core.fun.its1.40.ra

core.fun.its1.50.ra <- core_members(fun.its1.clean.ss.RA.rel, detection = 0, prevalence = 50/100)
core.fun.its1.50.ra

core.fun.its1.60.ra <- core_members(fun.its1.clean.ss.RA.rel, detection = 0, prevalence = 60/100)
core.fun.its1.60.ra

core.fun.its1.70.ra <- core_members(fun.its1.clean.ss.RA.rel, detection = 0, prevalence = 70/100)
core.fun.its1.70.ra # "SH181823.07FU_FR819718_refs" "SH216041.07FU_AF138287_refs" "SH211137.07FU_GU721261_reps" "SH194551.07FU_GQ280331_reps"

core.fun.its1.80.ra <- core_members(fun.its1.clean.ss.RA.rel, detection = 0, prevalence = 80/100)
core.fun.its1.80.ra #"SH194551.07FU_GQ280331_reps"

core.fun.its1.85.ra <- core_members(fun.its1.clean.ss.RA.rel, detection = 0, prevalence = 85/100)
core.fun.its1.85.ra 

core.fun.its1.90.ra <- core_members(fun.its1.clean.ss.RA.rel, detection = 0, prevalence = 90/100)
core.fun.its1.90.ra

core.fun.its1.95.ra <- core_members(fun.its1.clean.ss.RA.rel, detection = 0, prevalence = 95/100)
core.fun.its1.95.ra



### All samples
fun.clean.ss.RA

fun.clean.ss.rel <- microbiome::transform(fun.clean.ss, "compositional")

core.fun.40 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 40/100)
core.fun.40

core.fun.50 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 50/100)
core.fun.50

core.fun.60 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 60/100)
core.fun.60

core.fun.70 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 70/100)
core.fun.70 # "SH181823.07FU_FR819718_refs" "SH216041.07FU_AF138287_refs" "SH211137.07FU_GU721261_reps" "SH194551.07FU_GQ280331_reps"

core.fun.80 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 80/100)
core.fun.80 #"SH194551.07FU_GQ280331_reps"

core.fun.85 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 85/100)
core.fun.85 

core.fun.90 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 90/100)
core.fun.90

core.fun.95 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 95/100)
core.fun.95


tax_core_ra_fun <- subset(fun.list, OTU%in% core.fun.70.ra)
tax_core_fun <- subset(fun.list, OTU%in% core.fun.70)
write.table(tax_core_fun, "fungal 70% core.tsv", sep='\t', quote=F)

otu_fun.rel<-otu_table(fun.clean.ss.rel)
otu_fun.rel.core <- subset(otu_fun.rel, rownames(otu_fun.rel)%in%core.fun.70)
write.table(otu_fun.rel.core, "fungal abundance table of 70% core-RA.tsv", sep='\t', quote=F)

otu_fun.normal<-otu_table(fun.clean.nolog)
otu_fun.normal.core <- subset(otu_fun.normal, rownames(otu_fun.normal)%in%core.fun.70.ra)
write.table(otu_fun.normal.core, "fungal normalized abundance table of 70% core-RA.tsv", sep='\t', quote=F)

## Core's correlation with factors
### Abundance table

otu_bac.rel.core.t <- t(otu_bac.rel.core)
b.meta.abun <- merge(b.meta, otu_bac.rel.core.t, by = "row.names")

otu_fun.rel.core.t <- t(otu_fun.rel.core)
f.meta.abun <- merge(f.meta, otu_fun.rel.core.t, by = "row.names")

## Correlation between OTU abundance and each quantitative factors
cor_5.otu <- Hmisc::rcorr(b.meta.abun$BMI, b.meta.abun$DQ802087.1.1257, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$Age, b.meta.abun$DQ802087.1.1257, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$anti_CCP, b.meta.abun$DQ802087.1.1257, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$RA_factor, b.meta.abun$DQ802087.1.1257, type="spearman")
cor_5.otu$r
cor_5.otu$P


cor_5.otu <- Hmisc::rcorr(b.meta.abun$BMI, b.meta.abun$DQ808255.1.1393, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$Age, b.meta.abun$DQ808255.1.1393, type="spearman")
cor_5.otu$r #-0.1767938
cor_5.otu$P #0.04504181

cor_5.otu <- Hmisc::rcorr(b.meta.abun$anti_CCP, b.meta.abun$DQ808255.1.1393, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$RA_factor, b.meta.abun$DQ808255.1.1393, type="spearman")
cor_5.otu$r
cor_5.otu$P


cor_5.otu <- Hmisc::rcorr(b.meta.abun$BMI, b.meta.abun$CZBP01000074.1.1200, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$Age, b.meta.abun$CZBP01000074.1.1200, type="spearman")
cor_5.otu$r 
cor_5.otu$P 

cor_5.otu <- Hmisc::rcorr(b.meta.abun$anti_CCP, b.meta.abun$CZBP01000074.1.1200, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$RA_factor, b.meta.abun$CZBP01000074.1.1200, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$Total_Cholesterol, b.meta.abun$CZBP01000074.1.1200, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$BMI, b.meta.abun$DQ799508.1.1380, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$Age, b.meta.abun$DQ799508.1.1380, type="spearman")
cor_5.otu$r 
cor_5.otu$P 

cor_5.otu <- Hmisc::rcorr(b.meta.abun$anti_CCP, b.meta.abun$DQ799508.1.1380, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$RA_factor, b.meta.abun$DQ799508.1.1380, type="spearman")
cor_5.otu$r
cor_5.otu$P


cor_5.otu <- Hmisc::rcorr(b.meta.abun$BMI, b.meta.abun$DQ823861.1.1313, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$Age, b.meta.abun$DQ823861.1.1313, type="spearman")
cor_5.otu$r 
cor_5.otu$P 

cor_5.otu <- Hmisc::rcorr(b.meta.abun$anti_CCP, b.meta.abun$DQ823861.1.1313, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(b.meta.abun$RA_factor, b.meta.abun$DQ823861.1.1313, type="spearman")
cor_5.otu$r
cor_5.otu$P

## Fungi
cor_5.otu <- Hmisc::rcorr(f.meta.abun$anti_CCP, f.meta.abun$SH181823.07FU_FR819718_refs, type="spearman")
cor_5.otu$r
cor_5.otu$P


cor_5.otu <- Hmisc::rcorr(f.meta.abun$anti_CCP, f.meta.abun$SH216041.07FU_AF138287_refs, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(f.meta.abun$anti_CCP, f.meta.abun$SH211137.07FU_GU721261_reps, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(f.meta.abun$anti_CCP, f.meta.abun$SH194551.07FU_GQ280331_reps, type="spearman")
cor_5.otu$r
cor_5.otu$P


cor_5.otu <- Hmisc::rcorr(f.meta.abun$RA_factor, f.meta.abun$SH181823.07FU_FR819718_refs, type="spearman")
cor_5.otu$r
cor_5.otu$P


cor_5.otu <- Hmisc::rcorr(f.meta.abun$RA_factor, f.meta.abun$SH216041.07FU_AF138287_refs, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(f.meta.abun$RA_factor, f.meta.abun$SH211137.07FU_GU721261_reps, type="spearman")
cor_5.otu$r
cor_5.otu$P

cor_5.otu <- Hmisc::rcorr(f.meta.abun$RA_factor, f.meta.abun$SH194551.07FU_GQ280331_reps, type="spearman")
cor_5.otu$r
cor_5.otu$P

## Comparison in group
max.diversity <- aggregate(b.meta.abun$CZBP01000074.1.1200, by = list(b.meta.abun$Duration), max)

colnames(max.diversity) <- c("Duration", "MaxDiversity")


b.meta.abun$Duration <- as.factor(b.meta.abun$Duration)
##Kruskal-Wallis test
kw<-kruskal.test(DQ799508.1.1380 ~ Duration, data = b.meta.abun)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(DQ823861.1.1313 ~ Duration,
              data=b.meta.abun,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
names(dunn)[1] <- "Duration"
hsd1 <- merge(dunn,max.diversity, by = 'Duration')
p <- ggplot(data = b.meta.abun, aes(x=Duration, y=DQ799508.1.1380, fill = Duration)) + geom_boxplot(width = 0.8) +
  theme_bw()+#geom_text(data=hsd1,aes(x=Duration,y=log(MaxDiversity), label=Letter), vjust=-0.5,size=5) +
  geom_point(position='jitter',shape=1, alpha=.5)+
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("RA factor\n") + theme(aspect.ratio=1.5)+
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

dev.off()

### 

##Fungi
core.fun.99 <- core_members(fun.clean.ss.field.rel, detection = 0, prevalence = 99/100)
core.fun.99

core.fun <- subset(OTU_id.list, OTU_id.list$OTU%in% core.fun.99)



### Correlation between core OTU abundance and chemical properties
#2017 full data
df.otu.bac.clean.ss
df.otu.arch.clean.ss
df.fun.clean.ss

df.otu.core.bac <- subset(df.otu.bac.clean.ss, rownames(df.otu.bac.clean.ss)%in% core.bac$OTU_id)
df.otu.core.arch <- subset(df.otu.arch.clean.ss, rownames(df.otu.arch.clean.ss)%in% core.arch$OTU_id)
df.otu.core.fun <- subset(df.fun.clean.ss, rownames(df.fun.clean.ss)%in% core.fun$OTU_id)


b.meta.all
f.meta.all

cor_table.core.bac<-cbind(b.meta[,7:18], t(df.otu.core.bac))
summary(cor_table.core.bac)

cor_table.core.arch<-cbind(b.meta[,7:18], t(df.otu.core.arch))
summary(cor_table.core.arch)

cor_table.core.fun<-cbind(f.meta[,7:18], t(df.otu.core.fun))
summary(cor_table.core.fun)


##Correlation
mat.cor_table.core.bac <- as.matrix(cor_table.core.bac)
mat.cor_table.core.arch <- as.matrix(cor_table.core.arch)
mat.cor_table.core.fun <- as.matrix(cor_table.core.fun)


library(Hmisc)

cor_5 <- rcorr(mat.cor_table.core.bac, type="spearman")
M <- cor_5$r
p_mat <- cor_5$P

library(xlsx)
write.xlsx(M, "Bac_core_correlation_rho.xlsx")
write.xlsx(p_mat, "Bac_core_correlation_p value.xlsx")

corrplot::corrplot(M, type = "upper", 
                   p.mat = p_mat, sig.level = 0.05, insig = "blank", tl.col="black", col=brewer.pal(n=10, name="RdYlBu"))

corrplot::corrplot(M, type = "upper", 
                   tl.col="black", col=brewer.pal(n=10, name="RdYlBu"))

dev.off()