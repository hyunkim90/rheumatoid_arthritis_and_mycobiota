### 

# Set colors for plotting

### Bacteria

#### get unassigned vectors
#### get CP and MT phyloseq obj and vectors
# (1) CP
phy.cp <- subset_taxa(phy, Class == "D_2__Chloroplast") ## just confirming
phy.cp <- subset_taxa(phy, Order == "D_3__Chloroplast")
vec.cp <- rownames(otu_table(phy.cp))
length(rownames(otu_table(phy.cp))) ## 130 otus of CP
vec.cp

# (2) MT
# phy.mt <- subset_taxa(phy, Family == "D_4__Mitochondria")
# phy.mt <- subset_taxa(phy, Order == "D_3__Rickettsiales") #337
# vec.mt <- rownames(otu_table(phy.mt))
# tax_table(phy.mt)
# length(rownames(otu_table(phy.mt))) ## 337 otus of CP

# (3) Unassigned
unique(tax_table(phy)[,'Kingdom']) ## only bacteria, then no need to exclude
phy.un <- subset_taxa(phy, Kingdom == "Unassigned")
vec.un <- rownames(otu_table(phy.un))
tax_table(phy.un)
length(rownames(otu_table(phy.un))) ## 1 otus

### exclude those vectors
phy #1346 taxa and 129 samples

sample_names(phy)
sample_variables(phy)

### pop taxa application
## get rid of CP and MT otus
### pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### let's do it!!!

phy.clean <- pop_taxa(phy, c(vec.un,vec.cp))

##before clean-up
phy 
sum(otu_table(phy)) ##1837395

##after clean-up
phy.clean
sum(otu_table(phy.clean)) #1835114

# checking procedure of whether the MT and CP otus are cleaned
# taxa_names(subset_taxa(phy, Order == "D_3__Rickettsiales"))  # before 337
# taxa_names(subset_taxa(phy, Order=="D_3__Chloroplast")) # before 130
# 
# taxa_names(subset_taxa(phy.clean , Order == "D_3__Rickettsiales")) # after 0
# taxa_names(subset_taxa(phy.clean , Family=="D_4__Mitochondria")) # after 0
# taxa_names(subset_taxa(phy.clean , Order == "D_3__Chloroplast")) # after 0

tax_table(phy.clean)

#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(phy.two)[,colnames(tax_table(phy.two))] <- gsub(tax_table(phy.two)[,colnames(tax_table(phy.two))],pattern="[A-Z]_[0-9]__",replacement="")
# phy.test4 <- phy.two %>% psmelt() 
# phy.test4

tax_table(phy.clean)[,colnames(tax_table(phy.clean))] <- gsub(tax_table(phy.clean)[,colnames(tax_table(phy.clean))],pattern="[A-Z]_[0-9]__",replacement="")

#' sample_data(phy.clean)$SampleID <- factor(sample_data(phy.clean)$SampleID, levels =target_PAB)


tax_table(phy.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

phy.clean    ## 28330
str(phy.clean)
otu_table(phy.clean)

phy.clean  ## 28330


# ## fix it in phy.clean object!!! pop_taxa does the work
# phy.clean <- pop_taxa(phy.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(phy.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement 
phy.clean.otu <- otu_table(phy.clean)
head(phy.clean.otu)
df.clean.otu <- data.frame(phy.clean.otu)
dim(df.clean.otu)
df.clean.otu$total <- apply(df.clean.otu, 1, sum)
head(df.clean.otu)
df.clean.otu <- tibble::rownames_to_column(df.clean.otu, var = 'OTU')


sample_names(phy.clean)


## Exclude false-positive OTUs
##2017 samples
# negative <- df.clean.otu[,c('OTU','Negative')]
# negative_0 <- subset(negative, Negative > 0) 
# neg.otu <- negative_0$OTU  ### 69 otu to be eliminated
# neg.otu
# length(neg.otu)
# tax_table(phy)[neg.otu]
# library(xlsx)
# write.xlsx(tax_table(phy)[neg.otu], 'neg.otu_for_fungi_taxonomy.xlsx')

## Remove reads with over 464 bp
bac.seq <- read.fasta(file = "dna-sequences_16s.fasta", as.string = TRUE, seqtype = "DNA")


#length
getLength(bac.seq)

## max length and min length
getLength(bac.seq)
min(getLength(bac.seq)) #39
max(getLength(bac.seq)) #447

##
bac.seq[which(getLength(bac.seq)<400)]
otu_less_400bp <- attr(bac.seq[which(getLength(bac.seq)<400)], "names")

phy.clean
phy.clean.ss <- pop_taxa(phy.clean,otu_less_400bp)
phy.clean.ss  ## 1343
sum(otu_table(phy.clean)) #1835114
sum(otu_table(phy.clean.ss)) #1832965


summarize_phyloseq(phy.clean.ss)


## Divide bacteria and archaea
arch.clean.ss <- subset_taxa(phy.clean.ss, Kingdom == "Archaea")
arch.clean.ss <- phyloseq::filter_taxa(arch.clean.ss, function(x) sum(x) != 0, TRUE)
sum(otu_table(arch.clean.ss)) #1579


bac.clean.ss <- subset_taxa(phy.clean.ss, Kingdom == "Bacteria") 
bac.clean.ss <- phyloseq::filter_taxa(bac.clean.ss, function(x) sum(x) != 0, TRUE) # 1339 OTUs
sum(otu_table(bac.clean.ss)) #1831386


###Fungal community
#### get unassigned vectors

# (3) Unassigned
unique(tax_table(fun)[,'Kingdom']) ## "k__Chromista", "k__Plantae", "Unassigned"
tax_table(fun)[,'Kingdom']
fun.un <- subset_taxa(fun, Kingdom %in% c("Unassigned","k__Chromista","k__Plantae"))
vec.un <- rownames(otu_table(fun.un))
tax_table(fun.un)
length(rownames(otu_table(fun.un))) ##  43

### exclude those vectors
fun  # 135 samples, 1629 taxa/ 1641 taxa

sample_names(fun)
sample_variables(fun)

### pop taxa application
## get rid of CP and MT otus
## pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### Clean it up!!!
fun.clean <- pop_taxa(fun, c(vec.un))

fun.clean #135 samples 1598 OTUs

#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(fun.two)[,colnames(tax_table(fun.two))] <- gsub(tax_table(fun.two)[,colnames(tax_table(fun.two))],pattern="[A-Z]_[0-9]__",replacement="")
# fun.test4 <- fun.two %>% psmelt() 
# fun.test4


tax_table(fun.clean)[,colnames(tax_table(fun.clean))] <- gsub(tax_table(fun.clean)[,colnames(tax_table(fun.clean))],pattern="[a-z]__",replacement="")

#' sample_data(fun.clean)$SampleID <- factor(sample_data(fun.clean)$SampleID, levels =target_PAB)

tax_table(fun.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

fun.clean 
str(fun.clean)
otu_table(fun.clean)

# ## fix it in fun.clean object!!! pop_taxa does the work
# fun.clean <- pop_taxa(fun.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(fun.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement 
fun.clean.otu <- otu_table(fun.clean)
head(fun.clean.otu)
df.clean.otu <- data.frame(fun.clean.otu)
dim(df.clean.otu)
df.clean.otu$total <- apply(df.clean.otu, 1, sum)
head(df.clean.otu)
df.clean.otu <- tibble::rownames_to_column(df.clean.otu, var = 'OTU')


# how about we delete major taxa in the negative sequence?
negative <- df.clean.otu[,c("OTU",'F204','F205','F206')]
negative
negative$total <- apply(negative[,-1], 1, sum)
negative_0 <- subset(negative, negative$total > 0) 
neg.otu <- negative_0$OTU  ### 58 otu to be eliminated
neg.otu
length(neg.otu)
tax_table(fun.clean)[neg.otu]
library(xlsx)
write.xlsx(tax_table(fun.clean)[neg.otu], 'fungal neg.otu_taxonomy_old DB.xlsx')

#### substract reads of false-positive OTUs
fun.clean.otu<-otu_table(fun.clean)

calib1 <- fun.clean.otu[which(rownames(fun.clean.otu)==neg.otu[1])]-negative_0$total[which(negative_0$OTU == neg.otu[1])]
calib2 <- fun.clean.otu[which(rownames(fun.clean.otu)==neg.otu[2])]-negative_0$total[which(negative_0$OTU == neg.otu[2])]
calib3 <- fun.clean.otu[which(rownames(fun.clean.otu)==neg.otu[3])]-negative_0$total[which(negative_0$OTU == neg.otu[3])]
calib4 <- fun.clean.otu[which(rownames(fun.clean.otu)==neg.otu[4])]-negative_0$total[which(negative_0$OTU == neg.otu[4])]
calib5 <- fun.clean.otu[which(rownames(fun.clean.otu)==neg.otu[5])]-negative_0$total[which(negative_0$OTU == neg.otu[5])]

calib.combined <- rbind(calib1,calib2, calib3,calib4,calib5)
calib.combined[which(calib.combined<0)] <- 0

fun.clean.otu.calib <- subset(fun.clean.otu, rownames(fun.clean.otu)%in%rownames(fun.clean.otu)[-which(rownames(fun.clean.otu)%in%negative_0$OTU)])
fun.clean.otu.calib<-rbind(fun.clean.otu.calib,calib.combined)

otu_table(fun.clean) <- otu_table(fun.clean.otu.calib, taxa_are_rows = T)
## checking oats
# oats <- df.clean.otu[,c('OTU','F_Oats1','F_Oats2','F_Oats3')]
# oats$total <- apply(negative[,-1], 1, sum)
# oats_0 <- subset(oats, oats$total > 0) 
# oats_0$OTU
# oats.otu <- oats_0$OTU
# 
# intersect(oats.otu, neg.otu)
# setdiff(oats.otu, neg.otu)
# union(oats.otu, neg.otu)

## exclude all negative taxa
# fun.clean <- pop_taxa(fun.clean, negative_0$OTU)
# fun.clean  ### 1629  -> 1624 otu
# class(fun.clean)

sum(otu_table(fun)) ##3154736
sum(otu_table(fun.clean)) ## 2983621 # calibrated 3154052 (2020 DB)/3051164(2017 DB)

### OTUs in positive samples
fun.clean.positive <- subset_samples(fun.clean, Diagnosis == "Positive")
fun.clean.positive <- phyloseq::filter_taxa(fun.clean.positive, function(x) sum(x) != 0, TRUE)

tax_table(fun.clean.positive)
length(rownames(tax_table(fun.clean.positive)))
otu_table(fun.clean.positive)
## Phew!!
## only get the clinical samples 
remains<- f.map$SampleID[-which(f.map$Diagnosis == "Positive" |f.map$Diagnosis == "Negative")]

fun.clean.ss <- subset_samples(fun.clean, SampleID %in% remains)
fun.clean.ss <- phyloseq::filter_taxa(fun.clean.ss, function(x) sum(x) != 0, TRUE) #
fun.clean.ss  ## 1623 #1626 #1595
sum(otu_table(fun.clean.ss)) # 2911517 #3076932 #2974044



### Fungal ITS1
###Fungal community
#### get unassigned vectors

# (3) Unassigned
fun.its1 <- fun
unique(tax_table(fun.its1)[,'Kingdom']) ## "k__Chromista", "k__Plantae", "Unassigned"
tax_table(fun.its1)[,'Kingdom']
fun.its1.un <- subset_taxa(fun.its1, Kingdom %in% c("Unassigned","k__Rhizaria","k__Metazoa"))
vec.un <- rownames(otu_table(fun.its1.un))
tax_table(fun.its1.un)
length(rownames(otu_table(fun.its1.un))) ##  2

### exclude those vectors
fun.its1  # 135 samples, 1412 taxa

sample_names(fun.its1)
sample_variables(fun.its1)

### pop taxa application
## get rid of CP and MT otus
## pop_taxa fun.its1ction
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### Clean it up!!!
fun.its1.clean <- pop_taxa(fun.its1, c(vec.un))

fun.its1.clean #129 samples 1410 OTUs

#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(fun.its1.two)[,colnames(tax_table(fun.its1.two))] <- gsub(tax_table(fun.its1.two)[,colnames(tax_table(fun.its1.two))],pattern="[A-Z]_[0-9]__",replacement="")
# fun.its1.test4 <- fun.its1.two %>% psmelt() 
# fun.its1.test4


tax_table(fun.its1.clean)[,colnames(tax_table(fun.its1.clean))] <- gsub(tax_table(fun.its1.clean)[,colnames(tax_table(fun.its1.clean))],pattern="[a-z]__",replacement="")

#' sample_data(fun.its1.clean)$SampleID <- factor(sample_data(fun.its1.clean)$SampleID, levels =target_PAB)

tax_table(fun.its1.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

fun.its1.clean 
str(fun.its1.clean)
otu_table(fun.its1.clean)

# ## fix it in fun.its1.clean object!!! pop_taxa does the work
# fun.its1.clean <- pop_taxa(fun.its1.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(fun.its1.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement 
fun.its1.clean.otu <- otu_table(fun.its1.clean)
head(fun.its1.clean.otu)
df.clean.otu <- data.frame(fun.its1.clean.otu)
dim(df.clean.otu)
df.clean.otu$total <- apply(df.clean.otu, 1, sum)
head(df.clean.otu)
df.clean.otu <- tibble::rownames_to_column(df.clean.otu, var = 'OTU')


# how about we delete major taxa in the negative sequence?
negative <- df.clean.otu[,c("OTU",'F204','F205','F206')]
negative
negative$total <- apply(negative[,-1], 1, sum)
negative_0 <- subset(negative, negative$total > 0) 
neg.otu <- negative_0$OTU  ### 58 otu to be eliminated
neg.otu
length(neg.otu)
tax_table(fun.its1.clean)[neg.otu]
library(xlsx)
write.xlsx(tax_table(fun.its1.clean)[neg.otu], 'fungal neg.otu_taxonomy_its1.xlsx')

#### substract reads of false-positive OTUs
fun.its1.clean.otu<-otu_table(fun.its1.clean)

calib1 <- fun.its1.clean.otu[which(rownames(fun.its1.clean.otu)==neg.otu[1])]-negative_0$total[which(negative_0$OTU == neg.otu[1])]
calib2 <- fun.its1.clean.otu[which(rownames(fun.its1.clean.otu)==neg.otu[2])]-negative_0$total[which(negative_0$OTU == neg.otu[2])]
calib3 <- fun.its1.clean.otu[which(rownames(fun.its1.clean.otu)==neg.otu[3])]-negative_0$total[which(negative_0$OTU == neg.otu[3])]
calib4 <- fun.its1.clean.otu[which(rownames(fun.its1.clean.otu)==neg.otu[4])]-negative_0$total[which(negative_0$OTU == neg.otu[4])]
#calib5 <- fun.its1.clean.otu[which(rownames(fun.its1.clean.otu)==neg.otu[5])]-negative_0$total[which(negative_0$OTU == neg.otu[5])]

calib.combined <- rbind(calib1,calib2, calib3,calib4)
calib.combined[which(calib.combined<0)] <- 0

fun.its1.clean.otu.calib <- subset(fun.its1.clean.otu, rownames(fun.its1.clean.otu)%in%rownames(fun.its1.clean.otu)[-which(rownames(fun.its1.clean.otu)%in%negative_0$OTU)])
fun.its1.clean.otu.calib<-rbind(fun.its1.clean.otu.calib,calib.combined)

otu_table(fun.its1.clean) <- otu_table(fun.its1.clean.otu.calib, taxa_are_rows = T)
## checking oats
# oats <- df.clean.otu[,c('OTU','F_Oats1','F_Oats2','F_Oats3')]
# oats$total <- apply(negative[,-1], 1, sum)
# oats_0 <- subset(oats, oats$total > 0) 
# oats_0$OTU
# oats.otu <- oats_0$OTU
# 
# intersect(oats.otu, neg.otu)
# setdiff(oats.otu, neg.otu)
# union(oats.otu, neg.otu)

## exclude all negative taxa
# fun.its1.clean <- pop_taxa(fun.its1.clean, negative_0$OTU)
# fun.its1.clean  ### 1629  -> 1624 otu
# class(fun.its1.clean)

sum(otu_table(fun.its1)) ##3031401
sum(otu_table(fun.its1.clean)) ## 3026565

### OTUs in positive samples
fun.its1.clean.positive <- subset_samples(fun.its1.clean, Diagnosis == "Positive")
sample_names(fun.its1.clean.positive)
fun.its1.clean.positive <- phyloseq::filter_taxa(fun.its1.clean.positive, function(x) sum(x) != 0, TRUE)

tax_table(fun.its1.clean.positive)
length(rownames(tax_table(fun.its1.clean.positive)))
otu_table(fun.its1.clean.positive)
## Phew!!
## only get the clinical samples 
remains<- f.map$SampleID[-which(f.map$Diagnosis == "Positive" |f.map$Diagnosis == "Negative")]

fun.its1.clean.ss <- subset_samples(fun.its1.clean, SampleID %in% remains)
fun.its1.clean.ss <- phyloseq::filter_taxa(fun.its1.clean.ss, function(x) sum(x) != 1, TRUE) #remove singleton
fun.its1.clean.ss  ## 1407 
sum(otu_table(fun.its1.clean.ss)) # 2947365

summarize_phyloseq(fun.its1.clean.ss)



##### get rid of otu of less than 100 reads

fun.seq <- read.fasta(file = "dna-sequences_its1.fasta", as.string = TRUE, seqtype = "DNA")
fun.seq$SH216041.07FU_AF138287_refs
otu_less_than_100bp <- attr(fun.seq[which(getLength(fun.seq)<100)], "names")
fun.clean.ss
fun.clean.ss <- pop_taxa(fun.clean.ss,otu_less_than_100bp)
fun.clean.ss ## 9437
sum(otu_table(fun.clean.ss)) ## 12191641



## Designating OTU id
bac.list <- bac.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
fun.list <- fun.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
fun.list.its1 <- fun.its1.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))

bac.list$number <- paste0('B',1:dim(bac.list)[1])
bac.list

bac.list$OTU_id <- ifelse(is.na(bac.list$Genus),ifelse(is.na(bac.list$Family),paste0(bac.list$number,'_o_',bac.list$Order),paste0(bac.list$number,'_f_',bac.list$Family)),paste0(bac.list$number,'_',bac.list$Genus))
bac.list$OTU_id

fun.list$number <- paste0('F',1:dim(fun.list)[1])
fun.list

fun.list$OTU_id <- ifelse(is.na(fun.list$Genus),ifelse(is.na(fun.list$Family),paste0(fun.list$number,'_o_',fun.list$Order),paste0(fun.list$number,'_f_',fun.list$Family)),paste0(fun.list$number,'_',fun.list$Genus))
fun.list$OTU_id

fun.list.its1$number <- paste0('F',1:dim(fun.list.its1)[1])
fun.list.its1

fun.list.its1$OTU_id <- ifelse(is.na(fun.list.its1$Genus),ifelse(is.na(fun.list.its1$Family),paste0(fun.list.its1$number,'_o_',fun.list.its1$Order),paste0(fun.list.its1$number,'_f_',fun.list.its1$Family)),paste0(fun.list.its1$number,'_',fun.list.its1$Genus))
fun.list.its1$OTU_id

bac.list

fun.list

fun.list.candida <- subset(fun.list, Genus == "Candida")
fun.list.asper <- subset(fun.list, Genus == "Aspergillus")

write.csv(fun.list.candida,"candida tax.csv")
write.csv(fun.list.asper,"aspergillus tax.csv")

otu.list <- rbind(bac.list, arch.list, fun.list)
dim(otu.list)
# write.xlsx(otu.list,'otu_id_book.xlsx')
fun.list.its1.OTU_id <-fun.list.its1[c('OTU','OTU_id')]


OTU_id.list <- rbind(bac.list[c('OTU','OTU_id')],arch.list[c('OTU','OTU_id')],fun.list[c('OTU','OTU_id')])
OTU_id.list$OTU_id
