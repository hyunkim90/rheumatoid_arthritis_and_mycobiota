#####LEFse analysis ####
### Generating input file
library(yingtools2)
requireNamespace(c("phyloseq","data.table"),quietly=TRUE)
pkgs <- c("splines","stats4","survival","mvtnorm","modeltools","coin","MASS")
missing.pkgs <- setdiff(pkgs,installed.packages()[,"Package"])
if (length(missing.pkgs)>0) {
  warning("YTWarning: R packages are needed for the LEFSE scripts to work: ",paste(missing.pkgs,collapse=", "))
}

### Divide above and belowground compartments
map <- read.table(file = 'Metadata_bac.tsv', sep = '\t', header = TRUE)
rownames(map) <- map$SampleID

map <- subset(map, map$Replication != "Negative")
nrow(map)

sample_data(bac.clean.ss) <- sample_data(map) 

(filt.sample <- sample_sums(bac.clean.ss) > 0)
sum(sample_sums(bac.clean.ss) <= 0)  ## 1 sample discarded
bac.clean.ss.f <- prune_samples(filt.sample, bac.clean.ss)
bac.clean.ss.f 


bac.meta<-sample_data(bac.clean.ss.f)
bac.meta$Division <- ifelse(bac.meta$Compartment %in% c("Leaf","Stem","Seed"),"Aboveground","Belowground")
sample_data(bac.clean.ss.f)<-sample_data(bac.meta)

bac.clean.ss.f.merge <- merge_samples(bac.clean.ss.f, "Replication")
meta.data.bac<-read.csv("fun.meta.merged.csv")
rownames(meta.data.bac) <-meta.data.bac$Sample

sample_data(bac.clean.ss.f.merge)<-sample_data(meta.data.bac)


####Bacterial OTU
phy =bac.clean.ss.f.merge
idtable <- bac.list[c("OTU","OTU_id")]
class = "Division"
subclass="Compartment"
subject=NA
anova.alpha=0.05
wilcoxon.alpha=0.05
lda.cutoff=2.0
wilcoxon.within.subclass=FALSE
one.against.one=FALSE
mult.test.correction=0
make.lefse.plots=FALSE
by_otus=FALSE
levels=phyloseq::rank_names(phy)


keepvars <- c(class,subclass,subject,"sample")
keepvars <- unique(keepvars[!is.na(keepvars)])
samp <- get.samp(phy)[,keepvars]
if (by_otus) { #perform by otu only
  otu <- get.otu.melt(phy,sample_data=FALSE)
  otu.levels <- otu %>% mutate(taxon=otu) %>%
    group_by(sample,taxon) %>% summarise(pctseqs=sum(pctseqs)) %>%
    mutate(taxon=gsub(" ","_",taxon))
  names(otu.levels)[2] <- "OTU"
  otu.levels <- merge(otu.levels, idtable, by = "OTU")
  otu.levels <- otu.levels[-c(1)]
  names(otu.levels)[3] <- "taxon"
  
} else { #divide by taxonomy
  otu <- get.otu.melt(phy,sample_data=FALSE)
  otu.list <- lapply(1:length(levels),function(i) {
    lvls <- levels[1:i]
    lvl <- levels[i]
    otu.level <- otu
    otu.level$taxon <- do.call(paste,c(lapply(lvls,function(l) otu[[l]]),sep="|"))
    otu.level$rank <- lvl
    otu.level2 <- otu.level %>% group_by(sample,taxon,rank) %>% summarise(pctseqs=sum(pctseqs)) %>% ungroup()
    return(otu.level2)
  })
  otu.levels <- bind_rows(otu.list) %>%
    mutate(taxon=gsub(" ","_",taxon))
}
otu.tbl <- otu.levels %>%
  dcast(sample~taxon,value.var="pctseqs",fill=0) %>%
  left_join(samp,by="sample") %>%
  select_(.dots=c(keepvars,lazyeval::interp(~everything())))
if (is.na(subject) | subject!="sample") {
  otu.tbl <- otu.tbl %>% select(-sample)
}
tbl <- otu.tbl %>% t()
write.table(tbl,"lefse_bacterial OTU_above and below_2.txt",quote=FALSE,sep="\t",col.names=FALSE)

write.table(tbl,"lefse_bacterial taxa_above and below_2.txt",quote=FALSE,sep="\t",col.names=FALSE)

###Fungal OTU
### Divide above and belowground compartments
fun.meta<-sample_data(fun.clean.ss2)
fun.meta$Division <- ifelse(fun.meta$Compartment %in% c("Leaf","Stem","Seed"),"Aboveground","Belowground")
fun.meta$Compartment[which(fun.meta$Compartment == "Soil")] <- "Bulk_soil"

sample_data(fun.clean.ss2)<-sample_data(fun.meta)

#Merge replicates
fun.clean.ss2.merge <- merge_samples(fun.clean.ss2, "Replication")
fun.meta.merged<-sample_data(fun.clean.ss2.merge)
write.csv(fun.meta.merged, "fun.meta.merged.csv")
meta.data.fun<-read.csv("fun.meta.merged.csv")
rownames(meta.data.fun) <-meta.data.fun$Sample
class(meta.data.fun)
sample_data(fun.clean.ss2.merge)<-sample_data(meta.data.fun)

## all replicate
phy = fun.clean.ss2.merge #fun.clean.ss
idtable <- fun.list[c("OTU","OTU_id")]
class = "Division"
subclass="Compartment"
subject=NA
anova.alpha=0.05
wilcoxon.alpha=0.05
lda.cutoff=2.0
wilcoxon.within.subclass=FALSE
one.against.one=FALSE
mult.test.correction=0
make.lefse.plots=FALSE
by_otus=F
levels=phyloseq::rank_names(phy)

keepvars <- c(class,subclass,subject,"sample")
keepvars <- unique(keepvars[!is.na(keepvars)])
samp <- get.samp(phy)[,keepvars]
if (by_otus) { #perform by otu only
  otu <- get.otu.melt(phy,sample_data=FALSE)
  otu.levels <- otu %>% mutate(taxon=otu) %>%
    group_by(sample,taxon) %>% summarise(pctseqs=sum(pctseqs)) %>%
    mutate(taxon=gsub(" ","_",taxon))
  names(otu.levels)[2] <- "OTU"
  otu.levels <- merge(otu.levels, idtable, by = "OTU")
  otu.levels <- otu.levels[-c(1)]
  names(otu.levels)[3] <- "taxon"
  
} else { #divide by taxonomy
  otu <- get.otu.melt(phy,sample_data=FALSE)
  otu.list <- lapply(1:length(levels),function(i) {
    lvls <- levels[1:i]
    lvl <- levels[i]
    otu.level <- otu
    otu.level$taxon <- do.call(paste,c(lapply(lvls,function(l) otu[[l]]),sep="|"))
    otu.level$rank <- lvl
    otu.level2 <- otu.level %>% group_by(sample,taxon,rank) %>% summarise(pctseqs=sum(pctseqs)) %>% ungroup()
    return(otu.level2)
  })
  otu.levels <- bind_rows(otu.list) %>%
    mutate(taxon=gsub(" ","_",taxon))
}
otu.tbl <- otu.levels %>%
  dcast(sample~taxon,value.var="pctseqs",fill=0) %>%
  left_join(samp,by="sample") %>%
  select_(.dots=c(keepvars,lazyeval::interp(~everything())))
if (is.na(subject) | subject!="sample") {
  otu.tbl <- otu.tbl %>% select(-sample)
}
tbl <- otu.tbl %>% t()
write.table(tbl,"lefse_fungal OTU_above and below_2.txt",quote=FALSE,sep="\t",col.names=FALSE)

write.table(tbl,"lefse_fungal taxa_above and below_2.txt",quote=FALSE,sep="\t",col.names=FALSE)
