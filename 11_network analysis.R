### Network

## Normalized abundance
bac.clean.nolog
fun.clean.nolog

## extract abundance table
otu.bac <- otu_table(bac.clean.nolog)
otu.fun <- otu_table(fun.clean.nolog)[,sample.order]
# otu.fun <- otu_table(fun.its1.clean.nolog)

## Unifying Sample names
sample_names(otu.fun) <- sample_names(otu.bac)
head(otu.fun)

## attach otu id in otu table 
otu.bac <- data.frame(otu.bac)
otu.fun <- data.frame(otu.fun)
otu.bac$OTU <- rownames(otu.bac)
otu.fun$OTU <- rownames(otu.fun)

otu.bac.id <- merge(otu.bac, bac.list.OTU_id, by = "OTU")
otu.fun.id <- merge(otu.fun, fun.list.OTU_id, by = "OTU")
# otu.fun.id <- merge(otu.fun, fun.list.its1.OTU_id, by = "OTU")

## merge otu table
otu.tab.merged <- rbind(otu.bac.id, otu.fun.id)
rownames(otu.tab.merged) <- otu.tab.merged$OTU_id

##remove column 'OTU' 
otu.tab.merged <- otu.tab.merged[-c(1)]

##remove column 'OTU_id'
otu.tab.merged <- otu.tab.merged[-c(130)]

#write.table(otu.tab.merged, "otu.tab.merged.tsv", sep='\t', quote=F)

## divide otu tables by diagnosis (control and RA) 
control.subject<-b.meta$SampleID[which(b.meta$Diagnosis == "control")]
RA.subject<-b.meta$SampleID[which(b.meta$Diagnosis == "RA")]

## Control
otu.tab.merged.control <- subset(t(otu.tab.merged), rownames(t(otu.tab.merged)) %in% control.subject)
otu.tab.merged.control <- t(otu.tab.merged.control)

remained.OTUs.control<-rownames(otu.tab.merged.control)[which(rowSums(otu.tab.merged.control) != 0)]

otu.tab.merged.control <- subset(otu.tab.merged.control, rownames(otu.tab.merged.control) %in% remained.OTUs.control)

nrow(otu.tab.merged.control) #1302 (ITS2), 1276 (ITS1)
write.table(otu.tab.merged.control, "otu.tab.merged.control.tsv", sep='\t', quote=F)
# write.table(otu.tab.merged.control, "otu.tab.merged.control_its1.tsv", sep='\t', quote=F)

##RA (all samples)
otu.tab.merged.RA <- subset(t(otu.tab.merged), rownames(t(otu.tab.merged)) %in% RA.subject)
otu.tab.merged.RA <- t(otu.tab.merged.RA)
remained.OTUs.RA<-rownames(otu.tab.merged.RA)[which(rowSums(otu.tab.merged.RA) != 0)]
otu.tab.merged.RA <- subset(otu.tab.merged.RA, rownames(otu.tab.merged.RA) %in% remained.OTUs.RA)
nrow(otu.tab.merged.RA)

# write.table(otu.tab.merged.RA, "otu.tab.merged.RA_its1.tsv", sep='\t', quote=F)

## Random sampling on RA samples
library(MASS)
set.seed(123)

otu.tab.merged.RA <- subset(t(otu.tab.merged), rownames(t(otu.tab.merged)) %in% RA.subject)
dim(otu.tab.merged.RA) #nrow 99/ ncol 2933
simple.random.sampling <- sample(1:nrow(otu.tab.merged.RA), 30)

otu.tab.merged.RA_srs <- otu.tab.merged.RA[simple.random.sampling, ]
dim(otu.tab.merged.RA_srs) #nrow 30/ ncol 2745 (ITS1) // 30 2933 (ITS2)

otu.tab.merged.RA_srs <- t(otu.tab.merged.RA_srs)
remained.OTUs.RA<-rownames(otu.tab.merged.RA_srs)[which(rowSums(otu.tab.merged.RA_srs) != 0)]
otu.tab.merged.RA_srs <- subset(otu.tab.merged.RA_srs, rownames(otu.tab.merged.RA_srs) %in% remained.OTUs.RA)
nrow(otu.tab.merged.RA_srs) #1428 (ITS2) 
write.table(otu.tab.merged.RA_srs, "otu.tab.merged.RA_srs.tsv", sep='\t', quote=F)
# write.table(otu.tab.merged.RA_srs, "otu.tab.merged.RA_srs_its1.tsv", sep='\t', quote=F)









# ### control vs non-treated RA patients
# 
# ## subsampling 20 control subject
# otu.tab.merged.control <- subset(t(otu.tab.merged), rownames(t(otu.tab.merged)) %in% control.subject)
# dim(otu.tab.merged.control) #nrow 30/ ncol 2933
# simple.random.sampling.con <- sample(1:nrow(otu.tab.merged.control), 20)
# 
# otu.tab.merged.control_srs <- otu.tab.merged.control[simple.random.sampling.con, ]
# dim(otu.tab.merged.control_srs) #nrow 20/ ncol 2933
# 
# 
# otu.tab.merged.control_srs <- t(otu.tab.merged.control_srs)
# write.table(otu.tab.merged.control_srs, "otu.tab.merged.control_srs.tsv", sep='\t', quote=F)
# 
# 
# ## non-treated patients
# otu.tab.merged.nontreated <- subset(t(otu.tab.merged), rownames(t(otu.tab.merged)) %in% non.treated.RA)
# write.table(t(otu.tab.merged.nontreated), "otu.tab.merged.nontratedRA_srs.tsv", sep='\t', quote=F)
# ### Random sampling with subjects not treated with MTX
# length(b.meta$SampleID[which(b.meta$MTX == 0 & b.meta$TNF_inhibitor == 0)]) ##RA effect
# length(b.meta$SampleID[which(b.meta$MTX == 1 & b.meta$TNF_inhibitor == 0)])
# non.treated.RA <-b.meta$SampleID[which(b.meta$MTX == 0 & b.meta$TNF_inhibitor == 0)]
# b.meta
# 
# 
# 
# 
# ## Treatment methods
# ## TNF inhibitor
# otu.tab.merged.TNF <- subset(t(otu.tab.merged), rownames(t(otu.tab.merged)) %in% TNF_subject)
# otu.tab.merged.TNF <- t(otu.tab.merged.TNF)
# write.table(otu.tab.merged.TNF, "otu.tab.merged.TNF.tsv", sep='\t', quote=F)
# 
# 
# ##MTX
# otu.tab.merged.MTX <- subset(t(otu.tab.merged), rownames(t(otu.tab.merged)) %in% MTX_subject)
# dim(otu.tab.merged.MTX) #nrow 99/ ncol 2933
# simple.random.sampling <- sample(1:nrow(otu.tab.merged.MTX), 20)
# 
# otu.tab.merged.MTX_srs <- otu.tab.merged.MTX[simple.random.sampling, ]
# dim(otu.tab.merged.MTX_srs) #nrow 30/ ncol 2933
# 
# 
# otu.tab.merged.MTX_srs <- t(otu.tab.merged.MTX_srs)
# write.table(otu.tab.merged.MTX_srs, "otu.tab.merged.MTX.tsv", sep='\t', quote=F)
