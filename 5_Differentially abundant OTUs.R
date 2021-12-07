## Defining daOTUs
#Phyloseq files
bac.clean.ss

fun.clean.ss


zero.clean.filt <- phyloseq::filter_taxa(bac.clean.ss, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(zero.clean.filt) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(zero.clean.filt))

obj <- phyloseq_to_metagenomeSeq(zero.clean.filt)
obj 

## fitZig sample
obj <-  cumNorm(obj, p = cumNormStatFast(obj))
normFactor <-  normFactors(obj)
normFactor <-  log2(normFactor/median(normFactor) + 1)
settings <-  zigControl(maxit = 30, verbose = TRUE)

Type  <-  pData(obj)$Diagnosis

mod  <-  model.matrix(~Type)
colnames(mod)  <-  levels(Type)
colnames(mod)

res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)


zigFit = res$fit
finalMod= res$fit$design


contrast.matrix.1 =makeContrasts(control - RA, levels = finalMod)


fit2_diagnosis = contrasts.fit(zigFit, contrasts=contrast.matrix.1) 

fit3.1 = eBayes(fit2_diagnosis)

topTable(fit3.1, coef="control - RA")

res.1 <- topTable(fit3.1,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =1)
head(res.1)
dim(res.1)

## Taxonomy
tax_bac <- tax_table(bac.clean.ss) 

res.tax_bac <- merge(res.1,tax_bac, by = "row.names")

write.xlsx(res.tax_bac, 'bac_daOTU_diagnosis.xlsx')


## MA plot
log2AverageAbundance <- psmelt(bac.clean.ss) %>% group_by(OTU) %>% summarise(log2AverageAbundance=log2(mean(Abundance)))
log2AverageAbundance
Ta <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% select(Phylum,Class,Order,Family, Genus)
Ta <- unique(Ta)
Ta
resSig = res[!is.na(res$adj.P.Val), ]
resSig = data.frame(resSig)
head(resSig)
resSig <- tibble::rownames_to_column(resSig, 'OTU')
resSig <- left_join(resSig, Ab,by= c('OTU','OTU'))
resSig <- left_join(resSig,Ta,by=c('OTU','OTU'))


#### fungal community ###
zero.clean.filt <- phyloseq::filter_taxa(fun.clean.ss, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(zero.clean.filt) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(zero.clean.filt))

obj <- phyloseq_to_metagenomeSeq(zero.clean.filt)

## fitZig sample
obj <-  cumNorm(obj, p = cumNormStatFast(obj))
normFactor <-  normFactors(obj)
normFactor <-  log2(normFactor/median(normFactor) + 1)
settings <-  zigControl(maxit = 30, verbose = TRUE)

Type  <-  pData(obj)$Diagnosis

mod  <-  model.matrix(~Type)
colnames(mod)  <-  levels(Type)
colnames(mod)

res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)

class(res)

zigFit = res$fit
finalMod= res$fit$design


contrast.matrix.1 =makeContrasts(control - RA, levels = finalMod)


fit2_diagnosis = contrasts.fit(zigFit, contrasts=contrast.matrix.1) 

fit3.1 = eBayes(fit2_diagnosis)

topTable(fit3.1, coef="control - RA")

res.1 <- topTable(fit3.1,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =1)
head(res.1)
dim(res.1)

## Taxonomy
tax_fun <- tax_table(fun.clean.ss) 

res.tax_fun <- merge(res.1,tax_fun, by = "row.names")

write.xlsx(res.tax_fun, 'fun_daOTU_diagnosis.xlsx')


## MA plot
log2AverageAbundance <- psmelt(fun.clean.ss) %>% group_by(OTU) %>% summarise(log2AverageAbundance=log2(mean(Abundance)))
log2AverageAbundance
Ta <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% select(Phylum,Class,Order,Family, Genus)
Ta <- unique(Ta)
Ta
resSig = res[!is.na(res$adj.P.Val), ]
resSig = data.frame(resSig)
head(resSig)
resSig <- tibble::rownames_to_column(resSig, 'OTU')
resSig <- left_join(resSig, Ab,by= c('OTU','OTU'))
resSig <- left_join(resSig,Ta,by=c('OTU','OTU'))