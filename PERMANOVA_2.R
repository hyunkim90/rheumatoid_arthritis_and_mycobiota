## Effect of MTX
b.meta.edited.2 <- read.table("bacteria metadata_concentration.txt", sep = '\t', header =T)

row.names(b.meta.edited.2) <- b.meta.edited.2$SampleID 
b.meta<-b.meta.edited.2

b.meta$exercise_light <- as.character(b.meta$exercise_light)
b.meta$exercise_oderate <- as.character(b.meta$exercise_oderate)
b.meta$exercise_heavy <- as.character(b.meta$exercise_heavy)
b.meta$active_disease <- as.character(b.meta$active_disease)
b.meta$Duration <- as.character(b.meta$Duration)
b.meta$TNF_inhibitor <- as.character(b.meta$TNF_inhibitor)
b.meta$GC <- as.character(b.meta$GC)
b.meta$MTX <- as.character(b.meta$MTX)
b.meta$seropositive_RA <- as.character(b.meta$seropositive_RA)
b.meta$periodental_disease <-as.character(b.meta$periodental_disease)
b.meta$gu_edea_bleeding <-as.character(b.meta$gu_edea_bleeding)
b.meta$sore_on_the_tongue_hollow <-as.character(b.meta$sore_on_the_tongue_hollow)


b.meta$active_disease_2 <- b.meta$active_disease
b.meta$active_disease_2<-as.character(b.meta$active_disease_2)
b.meta$active_disease_2[which(b.meta$active_disease_2 %in% c("0", "1"))] <- "low"
b.meta$active_disease_2[which(b.meta$active_disease_2 %in% c("2", "3"))] <- "high"


b.meta$Ds_Duration <- as.numeric(as.character(b.meta$Ds_Duration))
b.meta$Ds_Duration[is.na(b.meta$Ds_Duration)] <- 0

b.meta$MTX_dose <- as.numeric(as.character(b.meta$MTX_dose))

b.meta$Diagnosis <- factor(b.meta$Diagnosis, levels = order.diagnosis)
b.meta$active_disease <- factor(b.meta$active_disease, levels = order.active_disease)

sample_data(bac.clean.log) <- sample_data(b.meta)
sample_data(bac.clean.ss) <- sample_data(b.meta)
sample_data(bac.clean.nolog) <-sample_data(b.meta)

## RA
bac.clean.ss.RA <- subset_samples(bac.clean.ss, Diagnosis == "RA") 
bac.clean.ss.RA <- phyloseq::filter_taxa(bac.clean.ss.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs

bac.clean.ss.RA.rel <- microbiome::transform(bac.clean.ss.RA, "compositional")

bac.clean.log.RA <- subset_samples(bac.clean.log, Diagnosis == "RA") 
bac.clean.log.RA <- phyloseq::filter_taxa(bac.clean.log.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs

bac.clean.nolog.RA <- subset_samples(bac.clean.nolog, Diagnosis == "RA") 
bac.clean.nolog.RA <- phyloseq::filter_taxa(bac.clean.nolog.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs


##Women
bac.clean.log.W <- subset_samples(bac.clean.log, Gender == "2")
bac.clean.log.W <- phyloseq::filter_taxa(bac.clean.log.W, function(x) sum(x) != 0, TRUE)


fun.clean.log.W <- subset_samples(fun.clean.log, Gender == "2")
fun.clean.log.W <- phyloseq::filter_taxa(fun.clean.log.W, function(x) sum(x) != 0, TRUE)


bac.clean.nolog.W <- subset_samples(bac.clean.log, Gender == "2")
bac.clean.nolog.W <- phyloseq::filter_taxa(bac.clean.nolog.W, function(x) sum(x) != 0, TRUE)


bac.clean.ss.rel.W <- subset_samples(bac.clean.ss.rel, Gender == "2")
bac.clean.ss.rel.W <- phyloseq::filter_taxa(bac.clean.ss.rel.W, function(x) sum(x) != 0, TRUE)




##RA and Women
bac.clean.log.RA.W <- subset_samples(bac.clean.log.RA, Gender == "2")
bac.clean.log.RA.W <- phyloseq::filter_taxa(bac.clean.log.RA.W, function(x) sum(x) != 0, TRUE)


bac.clean.nolog.RA.W <- subset_samples(bac.clean.nolog.RA, Gender == "2")
bac.clean.nolog.RA.W <- phyloseq::filter_taxa(bac.clean.nolog.RA.W, function(x) sum(x) != 0, TRUE)


bac.clean.ss.RA.rel.W <- subset_samples(bac.clean.ss.RA.rel, Gender == "2")
bac.clean.ss.RA.rel.W <- phyloseq::filter_taxa(bac.clean.ss.RA.rel.W, function(x) sum(x) != 0, TRUE)

b.meta
b.meta.ra <- subset(b.meta, Diagnosis == "RA")
b.meta.ra.women <- subset(b.meta.ra, Gender == "2")
b.meta.women <- subset(b.meta, Gender == "2")

### add periodental and periodental score
b.meta.periodental<-read.table("bacteria metadata_periodental.txt", sep = '\t', header =T)
rownames(b.meta.periodental) <- b.meta.periodental$SampleID
b.meta$periodental <- b.meta.periodental$periodental
head(b.meta)

b.meta$periodental_score <- b.meta.periodental$periodental_score
head(b.meta)

b.meta$periodental_score <- as.character(b.meta$periodental_score)
b.meta$periodental <- as.character(b.meta$periodental)


### add other variables
b.meta.other<-read.table("bacteria metadata_other.txt", sep = '\t', header =T)
rownames(b.meta.other) <- b.meta.other$SampleID
b.meta$CRP <- b.meta.other$CRP
b.meta$WBC <- b.meta.other$WBC
b.meta$ESR <- b.meta.other$ESR
b.meta$Hb <- b.meta.other$Hb
b.meta$Plt <- b.meta.other$Plt
b.meta$BUN <- b.meta.other$BUN
b.meta$Cr <- b.meta.other$Cr
b.meta$Total_Cholesterol <- b.meta.other$Total_cholesterol
b.meta$Triglyceride <- b.meta.other$Triglyceride
b.meta$HDL <- b.meta.other$HDL
head(b.meta)

b.meta$periodental_score <- b.meta.periodental$periodental_score
head(b.meta)


##PERMANOVA

##All
b.otu <- otu_table(bac.clean.ss.rel)
b.otu <- otu_table(bac.clean.nolog)
b.otu <- otu_table(bac.clean.log)

b.permanova <- adonis(formula = t(b.otu) ~ (MTX+MTX_dose), data = b.meta, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (periodental_disease+gu_edea_bleeding+sore_on_the_tongue_hollow), data = b.meta, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (periodental+periodental_score), data = b.meta, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (Age+BMI+Total_cholesterol+Duration+ HDL+Triglyceride+RA_factor+anti_CCP+CRP+ESR), data = b.meta, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (Diagnosis), data = b.meta, permutations=9999, method = "bray")
b.permanova

#Women
b.otu <- otu_table(bac.clean.ss.rel.W)
b.otu <- otu_table(bac.clean.nolog.W)
b.otu <- otu_table(bac.clean.log.W)

b.permanova <- adonis(formula = t(b.otu) ~ (MTX+MTX_dose), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (enstruation), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (periodental_disease+gu_edea_bleeding+sore_on_the_tongue_hollow), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (periodental+periodental_score), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (Age+BMI+Total_cholesterol+Duration+ HDL+Triglyceride+RA_factor+anti_CCP+CRP+ESR), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (Diagnosis), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

##Fungi

f.otu <- otu_table(fun.clean.log.W)

f.permanova <- adonis(formula = t(f.otu) ~ (Diagnosis+Age+BMI+Total_cholesterol+Duration+ HDL+Triglyceride+RA_factor+anti_CCP+CRP+ESR), data = f.meta.women, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu) ~ (Diagnosis), data = f.meta.women, permutations=9999, method = "bray")
f.permanova

f.otu <- otu_table(fun.clean.log)
f.permanova <- adonis(formula = t(f.otu) ~ (Age+BMI+Total_cholesterol+Duration+ HDL+Triglyceride+RA_factor+anti_CCP+CRP+ESR), data = f.meta, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu) ~ (Diagnosis), data = f.meta, permutations=9999, method = "bray")
f.permanova


## RA
b.otu.ra <- otu_table(bac.clean.ss.RA.rel)
b.otu.ra <- otu_table(bac.clean.nolog.RA)
b.otu.ra <- otu_table(bac.clean.log.RA)

b.permanova <- adonis(formula = t(b.otu.ra) ~ (MTX+MTX_dose), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu.ra) ~ (periodental_disease+gu_edea_bleeding+sore_on_the_tongue_hollow), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu.ra) ~ (periodental+periodental_score), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu.ra) ~ (Age+BMI+Total_cholesterol+Duration+ HDL+Triglyceride+RA_factor+anti_CCP+CRP+ESR), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova


## RA and women
b.otu.ra <- otu_table(bac.clean.ss.RA.rel.W)
b.otu.ra <- otu_table(bac.clean.nolog.RA.W)
b.otu.ra <- otu_table(bac.clean.log.RA.W)


b.permanova <- adonis(formula = t(b.otu.ra) ~ (enstruation), data = b.meta.ra.women, permutations=9999, method = "bray")
b.permanova


b.permanova <- adonis(formula = t(b.otu.ra) ~ (MTX+MTX_dose), data = b.meta.ra.women, permutations=9999, method = "bray")
b.permanova


b.permanova <- adonis(formula = t(b.otu.ra) ~ (periodental_disease+gu_edea_bleeding+sore_on_the_tongue_hollow), data = b.meta.ra.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu.ra) ~ (periodental+periodental_score), data = b.meta.ra.women, permutations=9999, method = "bray")
b.permanova


## remove samples which do not have blood test results
## Lipid
rm1<-b.meta.ra$SampleID[-which(b.meta.ra$Total_Cholesterol == 0)]
rm2<-b.meta.ra$SampleID[-which(b.meta.ra$Triglyceride == 0)]
rm3<-b.meta.ra$SampleID[-which(b.meta.ra$HDL == 0)]
rm4<-b.meta.ra$SampleID[-which(b.meta.ra$CRP == 0)]
rm5<-b.meta.ra$SampleID[-which(b.meta.ra$WBC == 0)]
rm6<-b.meta.ra$SampleID[-which(b.meta.ra$ESR == 0)]
rm7<-b.meta.ra$SampleID[-which(b.meta.ra$Hb == 0)]
rm8<-b.meta.ra$SampleID[-which(b.meta.ra$BUN == 0)]
rm9<-b.meta.ra$SampleID[-which(b.meta.ra$Cr == 0)]

b.otu.ra <- otu_table(bac.clean.ss.RA.rel.W)
b.otu.ra <- otu_table(bac.clean.nolog.RA.W)
b.otu.ra <- otu_table(bac.clean.log.RA.W)

length(b.meta$SampleID[which(b.meta$Total_Cholesterol>200)])

### contribution of total cholesterol
b.otu.ra.t <- subset(t(b.otu.ra), rownames(t(b.otu.ra)) %in% as.character(rm1))
b.meta.ra.edit <- subset(b.meta.ra.women, rownames(b.meta.ra.women)%in%rm1)

b.permanova <- adonis(formula = b.otu.ra.t ~ (Total_Cholesterol), data = b.meta.ra.edit, permutations=9999, method = "bray")
b.permanova

### contribution of Triglyceride
b.otu.ra.t <- subset(t(b.otu.ra), rownames(t(b.otu.ra)) %in% as.character(rm2))
b.meta.ra.edit <- subset(b.meta.ra, rownames(b.meta.ra)%in%rm2)

b.permanova <- adonis(formula = b.otu.ra.t ~ (Triglyceride), data = b.meta.ra.edit, permutations=9999, method = "bray")
b.permanova

### contribution of HDL
b.otu.ra.t <- subset(t(b.otu.ra), rownames(t(b.otu.ra)) %in% as.character(rm3))
b.meta.ra.edit <- subset(b.meta.ra.women, rownames(b.meta.ra.women)%in%rm3)

b.permanova <- adonis(formula = b.otu.ra.t ~ (HDL), data = b.meta.ra.edit, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = b.otu.ra.t ~ (Total_Cholesterol+HDL+Triglyceride), data = b.meta.ra.edit, permutations=9999, method = "bray")
b.permanova

### contribution of CRP
b.otu.ra.t <- subset(t(b.otu.ra), rownames(t(b.otu.ra)) %in% as.character(rm4))
b.meta.ra.edit <- subset(b.meta.ra, rownames(b.meta.ra)%in%rm4)

b.permanova <- adonis(formula = b.otu.ra.t ~ (CRP), data = b.meta.ra.edit, permutations=9999, method = "bray")
b.permanova

### contribution of wbc
b.otu.ra.t <- subset(t(b.otu.ra), rownames(t(b.otu.ra)) %in% as.character(rm5))
b.meta.ra.edit <- subset(b.meta.ra, rownames(b.meta.ra)%in%rm5)

b.permanova <- adonis(formula = b.otu.ra.t ~ (WBC), data = b.meta.ra.edit, permutations=9999, method = "bray")
b.permanova

### contribution of Hb
b.otu.ra.t <- subset(t(b.otu.ra), rownames(t(b.otu.ra)) %in% as.character(rm7))
b.meta.ra.edit <- subset(b.meta.ra, rownames(b.meta.ra)%in%rm7)

b.permanova <- adonis(formula = b.otu.ra.t ~ (Hb), data = b.meta.ra.edit, permutations=9999, method = "bray")
b.permanova


### contribution of ESR
b.otu.ra.t <- subset(t(b.otu.ra), rownames(t(b.otu.ra)) %in% as.character(rm6))
b.meta.ra.edit <- subset(b.meta.ra, rownames(b.meta.ra)%in%rm6)

b.permanova <- adonis(formula = b.otu.ra.t ~ (ESR), data = b.meta.ra.edit, permutations=9999, method = "bray")
b.permanova


### contribution of Plt
b.otu.ra.t <- subset(t(b.otu.ra), rownames(t(b.otu.ra)) %in% as.character(rm5))
b.meta.ra.edit <- subset(b.meta.ra, rownames(b.meta.ra)%in%rm5)

b.permanova <- adonis(formula = b.otu.ra.t ~ (Plt), data = b.meta.ra.edit, permutations=9999, method = "bray")
b.permanova

### contribution of BUN
b.otu.ra.t <- subset(t(b.otu.ra), rownames(t(b.otu.ra)) %in% as.character(rm8))
b.meta.ra.edit <- subset(b.meta.ra, rownames(b.meta.ra)%in%rm8)

b.permanova <- adonis(formula = b.otu.ra.t ~ (BUN), data = b.meta.ra.edit, permutations=9999, method = "bray")
b.permanova


### contribution of Cr
b.otu.ra.t <- subset(t(b.otu.ra), rownames(t(b.otu.ra)) %in% as.character(rm9))
b.meta.ra.edit <- subset(b.meta.ra, rownames(b.meta.ra)%in%rm9)

b.permanova <- adonis(formula = b.otu.ra.t ~ (Cr), data = b.meta.ra.edit, permutations=9999, method = "bray")
b.permanova