## correlation between OTU abundance and RA-related factors
## Normalized abundance
bac.clean.ss
fun.clean.ss

## extract abundance table
otu.bac.rel <- otu_table(bac.clean.ss.rel)
otu.bac.rel <- data.frame(otu.bac.rel)
otu.fun.rel <- otu_table(fun.clean.ss.rel)
otu.fun.rel <- data.frame(otu.fun.rel)

##column 'OTU'
otu.bac.rel$OTU <- rownames(otu.bac.rel)
otu.fun.rel$OTU <- rownames(otu.fun.rel)

## attach otu id in otu table 
otu.bac.rel.id <- merge(otu.bac.rel, bac.list.OTU_id, by = "OTU")
otu.fun.rel.id <- merge(otu.fun.rel, fun.list.OTU_id, by = "OTU")
rownames(otu.bac.rel.id) <- otu.bac.rel.id$OTU_id
rownames(otu.fun.rel.id) <- otu.fun.rel.id$OTU_id

##remove column 'OTU' 
otu.bac.rel.id <- otu.bac.rel.id[-c(1)]
otu.fun.rel.id <- otu.fun.rel.id[-c(1)]
##remove column 'OTU_id'
otu.bac.rel.id <- otu.bac.rel.id[-c(130)]
otu.fun.rel.id <- otu.fun.rel.id[-c(130)]

### attach metadata files
#bacteria
b.meta.numeric <- b.meta %>% select(Age, BMI, RA_factor, anti_CCP)

otu.bac.rel.meta <- merge(b.meta.numeric, t(otu.bac.rel.id), by = 'row.names')
head(otu.bac.rel.meta)
row.names(otu.bac.rel.meta) <- otu.bac.rel.meta$Row.names
otu.bac.rel.meta<-otu.bac.rel.meta[-c(1)]

cor_bac <- Hmisc::rcorr(as.matrix(otu.bac.rel.meta), type="spearman")
M.bac <- cor_bac$r
p_mat.bac <- cor_bac$P

cor.bac.OTU <- flattenCorrMatrix(M.bac,p_mat.bac)
cor.bac.OTU$sig <- gtools::stars.pval(cor.bac.OTU$p)

head(cor.bac.OTU)

RA_variables<-colnames(otu.bac.rel.meta[c(1:4)])
OTUnames.bac<-colnames(otu.bac.rel.meta[-c(1:4)])

cor.bac.OTU.trim <- subset(cor.bac.OTU, row%in%RA_variables)
cor.bac.OTU.trim <- subset(cor.bac.OTU.trim, column%in%OTUnames.bac)

cor.bac.OTU.trim$cor<- round(cor.bac.OTU.trim$cor, 5)
head(cor.bac.OTU.trim)


##Only significant
sig_cor.bac.OTU.trim <- subset(cor.bac.OTU.trim, p < 0.05)
sig_cor.bac.OTU.trim.trans <- dcast(sig_cor.bac.OTU.trim , row ~ column, value.var = "cor")
sig_cor.bac.OTU.trim.trans.p <- dcast(sig_cor.bac.OTU.trim , row ~ column, value.var = "p")
write.xlsx(t(sig_cor.bac.OTU.trim.trans), "relative abundance_siginificant correlation between bacterial OTU and quantitative factors.xlsx")




#Fungi
f.meta.numeric <- f.meta %>% select(Age, BMI, RA_factor, anti_CCP)
otu.fun.rel.meta <- merge(f.meta.numeric, t(otu.fun.rel.id), by = 'row.names')

row.names(otu.fun.rel.meta) <- otu.fun.rel.meta$Row.names
otu.fun.rel.meta<-otu.fun.rel.meta[-c(1)]

cor_fun <- Hmisc::rcorr(as.matrix(otu.fun.rel.meta), type="spearman")
M.fun <- cor_fun$r
p_mat.fun <- cor_fun$P

cor.fun.OTU <- flattenCorrMatrix(M.fun,p_mat.fun)
cor.fun.OTU$sig <- gtools::stars.pval(cor.fun.OTU$p)

head(cor.fun.OTU)

RA_variables<-colnames(otu.fun.rel.meta[c(1:4)])
OTUnames.fun<-colnames(otu.fun.rel.meta[-c(1:4)])

cor.fun.OTU.trim <- subset(cor.fun.OTU, row%in%RA_variables)
cor.fun.OTU.trim <- subset(cor.fun.OTU.trim, column%in%OTUnames.fun)

cor.fun.OTU.trim$cor<- round(cor.fun.OTU.trim$cor, 5)
head(cor.fun.OTU.trim)


##Only significant
sig_cor.fun.OTU.trim <- subset(cor.fun.OTU.trim, p < 0.05)
sig_cor.fun.OTU.trim.trans <- dcast(sig_cor.fun.OTU.trim , row ~ column, value.var = "cor")
sig_cor.fun.OTU.trim.trans.p <- dcast(sig_cor.fun.OTU.trim , row ~ column, value.var = "p")
write.xlsx(t(sig_cor.fun.OTU.trim.trans), "relative abundance_siginificant correlation between fungal OTU and quantitative factors.xlsx")



##Correlation
