### Finding common edge between or among groups


#Control (30) vs RA (30)
#Control
control_cor_df_padj <- read.table("2_0.4_edge_control.tsv", sep='\t', header =T)

control.pos <- droplevels(control_cor_df_padj[control_cor_df_padj$Cor>0,])
control.neg <- droplevels(control_cor_df_padj[control_cor_df_padj$Cor<0,])

control.pos.combine<- control.pos %>% tidyr::unite("Source_Target", Source:Target, remove = FALSE)
head(control.pos.combine$Source_Target)
length(control.pos.combine$Source_Target)#1281

control.neg.combine<- control.neg %>% tidyr::unite("Source_Target", Source:Target, remove = FALSE)
head(control.neg.combine$Source_Target)
length(control.neg.combine$Source_Target)#138

#RA
RA_cor_df_padj <- read.table("2_0.4_edge_RA_30.tsv", sep='\t', header =T)

RA.pos <- droplevels(RA_cor_df_padj[RA_cor_df_padj$Cor>0,])
RA.neg <- droplevels(RA_cor_df_padj[RA_cor_df_padj$Cor<0,])

RA.pos.combine<- RA.pos %>% tidyr::unite("Source_Target", Source:Target, remove = FALSE)
head(RA.pos.combine$Source_Target)
length(RA.pos.combine$Source_Target)#1540

RA.neg.combine<- RA.neg %>% tidyr::unite("Source_Target", Source:Target, remove = FALSE)
head(RA.neg.combine$Source_Target)
length(RA.neg.combine$Source_Target)#139

## Finding common edges
#Positive
common_pos_edge<-intersect(control.pos.combine$Source_Target, RA.pos.combine$Source_Target)
common_pos_edge

#Negative
common_neg_edge<-intersect(control.neg.combine$Source_Target, RA.neg.combine$Source_Target)


length(common_neg_edge)


##Positive to negative
common_pos.neg_edge<-intersect(control.pos.combine$Source_Target, RA.neg.combine$Source_Target)
common_pos.neg_edge

common_pos.neg_edge<-intersect(control.neg.combine$Source_Target, RA.pos.combine$Source_Target)
common_pos.neg_edge