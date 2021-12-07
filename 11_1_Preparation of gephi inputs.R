
### sparcc let's go!!
### making gephi input!!
### and also making edge and node list
library(reshape2)
library(tidyr)
## import edge list
## 

### Gephi input
# Control 
cor_mat <- read.table(file = 'control_median_correlation.tsv', sep = '\t', header = TRUE)
rownames(cor_mat) <- cor_mat$X
cor_mat <- cor_mat[-c(1)]

p_mat <- read.table(file = 'control_pvalues.tsv', sep = '\t', header = TRUE)
rownames(p_mat) <- p_mat$X
p_mat <- p_mat[-c(1)]

get_edge_node_withp <- function(cor_mat,p_mat,edgename,nodename,threshold){
  head(cor_mat)
  dim(cor_mat)
  
  any(is.na(cor_mat))
  cor_mat[upper.tri(cor_mat)] <- NA
  df_cor <- reshape2::melt(cor_mat, varnames = c('Source', 'Target'), na.rm = TRUE)
  
  dim(df_cor) ## 258840
  head(df_cor)
  unique(df_cor$variable)
  tail(df_cor)
  any(is.na(df_cor))
  
  # threshold <-  0.3
  df_cor_sig <- df_cor %>% filter(abs(value)>=threshold)  ## 0.3 is 22 edges   ## 0.2 is 124 edges  ## 0.25 is 48
  print(dim(df_cor_sig))  ## total 94 edges
  # colnames(df_cor_sig) <- c('Source','Target','Cor')
  
  head(p_mat)
  dim(p_mat)
  
  any(is.na(p_mat))
  p_mat[upper.tri(p_mat)] <- NA
  df_p <- reshape2::melt(p_mat, varnames = c('Source', 'Target'), na.rm = TRUE)
  
  dim(df_p) ## 4186
  head(df_p)
  tail(df_p)
  any(is.na(df_p))
  
  df_cor_p <- left_join(df_cor, df_p, by=c('OTU_id'='OTU_id','variable'='variable'))
  head(df_cor_p)
  
  colnames(df_cor_p) <- c('Source','Target','Cor','pseudo_p')
  df_cor_sig <- df_cor_p %>% filter(pseudo_p < 0.05) %>% filter(abs(Cor)>threshold)  ## 15
  dim(df_cor_sig)   ## 94
  head(df_cor_sig)
  df_cor_sig$positive <- df_cor_sig$Cor > 0
  
  write.table(df_cor_sig, file=edgename, quote=FALSE, sep='\t', row.names = F)
  v.source <- df_cor_sig$Source
  v.target <- df_cor_sig$Target
  v.id <- union(v.source,v.target)
  length(v.id)
  
  df.node <- data.frame(id=v.id)
  df.node$label <- df.node$id
  #df.node$sep <- df.node$id
  
  df.node$Bacteria <- ifelse(grepl("^B",df.node$id),"Bacteria","Fungi")
 
  write.table(df.node, file=nodename, quote=FALSE, sep='\t', row.names = F)
}



# get_edge_node_withoutp<- function(cor_mat,edgename,nodename,threshold){
#   head(cor_mat)
#   dim(cor_mat)
#   
#   any(is.na(cor_mat))
#   cor_mat[upper.tri(cor_mat)] <- NA
#   df_cor <- reshape2::melt(cor_mat, varnames = c('Source', 'Target'), na.rm = TRUE)
#   
#   dim(df_cor) ## 258840
#   head(df_cor)
#   unique(df_cor$variable)
#   tail(df_cor)
#   any(is.na(df_cor))
#   
#   # threshold <-  0.3
#   df_cor_sig <- df_cor %>% filter(abs(value)>=threshold)  ## 0.3 is 22 edges   ## 0.2 is 124 edges  ## 0.25 is 48
#   print(dim(df_cor_sig))  ## total 22 edges
#   colnames(df_cor_sig) <- c('Source','Target','Cor')
#   head(df_cor_sig)
#   df_cor_sig$positive <- df_cor_sig$Cor > 0
#   
#   write.table(df_cor_sig, file=edgename, quote=FALSE, sep='\t', row.names = F)
#   v.source <- df_cor_sig$Source
#   v.target <- df_cor_sig$Target
#   v.id <- union(v.source,v.target)
#   length(v.id)
#   
#   df.node <- data.frame(id=v.id)
#   df.node$label <- df.node$id
#   df.node$Bacteria <- ifelse(grepl("^B",df.node$id),TRUE,FALSE)
#   df.node <- df.node %>% left_join(fb.list4, by=c('id'='id'))
#   (df.node)
# 
#   write.table(df.node, file=nodename, quote=FALSE, sep='\t', row.names = F)
# }
# 
# df.node

get_gephi_input_withp <- function(keyword,front_number, threshold){
  # keyword <- 'mer.rel.ss.5'
  coreword <- paste0(keyword,'_median_correlation','.tsv')
  coreword
  cor_mat <- read.table(file = coreword, sep = '\t', header = TRUE)
  keyword
  pword <- paste0(keyword,'_pvalues','.tsv')
  p_mat <- read.table(file = pword, sep = '\t', header = TRUE)
  edgeword <- paste0(as.character(front_number),'_',as.character(threshold),'_','edge_',keyword,'.tsv')
  # threshold <- 0.2
  # front_number <- 1
  edgeword
  nodeword <- paste0(as.character(front_number),'_',as.character(threshold),'_','node_',keyword,'.tsv')
  nodeword
  
  get_edge_node_withp(cor_mat,p_mat,edgeword,nodeword,threshold)
}


#Control
get_gephi_input_withp('control',3,0.5)
get_gephi_input_withp('control',1,0.3)
get_gephi_input_withp('control',2,0.4)


#RA
get_gephi_input_withp('RA',1,0.5)
get_gephi_input_withp('RA',2,0.3)
get_gephi_input_withp('RA',3,0.4)



##RA 30 samples

get_gephi_input_withp('RA_30',3,0.5)
get_gephi_input_withp('RA_30',1,0.3)
get_gephi_input_withp('RA_30',2,0.4)


# ##RA 20 samples (non-treated)
# 
# get_gephi_input_withp('RA_20',3,0.5)
# get_gephi_input_withp('RA_20',1,0.3)
# get_gephi_input_withp('RA_20',2,0.4)
# 
# get_gephi_input_withp('control_20',3,0.5)
# get_gephi_input_withp('control_20',1,0.3)
# get_gephi_input_withp('control_20',2,0.4)
# 
# ### TNF inhibitor-treated subject
# 
# get_gephi_input_withp('TNF',3,0.5)
# get_gephi_input_withp('TNF',1,0.3)
# get_gephi_input_withp('TNF',2,0.4)
# 
# ### MTX inhibitor-treated subject
# 
# get_gephi_input_withp('MTX',3,0.5)
# get_gephi_input_withp('MTX',1,0.3)
# get_gephi_input_withp('MTX',2,0.4)
# 
# 
# 
# ### Network of ITS 1
# get_gephi_input_withp('control_its1',3,0.5)
# get_gephi_input_withp('control_its1',1,0.3)
# get_gephi_input_withp('control_its1',2,0.4)
# 
# get_gephi_input_withp('RA_its1',3,0.5)
# get_gephi_input_withp('RA_its1',1,0.3)
# get_gephi_input_withp('RA_its1',2,0.4)