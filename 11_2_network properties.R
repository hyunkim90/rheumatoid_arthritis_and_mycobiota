### Network properties
## read cor and p

control_cor_df_padj <- read.table("2_0.4_edge_control.tsv", sep='\t', header =T)
control_cor_df_padj <- read.table("2_0.4_edge_control_its1.tsv", sep='\t', header =T)

nodeattrib_control_combine <- data.frame(node=union(control_cor_df_padj$Source,control_cor_df_padj$Target))
nodeattrib_control_combine$kingdom <- 0

for (i in as.character(nodeattrib_control_combine$node))
{
  if (i %in% otu.bac.id$OTU_id == TRUE)
  {nodeattrib_control_combine[nodeattrib_control_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_control_combine[nodeattrib_control_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_control_combine) <- as.character(nodeattrib_control_combine$node)
nodeattrib_control_combine$kingdom



all_control_net <- graph_from_data_frame(control_cor_df_padj,direct=F, vertices=nodeattrib_control_combine)

## Number of nodes
length(V(all_control_net)) #716 (ITS2) / 370 (ITS1)

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_control_net)))) #349 / 251
length(grep("^F",names(V(all_control_net)))) #367 / 119


## Connections 
bb_occur_control <- droplevels(control_cor_df_padj[with(control_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_control) #379 /374

ff_occur_control <- droplevels(control_cor_df_padj[with(control_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_control) #617 /94

fb_occur_control <- droplevels(control_cor_df_padj[with(control_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_control) #387 /156



## Network properties
meta_degree <- sort(igraph::degree(all_control_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "41"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.04850213980029"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_control_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  5.78163545837797"
print(paste("mean clustering coefficient = ", transitivity(all_control_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.270780363778604"
print(paste("mean betweenness centrality = ", mean(betweenness(all_control_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality = 0.00606465960717704"
print(paste("mean closeness centrality = ", mean(closeness(all_control_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.0203369303504658"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_control_net, V(all_control_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors = 4.04850213980029"
##

net <- all_control_net
control_all_deg <- igraph::degree(net,mode="all")
control_all_betweenness <- betweenness(net, normalized = TRUE)
control_all_closeness <- closeness(net, normalized = TRUE)
control_all_transitivity <- transitivity(net, "local", vids = V(net))
names(control_all_transitivity)<- V(net)$name
control_all_transitivity[is.na(control_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
control_all_deg.1percent <- control_all_deg[control_all_deg >= quantile(control_all_deg,prob=1-n/100)]
length(control_all_deg.1percent) #8

control_all_betweenness.1percent <- control_all_betweenness[control_all_betweenness >= quantile(control_all_betweenness,prob=1-n/100)]
length(control_all_betweenness.1percent) #7

control_all_closeness.1percent <- control_all_closeness[control_all_closeness >= quantile(control_all_closeness,prob=1-n/100)]
length(control_all_closeness.1percent) #7

intersect(names(control_all_deg.1percent), names(control_all_betweenness.1percent))
intersect(names(control_all_deg.1percent), names(control_all_closeness.1percent))
intersect(names(control_all_betweenness.1percent), names(control_all_closeness.1percent))
#F87_Penicillium (ITS2 network)
#B17_f_Peptostreptococcaceae (ITS1 network)

### network properties of bacteria and fungi
df.control.degree<-data.frame(control_all_deg)
head(df.control.degree)
df.control.degree$Group <- "control"
names(df.control.degree)[1] <- c("Degree")

df.control.closeness<-data.frame(control_all_closeness)
head(df.control.closeness)
df.control.closeness$Group <- "control"
names(df.control.closeness)[1] <- c("Closeness")

df.control.betweenness<-data.frame(control_all_betweenness)
head(df.control.betweenness)
df.control.betweenness$Group <- "control"
names(df.control.betweenness)[1] <- c("Betweenness")

df.control.degree$kingdom <- ifelse(grepl("^B",rownames(df.control.degree)),'Bacteria', 'Fungi')
df.control.closeness$kingdom <- ifelse(grepl("^B",rownames(df.control.closeness)),'Bacteria', 'Fungi')
df.control.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.control.betweenness)),'Bacteria', 'Fungi')


print_degree <- function(deg){
  
  deg$kingdom <- factor(deg$kingdom, levels=c('Bacteria', 'Fungi'))  
  max.degree <- aggregate(deg$Degree, by = list(deg$kingdom), max)
  
  colnames(max.degree) <- c("kingdom", "maxdegree")
  
  # wilcoxon test
  x <- subset(deg, kingdom=='Bacteria')$Degree
  y <- subset(deg, kingdom=='Fungi')$Degree
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=deg, aes(x=kingdom, y=Degree, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}

p1 <- print_degree(df.control.degree)
p1


## Betweenness centrality
print_betweenness <- function(betw){
  
  betw$kingdom <- factor(betw$kingdom, levels=c('Bacteria','Fungi'))  
  max.betweenness <- aggregate(betw$Betweenness, by = list(betw$kingdom), max)
  
  colnames(max.betweenness) <- c("kingdom", "maxbetweenness")
  
  
  # wilcoxon test
  x <- subset(betw, kingdom=='Bacteria')$Betweenness
  y <- subset(betw, kingdom=='Fungi')$Betweenness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=betw, aes(x=kingdom, y=Betweenness, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Betweenness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p2 <- print_betweenness(df.control.betweenness)
p2

grid.arrange(p1, p2, p3,p4,p5,p6,p7,p8, ncol=4)


## Closeness centrality
print_closeness <- function(closen){
  
  closen$kingdom <- factor(closen$kingdom, levels=c('Bacteria','Fungi'))  
  max.closeness <- aggregate(closen$Closeness, by = list(closen$kingdom), max)
  
  colnames(max.closeness) <- c("kingdom", "maxcloseness")
  
  
  x <- subset(closen, kingdom=='Bacteria')$Closeness
  y <- subset(closen, kingdom=='Fungi')$Closeness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=closen, aes(x=kingdom, y=Closeness, fill=kingdom))+ geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Closeness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p3 <- print_closeness(df.control.closeness)
p3

gridExtra::grid.arrange(p1, p2, p3, ncol=3)

dev.off()



### Plotting hub
write.csv(df.control.net.properties,"Source data for network analysis and topological properties_control.csv")

head(df.control.degree)
head(df.control.betweenness)
df.control.degree$OTU_id <- rownames(df.control.degree)
df.control.betweenness$OTU_id <- rownames(df.control.betweenness)
df.control.closeness$OTU_id <- rownames(df.control.closeness)
df.control.net.properties <- merge(df.control.degree, df.control.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.control.net.properties <- merge(df.control.net.properties, df.control.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.control.net.properties$kingdom <- factor(df.control.net.properties$kingdom, levels = c('Bacteria', 'Fungi'))
n<-1
ggplot(df.control.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Betweenness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
geom_hline(yintercept=quantile(df.control.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.control.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

ggplot(df.control.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Closeness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.control.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.control.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)




### RA subgroup samples (randomly selected) - 30 samples
## read cor and p

RA_30_cor_df_padj <- read.table("2_0.4_edge_RA_30.tsv", sep = '\t', header = T)
RA_30_cor_df_padj <- read.table("2_0.4_edge_RA_its1.tsv", sep = '\t', header = T)

nodeattrib_RA_30_combine <- data.frame(node=union(RA_30_cor_df_padj$Source,RA_30_cor_df_padj$Target))
nodeattrib_RA_30_combine$kingdom <- 0

for (i in as.character(nodeattrib_RA_30_combine$node))
{
  if (i %in% otu.bac.id$OTU_id == TRUE)
  {nodeattrib_RA_30_combine[nodeattrib_RA_30_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_RA_30_combine[nodeattrib_RA_30_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_RA_30_combine) <- as.character(nodeattrib_RA_30_combine$node)
nodeattrib_RA_30_combine$kingdom



all_RA_30_net <- graph_from_data_frame(RA_30_cor_df_padj,direct=F, vertices=nodeattrib_RA_30_combine)

## Number of nodes
length(V(all_RA_30_net)) #801 / 441 (ITS1)

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_RA_30_net)))) #394 #319
length(grep("^F",names(V(all_RA_30_net)))) #407 #122


## Connections 
bb_occur_RA_30 <- droplevels(RA_30_cor_df_padj[with(RA_30_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_RA_30) #421 / 409

ff_occur_RA_30 <- droplevels(RA_30_cor_df_padj[with(RA_30_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_RA_30) #724 / 63

fb_occur_RA_30 <- droplevels(RA_30_cor_df_padj[with(RA_30_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_RA_30) #518 /169

bf_occur_RA_30 <- droplevels(RA_30_cor_df_padj[with(RA_30_cor_df_padj, grepl("^B",Source) & grepl("^F",Target)),])
nrow(bf_occur_RA_30) #0


## Network properties
meta_degree <- sort(igraph::degree(all_RA_30_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "33"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.19225967540574"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_RA_30_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  6.12956849677315"
print(paste("mean clustering coefficient = ", transitivity(all_RA_30_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.289254170755643"
print(paste("mean betweenness centrality = ", mean(betweenness(all_RA_30_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality = 0.00579664186975292"
print(paste("mean closeness centrality = ", mean(closeness(all_RA_30_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.020779265984292"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_RA_30_net, V(all_RA_30_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors = 4.19225967540574"
##

net <- all_RA_30_net
RA_30_all_deg <- igraph::degree(net,mode="all")
RA_30_all_betweenness <- betweenness(net, normalized = TRUE)
RA_30_all_closeness <- closeness(net, normalized = TRUE)
RA_30_all_transitivity <- transitivity(net, "local", vids = V(net))
names(RA_30_all_transitivity)<- V(net)$name
RA_30_all_transitivity[is.na(RA_30_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
RA_30_all_deg.1percent <- RA_30_all_deg[RA_30_all_deg >= quantile(RA_30_all_deg,prob=1-n/100)]
length(RA_30_all_deg.1percent) #9

RA_30_all_betweenness.1percent <- RA_30_all_betweenness[RA_30_all_betweenness >= quantile(RA_30_all_betweenness,prob=1-n/100)]
length(RA_30_all_betweenness.1percent) #8

RA_30_all_closeness.1percent <- RA_30_all_closeness[RA_30_all_closeness >= quantile(RA_30_all_closeness,prob=1-n/100)]
length(RA_30_all_closeness.1percent) #8

intersect(names(RA_30_all_deg.1percent), names(RA_30_all_betweenness.1percent))
intersect(names(RA_30_all_deg.1percent), names(RA_30_all_closeness.1percent))
intersect(names(RA_30_all_betweenness.1percent), names(RA_30_all_closeness.1percent))
#B3_f_Lachnospiraceae (ITS2)
#B3_f_Lachnospiraceae (ITS1)

### network properties of bacteria and fungi
df.RA_30.degree<-data.frame(RA_30_all_deg)
head(df.RA_30.degree)
df.RA_30.degree$Group <- "RA"
names(df.RA_30.degree)[1] <- c("Degree")

df.RA_30.closeness<-data.frame(RA_30_all_closeness)
head(df.RA_30.closeness)
df.RA_30.closeness$Group <- "RA"
names(df.RA_30.closeness)[1] <- c("Closeness")

df.RA_30.betweenness<-data.frame(RA_30_all_betweenness)
head(df.RA_30.betweenness)
df.RA_30.betweenness$Group <- "RA"
names(df.RA_30.betweenness)[1] <- c("Betweenness")

df.RA_30.degree$kingdom <- ifelse(grepl("^B",rownames(df.RA_30.degree)),'Bacteria', 'Fungi')
df.RA_30.closeness$kingdom <- ifelse(grepl("^B",rownames(df.RA_30.closeness)),'Bacteria', 'Fungi')
df.RA_30.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.RA_30.betweenness)),'Bacteria', 'Fungi')



### Plotting hub
df.RA_30.degree$OTU_id <- rownames(df.RA_30.degree)
df.RA_30.betweenness$OTU_id <- rownames(df.RA_30.betweenness)
df.RA_30.closeness$OTU_id <- rownames(df.RA_30.closeness)
df.RA_30.net.properties <- merge(df.RA_30.degree, df.RA_30.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.RA_30.net.properties <- merge(df.RA_30.net.properties, df.RA_30.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))

n<-1
ggplot(df.RA_30.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Betweenness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.RA_30.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.RA_30.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

ggplot(df.RA_30.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Closeness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.RA_30.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.RA_30.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

## kingdom comparison
print_closeness(df.RA_30.closeness)
print_betweenness(df.RA_30.betweenness)
print_degree(df.RA_30.degree)


# ####RA20
# ### RA subgroup samples (randomly selected)
# ## read cor and p
# 
# RA_20_cor_df_padj <- read.table("2_0.4_edge_RA_20.tsv", sep = '\t', header = T)
# 
# 
# nodeattrib_RA_20_combine <- data.frame(node=union(RA_20_cor_df_padj$Source,RA_20_cor_df_padj$Target))
# nodeattrib_RA_20_combine$kingdom <- 0
# 
# for (i in as.character(nodeattrib_RA_20_combine$node))
# {
#   if (i %in% otu.bac.id$OTU_id == TRUE)
#   {nodeattrib_RA_20_combine[nodeattrib_RA_20_combine$node==i,"kingdom"] <- "Bacteria"}
#   
#   else
#   {nodeattrib_RA_20_combine[nodeattrib_RA_20_combine$node==i,"kingdom"]<- "Fungi"}
# }
# 
# rownames(nodeattrib_RA_20_combine) <- as.character(nodeattrib_RA_20_combine$node)
# nodeattrib_RA_20_combine$kingdom
# 
# 
# 
# all_RA_20_net <- graph_from_data_frame(RA_20_cor_df_padj,direct=F, vertices=nodeattrib_RA_20_combine)
# 
# ## Number of nodes
# length(V(all_RA_20_net)) #801
# 
# ## Number of bacteria and fungi nodes
# length(grep("^B",names(V(all_RA_20_net)))) #401
# length(grep("^F",names(V(all_RA_20_net)))) #400
# 
# 
# ## Connections 
# bb_occur_RA_20 <- droplevels(RA_20_cor_df_padj[with(RA_20_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
# nrow(bb_occur_RA_20) #441
# 
# ff_occur_RA_20 <- droplevels(RA_20_cor_df_padj[with(RA_20_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
# nrow(ff_occur_RA_20) #716
# 
# fb_occur_RA_20 <- droplevels(RA_20_cor_df_padj[with(RA_20_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
# nrow(fb_occur_RA_20) #522
# 
# bf_occur_RA_20 <- droplevels(RA_20_cor_df_padj[with(RA_20_cor_df_padj, grepl("^B",Source) & grepl("^F",Target)),])
# nrow(bf_occur_RA_20) #0
# 
# 
# ## Network properties
# meta_degree <- sort(igraph::degree(all_RA_20_net,mode="all"),decr=T)
# print(c('max degree = ', max(meta_degree))) #"max degree = " "33"     
# print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.19225967540574"
# meta_degree <- as.data.frame(meta_degree)
# ggplot(meta_degree, aes(x=meta_degree)) + geom_density()
# 
# print(paste("average shortest path length = ", mean_distance(all_RA_20_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
# #"average shortest path length =  6.12956849677315"
# print(paste("mean clustering coefficient = ", transitivity(all_RA_20_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
# #"mean clustering coefficient =  0.289254170755643"
# print(paste("mean betweenness centrality = ", mean(betweenness(all_RA_20_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
# #"mean betweenness centrality = 0.00579664186975292"
# print(paste("mean closeness centrality = ", mean(closeness(all_RA_20_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
# #"mean closeness centrality =  0.020779265984292"
# print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_RA_20_net, V(all_RA_20_net)))))) #"mean number of neighbors =  223.182421227197"
# #"mean number of neighbors = 4.19225967540574"
# ##
# 
# net <- all_RA_20_net
# RA_20_all_deg <- igraph::degree(net,mode="all")
# RA_20_all_betweenness <- betweenness(net, normalized = TRUE)
# RA_20_all_closeness <- closeness(net, normalized = TRUE)
# RA_20_all_transitivity <- transitivity(net, "local", vids = V(net))
# names(RA_20_all_transitivity)<- V(net)$name
# RA_20_all_transitivity[is.na(RA_20_all_transitivity)] <- 0
# 
# 
# ## Defining hub OTUs
# n <- 1
# RA_20_all_deg.1percent <- RA_20_all_deg[RA_20_all_deg >= quantile(RA_20_all_deg,prob=1-n/100)]
# length(RA_20_all_deg.1percent) #9
# 
# RA_20_all_betweenness.1percent <- RA_20_all_betweenness[RA_20_all_betweenness >= quantile(RA_20_all_betweenness,prob=1-n/100)]
# length(RA_20_all_betweenness.1percent) #8
# 
# RA_20_all_closeness.1percent <- RA_20_all_closeness[RA_20_all_closeness >= quantile(RA_20_all_closeness,prob=1-n/100)]
# length(RA_20_all_closeness.1percent) #8
# 
# intersect(names(RA_20_all_deg.1percent), names(RA_20_all_betweenness.1percent))
# intersect(names(RA_20_all_deg.1percent), names(RA_20_all_closeness.1percent))
# intersect(names(RA_20_all_betweenness.1percent), names(RA_20_all_closeness.1percent))
# 
# ### network properties of bacteria and fungi
# df.RA_20.degree<-data.frame(RA_20_all_deg)
# head(df.RA_20.degree)
# df.RA_20.degree$Group <- "RA"
# names(df.RA_20.degree)[1] <- c("Degree")
# 
# df.RA_20.closeness<-data.frame(RA_20_all_closeness)
# head(df.RA_20.closeness)
# df.RA_20.closeness$Group <- "RA"
# names(df.RA_20.closeness)[1] <- c("Closeness")
# 
# df.RA_20.betweenness<-data.frame(RA_20_all_betweenness)
# head(df.RA_20.betweenness)
# df.RA_20.betweenness$Group <- "RA"
# names(df.RA_20.betweenness)[1] <- c("Betweenness")
# 
# df.RA_20.degree$kingdom <- ifelse(grepl("^B",rownames(df.RA_20.degree)),'Bacteria', 'Fungi')
# df.RA_20.closeness$kingdom <- ifelse(grepl("^B",rownames(df.RA_20.closeness)),'Bacteria', 'Fungi')
# df.RA_20.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.RA_20.betweenness)),'Bacteria', 'Fungi')
# 
# 
# 
# ### Plotting hub
# df.RA_20.degree$OTU_id <- rownames(df.RA_20.degree)
# df.RA_20.betweenness$OTU_id <- rownames(df.RA_20.betweenness)
# df.RA_20.closeness$OTU_id <- rownames(df.RA_20.closeness)
# df.RA_20.net.properties <- merge(df.RA_20.degree, df.RA_20.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
# df.RA_20.net.properties <- merge(df.RA_20.net.properties, df.RA_20.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
# 
# n<-1
# ggplot(df.RA_20.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
#   xlab('\n Degree')+
#   ylab("Betweenness centrality \n") +
#   geom_point(size=3, alpha=0.7) +
#   theme(aspect.ratio = 1)+
#   # ggtitle("Volcano Plot \n") +
#   theme(legend.text=element_text(size=13)) + 
#   theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
#   theme(legend.position="top") +
#   theme(legend.title=element_blank()) +
#   guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
#   guides(size=FALSE) +
#   #scale_x_continuous(breaks=seq(0,1,0.2))+
#   #scale_y_continuous(breaks=seq(-20,0,-5))+
#   scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
#   theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
#   theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
#   theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
#   theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
#   theme(panel.grid.major = element_blank())+
#   theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
#   geom_hline(yintercept=quantile(df.RA_20.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.RA_20.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)
# 
# ggplot(df.RA_20.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
#   xlab('\n Degree')+
#   ylab("Closeness centrality \n") +
#   geom_point(size=3, alpha=0.7) +
#   theme(aspect.ratio = 1)+
#   # ggtitle("Volcano Plot \n") +
#   theme(legend.text=element_text(size=13)) + 
#   theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
#   theme(legend.position="top") +
#   theme(legend.title=element_blank()) +
#   guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
#   guides(size=FALSE) +
#   #scale_x_continuous(breaks=seq(0,1,0.2))+
#   #scale_y_continuous(breaks=seq(-20,0,-5))+
#   scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
#   theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
#   theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
#   theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
#   theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
#   theme(panel.grid.major = element_blank())+
#   theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
#   geom_hline(yintercept=quantile(df.RA_20.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.RA_20.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)
# 
# ## kingdom comparison
# print_closeness(df.RA_20.closeness)
# print_betweenness(df.RA_20.betweenness)
# print_degree(df.RA_20.degree)

# ### Control 20
# ####RA20
# ### RA subgroup samples (randomly selected)
# ## read cor and p
# 
# control_20_cor_df_padj <- read.table("2_0.4_edge_control_20.tsv", sep = '\t', header = T)
# 
# 
# nodeattrib_control_20_combine <- data.frame(node=union(control_20_cor_df_padj$Source,control_20_cor_df_padj$Target))
# nodeattrib_control_20_combine$kingdom <- 0
# 
# for (i in as.character(nodeattrib_control_20_combine$node))
# {
#   if (i %in% otu.bac.id$OTU_id == TRUE)
#   {nodeattrib_control_20_combine[nodeattrib_control_20_combine$node==i,"kingdom"] <- "Bacteria"}
#   
#   else
#   {nodeattrib_control_20_combine[nodeattrib_control_20_combine$node==i,"kingdom"]<- "Fungi"}
# }
# 
# rownames(nodeattrib_control_20_combine) <- as.character(nodeattrib_control_20_combine$node)
# nodeattrib_control_20_combine$kingdom
# 
# 
# 
# all_control_20_net <- graph_from_data_frame(control_20_cor_df_padj,direct=F, vertices=nodeattrib_control_20_combine)
# 
# ## Number of nodes
# length(V(all_control_20_net)) #926
# 
# ## Number of bacteria and fungi nodes
# length(grep("^B",names(V(all_control_20_net)))) #530
# length(grep("^F",names(V(all_control_20_net)))) #396
# 
# 
# ## Connections 
# bb_occur_control_20 <- droplevels(control_20_cor_df_padj[with(control_20_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
# nrow(bb_occur_control_20) #441
# 
# ff_occur_control_20 <- droplevels(control_20_cor_df_padj[with(control_20_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
# nrow(ff_occur_control_20) #716
# 
# fb_occur_control_20 <- droplevels(control_20_cor_df_padj[with(control_20_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
# nrow(fb_occur_control_20) #522
# 
# bf_occur_control_20 <- droplevels(control_20_cor_df_padj[with(control_20_cor_df_padj, grepl("^B",Source) & grepl("^F",Target)),])
# nrow(bf_occur_control_20) #0
# 
# 
# ## Network properties
# meta_degree <- sort(igraph::degree(all_control_20_net,mode="all"),decr=T)
# print(c('max degree = ', max(meta_degree))) #"max degree = " "33"     
# print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.19225967540574"
# meta_degree <- as.data.frame(meta_degree)
# ggplot(meta_degree, aes(x=meta_degree)) + geom_density()
# 
# print(paste("average shortest path length = ", mean_distance(all_control_20_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
# #"average shortest path length =  6.12956849677315"
# print(paste("mean clustering coefficient = ", transitivity(all_control_20_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
# #"mean clustering coefficient =  0.289254170755643"
# print(paste("mean betweenness centrality = ", mean(betweenness(all_control_20_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
# #"mean betweenness centrality = 0.00579664186975292"
# print(paste("mean closeness centrality = ", mean(closeness(all_control_20_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
# #"mean closeness centrality =  0.020779265984292"
# print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_control_20_net, V(all_control_20_net)))))) #"mean number of neighbors =  223.182421227197"
# #"mean number of neighbors = 4.19225967540574"
# ##
# 
# net <- all_control_20_net
# control_20_all_deg <- igraph::degree(net,mode="all")
# control_20_all_betweenness <- betweenness(net, normalized = TRUE)
# control_20_all_closeness <- closeness(net, normalized = TRUE)
# control_20_all_transitivity <- transitivity(net, "local", vids = V(net))
# names(control_20_all_transitivity)<- V(net)$name
# control_20_all_transitivity[is.na(control_20_all_transitivity)] <- 0
# 
# 
# ## Defining hub OTUs
# n <- 1
# control_20_all_deg.1percent <- control_20_all_deg[control_20_all_deg >= quantile(control_20_all_deg,prob=1-n/100)]
# length(control_20_all_deg.1percent) #9
# 
# control_20_all_betweenness.1percent <- control_20_all_betweenness[control_20_all_betweenness >= quantile(control_20_all_betweenness,prob=1-n/100)]
# length(control_20_all_betweenness.1percent) #8
# 
# control_20_all_closeness.1percent <- control_20_all_closeness[control_20_all_closeness >= quantile(control_20_all_closeness,prob=1-n/100)]
# length(control_20_all_closeness.1percent) #8
# 
# intersect(names(control_20_all_deg.1percent), names(control_20_all_betweenness.1percent))
# intersect(names(control_20_all_deg.1percent), names(control_20_all_closeness.1percent))
# intersect(names(control_20_all_betweenness.1percent), names(control_20_all_closeness.1percent))
# 
# ### network properties of bacteria and fungi
# df.control_20.degree<-data.frame(control_20_all_deg)
# head(df.control_20.degree)
# df.control_20.degree$Group <- "RA"
# names(df.control_20.degree)[1] <- c("Degree")
# 
# df.control_20.closeness<-data.frame(control_20_all_closeness)
# head(df.control_20.closeness)
# df.control_20.closeness$Group <- "RA"
# names(df.control_20.closeness)[1] <- c("Closeness")
# 
# df.control_20.betweenness<-data.frame(control_20_all_betweenness)
# head(df.control_20.betweenness)
# df.control_20.betweenness$Group <- "RA"
# names(df.control_20.betweenness)[1] <- c("Betweenness")
# 
# df.control_20.degree$kingdom <- ifelse(grepl("^B",rownames(df.control_20.degree)),'Bacteria', 'Fungi')
# df.control_20.closeness$kingdom <- ifelse(grepl("^B",rownames(df.control_20.closeness)),'Bacteria', 'Fungi')
# df.control_20.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.control_20.betweenness)),'Bacteria', 'Fungi')
# 
# 
# 
# ### Plotting hub
# df.control_20.degree$OTU_id <- rownames(df.control_20.degree)
# df.control_20.betweenness$OTU_id <- rownames(df.control_20.betweenness)
# df.control_20.closeness$OTU_id <- rownames(df.control_20.closeness)
# df.control_20.net.properties <- merge(df.control_20.degree, df.control_20.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
# df.control_20.net.properties <- merge(df.control_20.net.properties, df.control_20.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
# 
# n<-1
# ggplot(df.control_20.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
#   xlab('\n Degree')+
#   ylab("Betweenness centrality \n") +
#   geom_point(size=3, alpha=0.7) +
#   theme(aspect.ratio = 1)+
#   # ggtitle("Volcano Plot \n") +
#   theme(legend.text=element_text(size=13)) + 
#   theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
#   theme(legend.position="top") +
#   theme(legend.title=element_blank()) +
#   guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
#   guides(size=FALSE) +
#   #scale_x_continuous(breaks=seq(0,1,0.2))+
#   #scale_y_continuous(breaks=seq(-20,0,-5))+
#   scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
#   theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
#   theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
#   theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
#   theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
#   theme(panel.grid.major = element_blank())+
#   theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
#   geom_hline(yintercept=quantile(df.control_20.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.control_20.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)
# 
# ggplot(df.control_20.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
#   xlab('\n Degree')+
#   ylab("Closeness centrality \n") +
#   geom_point(size=3, alpha=0.7) +
#   theme(aspect.ratio = 1)+
#   # ggtitle("Volcano Plot \n") +
#   theme(legend.text=element_text(size=13)) + 
#   theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
#   theme(legend.position="top") +
#   theme(legend.title=element_blank()) +
#   guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
#   guides(size=FALSE) +
#   #scale_x_continuous(breaks=seq(0,1,0.2))+
#   #scale_y_continuous(breaks=seq(-20,0,-5))+
#   scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
#   theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
#   theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
#   theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
#   theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
#   theme(panel.grid.major = element_blank())+
#   theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
#   geom_hline(yintercept=quantile(df.control_20.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.control_20.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)
# 
# ## kingdom comparison
# print_closeness(df.control_20.closeness)
# print_betweenness(df.control_20.betweenness)
# print_degree(df.control_20.degree)
# 

### Control vs RA
df.control.net.properties
df.RA_30.net.properties
head(df.net.properties)

df.net.properties <- rbind(df.control.net.properties,df.RA_30.net.properties)

write.csv(df.net.properties, "Source data for network analysis_control and RA.csv")

print_degree_diagnosis <- function(deg){
  
  deg$Group <- factor(deg$Group, levels=c('control', 'RA'))  
  max.degree <- aggregate(deg$Degree, by = list(deg$Group), max)
  
  colnames(max.degree) <- c("Group", "maxdegree")
  
  # wilcoxon test
  x <- subset(deg, Group=='control')$Degree
  y <- subset(deg, Group=='RA')$Degree
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=deg, aes(x=Group, y=Degree, fill=Group)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
print_betweenness_diagnosis <- function(deg){
  
  deg$Group <- factor(deg$Group, levels=c('control', 'RA'))  
  max.degree <- aggregate(deg$Betweenness, by = list(deg$Group), max)
  
  colnames(max.degree) <- c("Group", "maxdegree")
  
  # wilcoxon test
  x <- subset(deg, Group=='control')$Betweenness
  y <- subset(deg, Group=='RA')$Betweenness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=deg, aes(x=Group, y=Betweenness, fill=Group)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Betweenness \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
print_closeness_diagnosis <- function(deg){
  
  deg$Group <- factor(deg$Group, levels=c('control', 'RA'))  
  max.degree <- aggregate(deg$Closeness, by = list(deg$Group), max)
  
  colnames(max.degree) <- c("Group", "maxdegree")
  
  # wilcoxon test
  x <- subset(deg, Group=='control')$Closeness
  y <- subset(deg, Group=='RA')$Closeness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=deg, aes(x=Group, y=Closeness, fill=Group)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Closeness \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}

print_degree_diagnosis(df.net.properties)
print_betweenness_diagnosis(df.net.properties)
print_closeness_diagnosis(df.net.properties)

write.xlsx(df.net.properties, "Node properties - ITS2 based network.xlsx")
write.xlsx(df.net.properties, "Node properties - ITS1 based network.xlsx")
## Comparison kingdom and diagnosis
df.net.properties
df.net.properties.bac <- subset(df.net.properties, kingdom == "Bacteria")
df.net.properties.fun <- subset(df.net.properties, kingdom == "Fungi")


print_degree_diagnosis(df.net.properties.bac)
print_betweenness_diagnosis(df.net.properties.bac)
print_closeness_diagnosis(df.net.properties.bac)


print_degree_diagnosis(df.net.properties.fun)
print_betweenness_diagnosis(df.net.properties.fun)
print_closeness_diagnosis(df.net.properties.fun)



#### Candida
# df.net.properties
# 
# RA_30_cor_df_padj$Target[RA_30_cor_df_padj$Source == "F1_Candida"]
# #B4_Bifidobacterium               B30_.Ruminococcus..torques.group B74_Dorea
# RA_30_cor_df_padj$Cor[RA_30_cor_df_padj$Source == "F1_Candida"]
# # -0.4202 -0.4236 -0.4122
# 
# RA_30_cor_df_padj$Source[RA_30_cor_df_padj$Target == "B4_Bifidobacterium"]
# #B4_Bifidobacterium               B30_.Ruminococcus..torques.group B74_Dorea
# 
# RA_30_cor_df_padj$Source[RA_30_cor_df_padj$Target == "B30_.Ruminococcus..torques.group"]
# #B52_Collinsella                  B30_[Ruminococcus] torques group B2_Bifidobacterium
# RA_30_cor_df_padj$Source[RA_30_cor_df_padj$Target == "B2_Bifidobacterium"]
# #F10_Rhodotorula
# RA_30_cor_df_padj$Cor[RA_30_cor_df_padj$Target == "B2_Bifidobacterium"]
# 
# RA_30_cor_df_padj$Target[RA_30_cor_df_padj$Source == "B2_Bifidobacterium"]
# RA_30_cor_df_padj$Cor[RA_30_cor_df_padj$Source == "B2_Bifidobacterium"]
# RA_30_cor_df_padj$Cor[RA_30_cor_df_padj$Source == "F1_Candida"]
# # -0.4202 -0.4236 -0.4122
# 
# control_cor_df_padj$Target[control_cor_df_padj$Source == "F1_Candida"]
# #B13_Blautia B39_Blautia B36_Blautia B38_Blautia
# control_cor_df_padj$Cor[control_cor_df_padj$Source == "F1_Candida"]
# #0.4403 -0.4102  0.4427  0.4662


### All RA samples
## read cor and p
cor_RA <- read.table(file = 'RA_median_correlation.tsv', sep = '\t', header = TRUE)
rownames(cor_RA) <- cor_RA$OTU_id
cor_RA <- cor_RA[-c(1)]

p_RA <- read.table(file = 'RA_pvalues.tsv', sep = '\t', header = TRUE)
rownames(p_RA) <- p_RA$OTU_id
p_RA <- p_RA[-c(1)]

source("CorrDF.r")

RA_cor_df <- CorrDF(cor_RA, p_RA)
RA_cor_df$padj <- p.adjust(RA_cor_df$p, method="none")
RA_cor_df_padj <- RA_cor_df[which(abs(RA_cor_df$cor) >= 0.4),]
RA_cor_df_padj <- RA_cor_df_padj[which(RA_cor_df_padj$padj < 0.05),]


nodeattrib_RA_combine <- data.frame(node=union(RA_cor_df_padj$from,RA_cor_df_padj$to))
nodeattrib_RA_combine$kingdom <- 0

for (i in as.character(nodeattrib_RA_combine$node))
{
  if (i %in% otu.bac.id$OTU_id == TRUE)
  {nodeattrib_RA_combine[nodeattrib_RA_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_RA_combine[nodeattrib_RA_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_RA_combine) <- as.character(nodeattrib_RA_combine$node)
nodeattrib_RA_combine$kingdom



all_RA_net <- graph_from_data_frame(RA_cor_df_padj,direct=F, vertices=nodeattrib_RA_combine)

## Number of nodes
length(V(all_RA_net)) #128

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_RA_net)))) #91
length(grep("^F",names(V(all_RA_net)))) #37


## Connections 
bb_occur_RA <- droplevels(RA_cor_df_padj[with(RA_cor_df_padj, grepl("^B",from) & grepl("^B",to)),])
nrow(bb_occur_RA) #142

ff_occur_RA <- droplevels(RA_cor_df_padj[with(RA_cor_df_padj, grepl("^F",from) & grepl("^F",to)),])
nrow(ff_occur_RA) #27

fb_occur_RA <- droplevels(RA_cor_df_padj[with(RA_cor_df_padj, grepl("^F",from) & grepl("^B",to)),])
nrow(fb_occur_RA) #0

bf_occur_RA <- droplevels(RA_cor_df_padj[with(RA_cor_df_padj, grepl("^B",from) & grepl("^F",to)),])
nrow(bf_occur_RA) #0


## Network properties
meta_degree <- sort(igraph::degree(all_RA_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "8"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "2.640625"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_RA_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  1.22685185185185"
print(paste("mean clustering coefficient = ", transitivity(all_RA_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.883928571428571"
print(paste("mean betweenness centrality = ", mean(betweenness(all_RA_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality = 4.78455818022747e-05"
print(paste("mean closeness centrality = ", mean(closeness(all_RA_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.00802599418391875"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_RA_net, V(all_RA_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors = 2.640625"
##

net <- all_RA_net
RA_all_deg <- igraph::degree(net,mode="all")
RA_all_betweenness <- betweenness(net, normalized = TRUE)
RA_all_closeness <- closeness(net, normalized = TRUE)
RA_all_transitivity <- transitivity(net, "local", vids = V(net))
names(RA_all_transitivity)<- V(net)$name
RA_all_transitivity[is.na(RA_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
RA_all_deg.1percent <- RA_all_deg[RA_all_deg >= quantile(RA_all_deg,prob=1-n/100)]
length(RA_all_deg.1percent) #2

RA_all_betweenness.1percent <- RA_all_betweenness[RA_all_betweenness >= quantile(RA_all_betweenness,prob=1-n/100)]
length(RA_all_betweenness.1percent) #2

RA_all_closeness.1percent <- RA_all_closeness[RA_all_closeness >= quantile(RA_all_closeness,prob=1-n/100)]
length(RA_all_closeness.1percent) #2

intersect(names(RA_all_deg.1percent), names(RA_all_betweenness.1percent))
intersect(names(RA_all_deg.1percent), names(RA_all_closeness.1percent))
intersect(names(RA_all_betweenness.1percent), names(RA_all_closeness.1percent))

### network properties of bacteria and fungi
df.RA.degree<-data.frame(RA_all_deg)
head(df.RA.degree)
df.RA.degree$Group <- "RA"
names(df.RA.degree)[1] <- c("Degree")

df.RA.closeness<-data.frame(RA_all_closeness)
head(df.RA.closeness)
df.RA.closeness$Group <- "RA"
names(df.RA.closeness)[1] <- c("Closeness")

df.RA.betweenness<-data.frame(RA_all_betweenness)
head(df.RA.betweenness)
df.RA.betweenness$Group <- "RA"
names(df.RA.betweenness)[1] <- c("Betweenness")

df.RA.degree$kingdom <- ifelse(grepl("^B",rownames(df.RA.degree)),'Bacteria', 'Fungi')
df.RA.closeness$kingdom <- ifelse(grepl("^B",rownames(df.RA.closeness)),'Bacteria', 'Fungi')
df.RA.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.RA.betweenness)),'Bacteria', 'Fungi')



### Plotting hub
df.RA.degree$OTU_id <- rownames(df.RA.degree)
df.RA.betweenness$OTU_id <- rownames(df.RA.betweenness)
df.RA.closeness$OTU_id <- rownames(df.RA.closeness)
df.RA.net.properties <- merge(df.RA.degree, df.RA.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.RA.net.properties <- merge(df.RA.net.properties, df.RA.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))

n<-1
ggplot(df.RA.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Betweenness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.RA.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.RA.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

ggplot(df.RA.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Closeness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.RA.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.RA.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

## kingdom comparison
print_closeness(df.RA.closeness)
print_betweenness(df.RA.betweenness)
print_degree(df.RA.degree)


