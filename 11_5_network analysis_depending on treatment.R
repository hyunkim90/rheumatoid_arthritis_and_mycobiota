#### Network properties of network (20 samples)
### Network properties
## read cor and p

control_cor_df_padj <- read.table("2_0.4_edge_control.tsv", sep='\t', header =T)


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
length(V(all_control_net)) #579

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_control_net)))) #299
length(grep("^F",names(V(all_control_net)))) #280


## Connections 
bb_occur_control <- droplevels(control_cor_df_padj[with(control_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_control) #301

ff_occur_control <- droplevels(control_cor_df_padj[with(control_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_control) #472

fb_occur_control <- droplevels(control_cor_df_padj[with(control_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_control) #331



## Network properties
meta_degree <- sort(igraph::degree(all_control_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "41"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "3.81347150259067"
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

intersect(names(control_all_deg.1percent), names(control_all_betweenness.1percent)) #"B9_f_Lachnospiraceae" "F167_Hannaella"
intersect(names(control_all_deg.1percent), names(control_all_closeness.1percent)) #"F167_Hannaella"
intersect(names(control_all_betweenness.1percent), names(control_all_closeness.1percent))# "F63_Penicillium" "F167_Hannaella"  "F15_Candida"    

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
  scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
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
  scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.control.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.control.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)



### RA not treatedd 20 samples
## read cor and p

RA_cor_df_padj <- read.table("3_0.5_edge_RA_20.tsv", sep='\t', header =T)


nodeattrib_RA_combine <- data.frame(node=union(RA_cor_df_padj$Source,RA_cor_df_padj$Target))
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
length(V(all_RA_net)) #684

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_RA_net)))) #366
length(grep("^F",names(V(all_RA_net)))) #318


## Connections 
bb_occur_RA <- droplevels(RA_cor_df_padj[with(RA_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_RA) #487

ff_occur_RA <- droplevels(RA_cor_df_padj[with(RA_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_RA) #687

fb_occur_RA <- droplevels(RA_cor_df_padj[with(RA_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_RA) #567

bf_occur_RA <- droplevels(RA_cor_df_padj[with(RA_cor_df_padj, grepl("^B",Source) & grepl("^F",Target)),])
nrow(bf_occur_RA) #0


## Network properties
meta_degree <- sort(igraph::degree(all_RA_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "32"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "5.0906432748538"
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
length(RA_all_deg.1percent) #11

RA_all_betweenness.1percent <- RA_all_betweenness[RA_all_betweenness >= quantile(RA_all_betweenness,prob=1-n/100)]
length(RA_all_betweenness.1percent) #7

RA_all_closeness.1percent <- RA_all_closeness[RA_all_closeness >= quantile(RA_all_closeness,prob=1-n/100)]
length(RA_all_closeness.1percent) #7

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



### TNF inhibitor-treated subjects
## read cor and p

TNF_cor_df_padj <- read.table("3_0.5_edge_TNF.tsv", sep = '\t', header = T)


nodeattrib_TNF_combine <- data.frame(node=union(TNF_cor_df_padj$Source,TNF_cor_df_padj$Target))
nodeattrib_TNF_combine$kingdom <- 0

for (i in as.character(nodeattrib_TNF_combine$node))
{
  if (i %in% otu.bac.id$OTU_id == TRUE)
  {nodeattrib_TNF_combine[nodeattrib_TNF_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_TNF_combine[nodeattrib_TNF_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_TNF_combine) <- as.character(nodeattrib_TNF_combine$node)
nodeattrib_TNF_combine$kingdom



all_TNF_net <- graph_from_data_frame(TNF_cor_df_padj,direct=F, vertices=nodeattrib_TNF_combine)

## Number of nodes
length(V(all_TNF_net)) #697

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_TNF_net)))) #341
length(grep("^F",names(V(all_TNF_net)))) #356


## Connections 
bb_occur_TNF <- droplevels(TNF_cor_df_padj[with(TNF_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_TNF) #441

ff_occur_TNF <- droplevels(TNF_cor_df_padj[with(TNF_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_TNF) #716

fb_occur_TNF <- droplevels(TNF_cor_df_padj[with(TNF_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_TNF) #522

bf_occur_TNF <- droplevels(TNF_cor_df_padj[with(TNF_cor_df_padj, grepl("^B",Source) & grepl("^F",Target)),])
nrow(bf_occur_TNF) #0


## Network properties
meta_degree <- sort(igraph::degree(all_TNF_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "33"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.19225967540574"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_TNF_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  6.12956849677315"
print(paste("mean clustering coefficient = ", transitivity(all_TNF_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.289254170755643"
print(paste("mean betweenness centrality = ", mean(betweenness(all_TNF_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality = 0.00579664186975292"
print(paste("mean closeness centrality = ", mean(closeness(all_TNF_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.020779265984292"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_TNF_net, V(all_TNF_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors = 4.19225967540574"
##

net <- all_TNF_net
TNF_all_deg <- igraph::degree(net,mode="all")
TNF_all_betweenness <- betweenness(net, normalized = TRUE)
TNF_all_closeness <- closeness(net, normalized = TRUE)
TNF_all_transitivity <- transitivity(net, "local", vids = V(net))
names(TNF_all_transitivity)<- V(net)$name
TNF_all_transitivity[is.na(TNF_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
TNF_all_deg.1percent <- TNF_all_deg[TNF_all_deg >= quantile(TNF_all_deg,prob=1-n/100)]
length(TNF_all_deg.1percent) #9

TNF_all_betweenness.1percent <- TNF_all_betweenness[TNF_all_betweenness >= quantile(TNF_all_betweenness,prob=1-n/100)]
length(TNF_all_betweenness.1percent) #8

TNF_all_closeness.1percent <- TNF_all_closeness[TNF_all_closeness >= quantile(TNF_all_closeness,prob=1-n/100)]
length(TNF_all_closeness.1percent) #8

intersect(names(TNF_all_deg.1percent), names(TNF_all_betweenness.1percent)) #"B7_f_Lachnospiraceae"  "F70_o_Chaetothyriales" "F168_f_Didymellaceae" 
intersect(names(TNF_all_deg.1percent), names(TNF_all_closeness.1percent)) #"B7_f_Lachnospiraceae"
intersect(names(TNF_all_betweenness.1percent), names(TNF_all_closeness.1percent))#"B7_f_Lachnospiraceae" "B9_f_Lachnospiraceae" "B68_Faecalibacterium" "F43_Penicillium"

### network properties of bacteria and fungi
df.TNF.degree<-data.frame(TNF_all_deg)
head(df.TNF.degree)
df.TNF.degree$Group <- "RA"
names(df.TNF.degree)[1] <- c("Degree")

df.TNF.closeness<-data.frame(TNF_all_closeness)
head(df.TNF.closeness)
df.TNF.closeness$Group <- "RA"
names(df.TNF.closeness)[1] <- c("Closeness")

df.TNF.betweenness<-data.frame(TNF_all_betweenness)
head(df.TNF.betweenness)
df.TNF.betweenness$Group <- "RA"
names(df.TNF.betweenness)[1] <- c("Betweenness")

df.TNF.degree$kingdom <- ifelse(grepl("^B",rownames(df.TNF.degree)),'Bacteria', 'Fungi')
df.TNF.closeness$kingdom <- ifelse(grepl("^B",rownames(df.TNF.closeness)),'Bacteria', 'Fungi')
df.TNF.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.TNF.betweenness)),'Bacteria', 'Fungi')



### Plotting hub
df.TNF.degree$OTU_id <- rownames(df.TNF.degree)
df.TNF.betweenness$OTU_id <- rownames(df.TNF.betweenness)
df.TNF.closeness$OTU_id <- rownames(df.TNF.closeness)
df.TNF.net.properties <- merge(df.TNF.degree, df.TNF.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.TNF.net.properties <- merge(df.TNF.net.properties, df.TNF.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))

n<-1
ggplot(df.TNF.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
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
  geom_hline(yintercept=quantile(df.TNF.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.TNF.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

ggplot(df.TNF.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
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
  geom_hline(yintercept=quantile(df.TNF.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.TNF.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

## kingdom comparison
print_closeness(df.TNF.closeness)
print_betweenness(df.TNF.betweenness)
print_degree(df.TNF.degree)


####RA20
### RA subgroup samples (randomly selected)
## read cor and p

MTX_cor_df_padj <- read.table("3_0.5_edge_MTX.tsv", sep = '\t', header = T)


nodeattrib_MTX_combine <- data.frame(node=union(MTX_cor_df_padj$Source,MTX_cor_df_padj$Target))
nodeattrib_MTX_combine$kingdom <- 0

for (i in as.character(nodeattrib_MTX_combine$node))
{
  if (i %in% otu.bac.id$OTU_id == TRUE)
  {nodeattrib_MTX_combine[nodeattrib_MTX_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_MTX_combine[nodeattrib_MTX_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_MTX_combine) <- as.character(nodeattrib_MTX_combine$node)
nodeattrib_MTX_combine$kingdom



all_MTX_net <- graph_from_data_frame(MTX_cor_df_padj,direct=F, vertices=nodeattrib_MTX_combine)

## Number of nodes
length(V(all_MTX_net)) #720

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_MTX_net)))) #410
length(grep("^F",names(V(all_MTX_net)))) #310


## Connections 
bb_occur_MTX <- droplevels(MTX_cor_df_padj[with(MTX_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_MTX) #364

ff_occur_MTX <- droplevels(MTX_cor_df_padj[with(MTX_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_MTX) #600

fb_occur_MTX <- droplevels(MTX_cor_df_padj[with(MTX_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_MTX) #521

bf_occur_MTX <- droplevels(MTX_cor_df_padj[with(MTX_cor_df_padj, grepl("^B",Source) & grepl("^F",Target)),])
nrow(bf_occur_MTX) #0


## Network properties
meta_degree <- sort(igraph::degree(all_MTX_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "37"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.125"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_MTX_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  6.12956849677315"
print(paste("mean clustering coefficient = ", transitivity(all_MTX_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.289254170755643"
print(paste("mean betweenness centrality = ", mean(betweenness(all_MTX_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality = 0.00579664186975292"
print(paste("mean closeness centrality = ", mean(closeness(all_MTX_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.020779265984292"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_MTX_net, V(all_MTX_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors = 4.19225967540574"
##

net <- all_MTX_net
MTX_all_deg <- igraph::degree(net,mode="all")
MTX_all_betweenness <- betweenness(net, normalized = TRUE)
MTX_all_closeness <- closeness(net, normalized = TRUE)
MTX_all_transitivity <- transitivity(net, "local", vids = V(net))
names(MTX_all_transitivity)<- V(net)$name
MTX_all_transitivity[is.na(MTX_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
MTX_all_deg.1percent <- MTX_all_deg[MTX_all_deg >= quantile(MTX_all_deg,prob=1-n/100)]
length(MTX_all_deg.1percent) #8

MTX_all_betweenness.1percent <- MTX_all_betweenness[MTX_all_betweenness >= quantile(MTX_all_betweenness,prob=1-n/100)]
length(MTX_all_betweenness.1percent) #8

MTX_all_closeness.1percent <- MTX_all_closeness[MTX_all_closeness >= quantile(MTX_all_closeness,prob=1-n/100)]
length(MTX_all_closeness.1percent) #8

intersect(names(MTX_all_deg.1percent), names(MTX_all_betweenness.1percent)) #"F53_o_NA"             "B3_f_Lachnospiraceae" "F52_Pichia"  
intersect(names(MTX_all_deg.1percent), names(MTX_all_closeness.1percent)) #"F52_Pichia"
intersect(names(MTX_all_betweenness.1percent), names(MTX_all_closeness.1percent))#"F52_Pichia"       "F60_Cladosporium" "F202_Mrakia" 

### network properties of bacteria and fungi
df.MTX.degree<-data.frame(MTX_all_deg)
head(df.MTX.degree)
df.MTX.degree$Group <- "RA"
names(df.MTX.degree)[1] <- c("Degree")

df.MTX.closeness<-data.frame(MTX_all_closeness)
head(df.MTX.closeness)
df.MTX.closeness$Group <- "RA"
names(df.MTX.closeness)[1] <- c("Closeness")

df.MTX.betweenness<-data.frame(MTX_all_betweenness)
head(df.MTX.betweenness)
df.MTX.betweenness$Group <- "RA"
names(df.MTX.betweenness)[1] <- c("Betweenness")

df.MTX.degree$kingdom <- ifelse(grepl("^B",rownames(df.MTX.degree)),'Bacteria', 'Fungi')
df.MTX.closeness$kingdom <- ifelse(grepl("^B",rownames(df.MTX.closeness)),'Bacteria', 'Fungi')
df.MTX.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.MTX.betweenness)),'Bacteria', 'Fungi')



### Plotting hub
df.MTX.degree$OTU_id <- rownames(df.MTX.degree)
df.MTX.betweenness$OTU_id <- rownames(df.MTX.betweenness)
df.MTX.closeness$OTU_id <- rownames(df.MTX.closeness)
df.MTX.net.properties <- merge(df.MTX.degree, df.MTX.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.MTX.net.properties <- merge(df.MTX.net.properties, df.MTX.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))

n<-1
ggplot(df.MTX.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
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
  geom_hline(yintercept=quantile(df.MTX.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.MTX.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

ggplot(df.MTX.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
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
  geom_hline(yintercept=quantile(df.MTX.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.MTX.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

## kingdom comparison
print_closeness(df.MTX.closeness)
print_betweenness(df.MTX.betweenness)
print_degree(df.MTX.degree)


### Control vs RA 
df.control.net.properties$Group2 <- "Control"
df.RA.net.properties$Group2 <- "Not_treated"
df.TNF.net.properties$Group2 <- "TNF_treated"
df.MTX.net.properties$Group2 <- "MTX_treated"

df.net.properties <- rbind(df.control.net.properties,df.RA.net.properties,df.TNF.net.properties,df.MTX.net.properties)
print_degree_treatment <- function(deg){
  
  deg$Group2 <- factor(deg$Group2, levels=c('Control', 'Not_treated','MTX_treated','TNF_treated'))  
  max.degree <- aggregate(deg$Degree, by = list(deg$Group2), max)
  
  colnames(max.degree) <- c("Group2", "maxdegree")
  
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(Degree ~ Group2, data = deg)
  kw$p.value
  
  #library(FSA)
  DT = dunnTest(Degree ~ Group2, data = deg,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  names(dunn)[1] <- "Group2"
  hsd1 <- merge(dunn,max.degree, by = 'Group2')
  
  p<-ggplot(data=deg, aes(x=Group2, y=Degree, fill=Group2)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', kw$p.value))+geom_text(data=hsd1,aes(x=Group2,y=maxdegree,label=Letter),vjust=-1)+
    ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    #scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900")) + 
    theme(aspect.ratio = 1.2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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
print_betweenness_treatment <- function(deg){
  
  deg$Group2 <- factor(deg$Group2, levels=c('Control', 'Not_treated','MTX_treated','TNF_treated'))
  max.degree <- aggregate(deg$Betweenness, by = list(deg$Group2), max)
  
  colnames(max.degree) <- c("Group2", "maxdegree")
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(Betweenness ~ Group2, data = deg)
  kw$p.value
  
  # #library(FSA)
  DT = dunnTest(Betweenness ~ Group2, data = deg,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  # dunn<-cldList(P.adj ~ Comparison,
  #               data = PT,
  #               threshold = 0.05)
  # names(dunn)[1] <- "Group2"
  # hsd1 <- merge(dunn,max.degree, by = 'Group2')
  # 
  p<-ggplot(data=deg, aes(x=Group2, y=Betweenness, fill=Group2)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', kw$p.value))+
    ylab("Betweenness \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    #scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900")) + 
    theme(aspect.ratio = 1.2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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
print_closeness_treatment <- function(deg){
  
  deg$Group2 <- factor(deg$Group2, levels=c('Control', 'Not_treated','MTX_treated','TNF_treated'))
  max.degree <- aggregate(deg$Closeness, by = list(deg$Group2), max)
  
  colnames(max.degree) <- c("Group2", "maxdegree")
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(Closeness ~ Group2, data = deg)
  kw$p.value
  
  #library(FSA)
  DT = dunnTest(Closeness ~ Group2, data = deg,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  names(dunn)[1] <- "Group2"
  hsd1 <- merge(dunn,max.degree, by = 'Group2')
  
  p<-ggplot(data=deg, aes(x=Group2, y=Closeness, fill=Group2)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + geom_text(data=hsd1,aes(x=Group2,y=maxdegree,label=Letter),vjust=-1)+
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ',kw$p.value))+
    ylab("Closeness \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    #scale_fill_manual(values=c("control" = "#6699CC", "RA" = "#CC9900")) + 
    theme(aspect.ratio = 1.2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
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

print_degree_treatment(df.net.properties)
print_betweenness_treatment(df.net.properties)
print_closeness_treatment(df.net.properties)



## Comparison kingdom and diagnosis
df.net.properties
df.net.properties.bac <- subset(df.net.properties, kingdom == "Bacteria")
df.net.properties.fun <- subset(df.net.properties, kingdom == "Fungi")


print_degree_treatment(df.net.properties.bac)
print_betweenness_treatment(df.net.properties.bac)
print_closeness_treatment(df.net.properties.bac)


print_degree_treatment(df.net.properties.fun)
print_betweenness_treatment(df.net.properties.fun)
print_closeness_treatment(df.net.properties.fun)



#### Candida
df.net.properties

RA_30_cor_df_padj$Target[RA_30_cor_df_padj$Source == "F1_Candida"]
#B4_Bifidobacterium               B30_.Ruminococcus..torques.group B74_Dorea
RA_30_cor_df_padj$Cor[RA_30_cor_df_padj$Source == "F1_Candida"]
# -0.4202 -0.4236 -0.4122

RA_30_cor_df_padj$Source[RA_30_cor_df_padj$Target == "B4_Bifidobacterium"]
#B4_Bifidobacterium               B30_.Ruminococcus..torques.group B74_Dorea

RA_30_cor_df_padj$Source[RA_30_cor_df_padj$Target == "B30_.Ruminococcus..torques.group"]
#B52_Collinsella                  B30_[Ruminococcus] torques group B2_Bifidobacterium
RA_30_cor_df_padj$Source[RA_30_cor_df_padj$Target == "B2_Bifidobacterium"]
#F10_Rhodotorula
RA_30_cor_df_padj$Cor[RA_30_cor_df_padj$Target == "B2_Bifidobacterium"]

RA_30_cor_df_padj$Target[RA_30_cor_df_padj$Source == "B2_Bifidobacterium"]
RA_30_cor_df_padj$Cor[RA_30_cor_df_padj$Source == "B2_Bifidobacterium"]
RA_30_cor_df_padj$Cor[RA_30_cor_df_padj$Source == "F1_Candida"]
# -0.4202 -0.4236 -0.4122

control_cor_df_padj$Target[control_cor_df_padj$Source == "F1_Candida"]
#B13_Blautia B39_Blautia B36_Blautia B38_Blautia
control_cor_df_padj$Cor[control_cor_df_padj$Source == "F1_Candida"]
#0.4403 -0.4102  0.4427  0.4662








### Proportion of positive and negative associations
