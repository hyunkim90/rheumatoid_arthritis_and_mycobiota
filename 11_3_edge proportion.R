
## import the edge list


draw_edge_proportion <- function(filename){
  edge.mer <- read.table(file = filename, sep = '\t', header = TRUE)
  head(edge.mer)
  
  edge.mer$BF <- ifelse(grepl("^B",edge.mer$Source) & grepl("^F",edge.mer$Target) | grepl("^F",edge.mer$Source) & grepl("^B",edge.mer$Target),'BF',ifelse(grepl("^B",edge.mer$Source) & grepl("^B",edge.mer$Target),'B','F'))
  edge.mer
  
  edge.mer$PN <- ifelse(edge.mer$positive,'positive','negative')
  edge.mer$count <- 1
  edge.mer %>% arrange(desc(Cor))
  edge.mer %>% arrange(desc(BF))
  
  df.edge.mer <- edge.mer %>% group_by(BF,PN) %>% summarise(count=sum(count))
  
  df.edge.mer$Proportion <- df.edge.mer$count / sum(df.edge.mer$count)
  df.edge.mer
  print(df.edge.mer %>% filter(BF == 'B') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(BF == 'BF') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(BF == 'F') %>% mutate(ratio = Proportion / sum(Proportion)))
  
  p.edge.mer <- ggplot(df.edge.mer, aes(x=BF, y = Proportion, fill = PN)) + 
    geom_bar(stat="identity", width = 0.8, position = 'stack', colour="black") +
    #scale_fill_discrete() +
    scale_fill_manual(values = c('red2','blue2')) +
    
    xlab('')+
    ylab("Proportion of edges \n") +
    #ggtitle("Phylum Community Composition by Sample \n") +
    ## adjust positions
    guides(fill = guide_legend(ncol = 2,reverse = T))+theme(aspect.ratio = 1.5)+
    theme(legend.position="bottom",legend.title=element_blank()) +
    theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(size=15, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    scale_y_continuous(breaks=seq(0,1,0.1))+
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(),plot.background=element_blank())
  
  return(p.edge.mer)
}

## 30 samples
draw_edge_proportion("2_0.4_edge_control.tsv")
draw_edge_proportion("2_0.4_edge_RA_30.tsv")
draw_edge_proportion("2_0.4_edge_control_its1.tsv")
draw_edge_proportion("2_0.4_edge_RA_its1.tsv")
