## Correlation between RA-related variables and functional properties predicted using Tax4Fun2
predictedpathway <- read.table('pathway_prediction_wo_description.tsv', sep = '\t', header = T)
head(predictedpathway)
rownames(predictedpathway) <- predictedpathway[,1]

#remove pathway names
predictedpathway.2<- predictedpathway[,-1]
head(predictedpathway.2)
predictedpathway.3 <- t(predictedpathway.2)
head(predictedpathway.3)

b.meta
predictedpathway.4<-data.frame(predictedpathway.3)
head(predictedpathway.4)


b.meta.w.pathway <- merge(b.meta,predictedpathway.4, by = 'row.names')


##List of lipid metabolism
lipid<-pathwaylist$pathway[which(pathwaylist$level2 == "Lipid metabolism")]
pathwaylist$level1[which(pathwaylist$level2 == "Lipid metabolism")]

##
predictedpathway.4.lipid.t <- subset(t(predictedpathway.4), row.names(t(predictedpathway.4))%in%lipid)
head(predictedpathway.4.lipid.t)

predictedpathway.4.lipid <-t(predictedpathway.4.lipid.t)
predictedpathway.4.lipid <- data.frame(predictedpathway.4.lipid)

sapply(predictedpathway.4.lipid, class) ## soil properties are class 'factor'
indx.div <- sapply(predictedpathway.4.lipid, is.factor)
predictedpathway.4.lipid[indx.div] <- lapply(predictedpathway.4.lipid[indx.div], function(x) as.numeric(as.character(x)))
head(predictedpathway.4.lipid)
sapply(predictedpathway.4.lipid, class)


b.meta.w.pathway <- merge(b.meta.ra.edit,predictedpathway.4.lipid, by = 'row.names')

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00061, type="spearman")
cor_5$r
cor_5$P

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00062, type="spearman")
cor_5$r
cor_5$P

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00071, type="spearman")
cor_5$r
cor_5$P

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00072, type="spearman")
cor_5$r
cor_5$P

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00100, type="spearman")
cor_5$r
cor_5$P

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00120, type="spearman")
cor_5$r
cor_5$P

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00121, type="spearman")
cor_5$r
cor_5$P

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00140, type="spearman")
cor_5$r #-0.247529 #Steroid hormone biosynthesis
cor_5$P #0.02495664

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00561, type="spearman")
cor_5$r
cor_5$P

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00564, type="spearman")
cor_5$r
cor_5$P

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00565, type="spearman")
cor_5$r #-0.2502109 #Ether lipid metabolism
cor_5$P #0.02338017

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00590, type="spearman")
cor_5$r 
cor_5$P 

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00591, type="spearman")
cor_5$r 
cor_5$P 

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00592, type="spearman")
cor_5$r 
cor_5$P 

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko00600, type="spearman")
cor_5$r 
cor_5$P 

cor_5 <- Hmisc::rcorr(b.meta.w.pathway$Total_Cholesterol, b.meta.w.pathway$ko01040, type="spearman")
cor_5$r #-0.2193351 #Biosynthesis of unsaturated fatty acids
cor_5$P #0.04772018


ggplot(b.meta.w.pathway, aes(x=Total_Cholesterol, y=log(ko01040))) +
  xlab('\n Total cholesterol')+
  ylab("ko01040 \n") +
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
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred")

lmHeight = lm(Total_Cholesterol~ko01040, data = b.meta.w.pathway) #Create the linear regression
summary(lmHeight) #Adjusted R-squared:  0.05474 p-value: 0.01942



ggplot(b.meta.w.pathway, aes(x=Total_Cholesterol, y=log(ko00565))) +
  xlab('\n Total cholesterol')+
  ylab("ko00565 \n") +
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
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred")

lmHeight = lm(Total_Cholesterol~ko00565, data = b.meta.w.pathway) #Create the linear regression
summary(lmHeight) #Adjusted R-squared:  0.04045 p-value: 0.03878


ggplot(b.meta.w.pathway, aes(x=Total_Cholesterol, y=log(ko00140))) +
  xlab('\n Total cholesterol')+
  ylab("ko00140 \n") +
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
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred")

lmHeight = lm(Total_Cholesterol~ko00140, data = b.meta.w.pathway) #Create the linear regression
summary(lmHeight) #Adjusted R-squared:  0.0424  p-value: 0.03526
