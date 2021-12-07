##Calculate contribution of each factor to compositional varations estimated by beta diversity (PCoA)
order.diagnosis <- c("control", "RA")
order.active_disease <- c('control', '0', "1", '2','3')



b.meta <- sample_data(bac.clean.log)
b.meta <- data.frame(b.meta)

b.meta$exercise_light <- as.character(b.meta$exercise_light)
b.meta$exercise_oderate <- as.character(b.meta$exercise_oderate)
b.meta$exercise_heavy <- as.character(b.meta$exercise_heavy)
b.meta$active_disease <- as.character(b.meta$active_disease)
b.meta$Duration <- as.character(b.meta$Duration)
b.meta$TNF_inhibitor <- as.character(b.meta$TNF_inhibitor)
b.meta$GC <- as.character(b.meta$GC)
b.meta$MTX <- as.character(b.meta$MTX)
b.meta$seropositive_RA <- as.character(b.meta$seropositive_RA)


b.meta$Diagnosis <- factor(b.meta$Diagnosis, levels = order.diagnosis)
b.meta$active_disease <- factor(b.meta$active_disease, levels = order.active_disease)

sample_data(bac.clean.log) <- sample_data(b.meta)
sample_data(bac.clean.ss) <- sample_data(b.meta)
sample_data(bac.clean.nolog) <-sample_data(b.meta)

bac.clean.log
bac.clean.nolog
bac.clean.ss.rel

bray1.bac <-  ordinate(bac.clean.log, 'PCoA', 'bray')
sample_data(bac.clean.log) <- sample_data(b.meta)
#cap_uni <-  ordinate(bac.clean.log, 'CAP', 'unifrac', weight = TRUE, ~ Diagnosis)

write.csv(bray1.bac$vectors, "Supplementary Fig. A_PCoA.csv")
write.csv(b.meta, "Source data for PCoA and CAP_bac.csv")
write.csv(f.meta, "Source data for PCoA and CAP_fun.csv")
#bray3.bac <-  ordinate(bac.clean.log, 'PCoA', 'unifrac')

##Unconstrained PCoA
#Control vs RA
plot_ordination(bac.clean.log, bray1.bac, type = "samples", color='Diagnosis', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
#Treatment
plot_ordination(bac.clean.log, bray1.bac, type = "samples", color='treatment', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#Age
plot_ordination(bac.clean.ss.rel, bray1.bac, type = "samples", color='Age', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


## MTX
plot_ordination(bac.clean.log, bray1.bac, type = "samples", color='MTX', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### CAP (Canonical analysis of principals)
##Relative abundance
b.cap.diagnosis <- ordinate(bac.clean.ss.rel, "CAP", "bray", ~ Diagnosis)

perm_anova.ord <- anova.cca(b.cap.diagnosis)
perm_anova.ord2 <- permutest(b.cap.diagnosis)

## Plotting
# Diagnosis
plot.b.cap.type     <- plot_ordination(bac.clean.ss.rel, b.cap.diagnosis, type = "samples", color='Diagnosis', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.type


##Duration
b.cap.duration <- ordinate(bac.clean.ss.rel, "CAP", "bray", ~ Duration)

perm_anova.ord <- anova.cca(b.cap.duration)
perm_anova.ord2 <- permutest(b.cap.duration)

## Plotting
plot.b.cap.duration     <- plot_ordination(bac.clean.ss.rel, b.cap.duration, type = "samples", color='Duration', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.duration


##Active disease
b.cap.active <- ordinate(bac.clean.ss.rel, "CAP", "bray", ~ active_disease)

perm_anova.ord <- anova.cca(b.cap.active)
perm_anova.ord2 <- permutest(b.cap.active)

## Plotting
plot.b.cap.active     <- plot_ordination(bac.clean.ss.rel, b.cap.active, type = "samples", color='active_disease', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.active


##Age
b.cap.age <- ordinate(bac.clean.ss.rel, "CAP", "bray", ~ Age)

perm_anova.ord <- anova.cca(b.cap.age)
perm_anova.ord2 <- permutest(b.cap.age)

## Plotting
plot.b.cap.age     <- plot_ordination(bac.clean.ss.rel, b.cap.age, type = "samples", color='Age', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+ 
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.age

##Normalized abundance
b.cap.diagnosis <- ordinate(bac.clean.nolog, "CAP", "bray", ~ Diagnosis)

perm_anova.ord <- anova.cca(b.cap.diagnosis)
perm_anova.ord2 <- permutest(b.cap.diagnosis)

## Plotting
# Diagnosis
plot.b.cap.type     <- plot_ordination(bac.clean.nolog, b.cap.diagnosis, type = "samples", color='Diagnosis', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.type


##Duration
b.cap.duration <- ordinate(bac.clean.nolog, "CAP", "bray", ~ Duration)

perm_anova.ord <- anova.cca(b.cap.duration)
perm_anova.ord2 <- permutest(b.cap.duration)

## Plotting
plot.b.cap.duration     <- plot_ordination(bac.clean.nolog, b.cap.duration, type = "samples", color='Duration', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.duration


##Active disease
b.cap.active <- ordinate(bac.clean.nolog, "CAP", "bray", ~ active_disease)

perm_anova.ord <- anova.cca(b.cap.active)
perm_anova.ord2 <- permutest(b.cap.active)

## Plotting
plot.b.cap.active     <- plot_ordination(bac.clean.nolog, b.cap.active, type = "samples", color='active_disease', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.active


##Age
b.cap.age <- ordinate(bac.clean.nolog, "CAP", "bray", ~ Age)

perm_anova.ord <- anova.cca(b.cap.age)
perm_anova.ord2 <- permutest(b.cap.age)

## Plotting
plot.b.cap.age     <- plot_ordination(bac.clean.nolog, b.cap.age, type = "samples", color='Age', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.age

##Normalized abundance and log transformed
b.cap.diagnosis <- ordinate(bac.clean.log.2, "CAP", "bray", ~ Diagnosis)

b.cap.diagnosis$CCA$Xbar

write.csv(f.cap.diagnosis$CCA$Xbar,"Fig. 2A_fun.csv")

perm_anova.ord <- anova.cca(b.cap.diagnosis)
perm_anova.ord2 <- permutest(b.cap.diagnosis)

## Plotting
# Diagnosis
plot.b.cap.type     <- plot_ordination(bac.clean.log.2, b.cap.diagnosis, type = "samples", color='Blautia', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.type


##Duration
b.cap.duration <- ordinate(bac.clean.log, "CAP", "bray", ~ Duration)

perm_anova.ord <- anova.cca(b.cap.duration)
perm_anova.ord2 <- permutest(b.cap.duration)

## Plotting
plot.b.cap.duration     <- plot_ordination(bac.clean.log, b.cap.duration, type = "samples", color='Duration', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.duration


##Active disease
b.cap.active <- ordinate(bac.clean.log, "CAP", "bray", ~ active_disease)

perm_anova.ord <- anova.cca(b.cap.active)
perm_anova.ord2 <- permutest(b.cap.active)

## Plotting
plot.b.cap.active     <- plot_ordination(bac.clean.log, b.cap.active, type = "samples", color='active_disease', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.active


##Age
b.cap.age <- ordinate(bac.clean.log, "CAP", "bray", ~ Age)

perm_anova.ord <- anova.cca(b.cap.age)
perm_anova.ord2 <- permutest(b.cap.age)

## Plotting
plot.b.cap.age     <- plot_ordination(bac.clean.log, b.cap.age, type = "samples", color='Age', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.age


##Active disease
b.cap.active <- ordinate(bac.clean.log, "CAP", "bray", ~ TX)

perm_anova.ord <- anova.cca(b.cap.active)
perm_anova.ord2 <- permutest(b.cap.active)

## Plotting
plot.b.cap.active     <- plot_ordination(bac.clean.log, b.cap.active, type = "samples", color='TX', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.active


##MTX
b.cap.active <- ordinate(bac.clean.log, "CAP", "bray", ~ MTX)

perm_anova.ord <- anova.cca(b.cap.active)
perm_anova.ord2 <- permutest(b.cap.active)

## Plotting
plot.b.cap.active     <- plot_ordination(bac.clean.log, b.cap.active, type = "samples", color='MTX', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.active

## Only RA samples
bac.clean.log.RA

b.meta.ra$active_disease <- as.character(b.meta.ra$active_disease)
sample_data(bac.clean.log.RA) <- sample_data(b.meta.ra) 

#Active disease
b.cap.active <- ordinate(bac.clean.log.RA, "CAP", "bray", ~ active_disease)

perm_anova.ord <- anova.cca(b.cap.active)
perm_anova.ord2 <- permutest(b.cap.active)

## Plotting
# Diagnosis
plot.b.cap.active     <- plot_ordination(bac.clean.log.RA, b.cap.active, type = "samples", color='active_disease', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.active


#Duration
b.cap.duration <- ordinate(bac.clean.log.RA, "CAP", "bray", ~ Duration)

perm_anova.ord <- anova.cca(b.cap.duration)
perm_anova.ord2 <- permutest(b.cap.duration)

## Plotting
# Diagnosis
plot.b.cap.duration     <- plot_ordination(bac.clean.log.RA, b.cap.duration, type = "samples", color='Duration', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.duration


###PERMANOVA
## ALL
bac.clean.ss.rel
bac.clean.log
bac.clean.nolog

b.otu <- otu_table(bac.clean.log)
b.meta <- sample_data(bac.clean.log)
b.meta <- data.frame(b.meta)

##new variable for active disease
b.meta$active_disease_2 <- b.meta$active_disease
b.meta$active_disease_2<-as.character(b.meta$active_disease_2)
b.meta$active_disease_2[which(b.meta$active_disease_2 %in% c("0", "1"))] <- "low"
b.meta$active_disease_2[which(b.meta$active_disease_2 %in% c("2", "3"))] <- "high"

b.permanova <- adonis(formula = t(b.otu) ~ (Diagnosis+RA_factor+anti_CCP+seropositive_RA+GC+MTX+TNF_inhibitor+exercise_light+exercise_oderate+exercise_heavy+active_disease_2+Duration+Age+BMI+autoimmune_disease+diabetes), data = b.meta, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (TNF_inhibitor+TX), data = b.meta, permutations=9999, method = "bray")
b.permanova

##Only RA
bac.clean.ss.RA <- subset_samples(bac.clean.ss, Diagnosis == "RA") 
bac.clean.ss.RA <- phyloseq::filter_taxa(bac.clean.ss.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs

bac.clean.ss.RA.rel <- microbiome::transform(bac.clean.ss.RA, "compositional")
t(otu_table(bac.clean.ss.RA.rel))

bac.clean.log.RA <- subset_samples(bac.clean.log, Diagnosis == "RA") 
bac.clean.log.RA <- phyloseq::filter_taxa(bac.clean.log.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs

bac.clean.nolog.RA <- subset_samples(bac.clean.nolog, Diagnosis == "RA") 
bac.clean.nolog.RA <- phyloseq::filter_taxa(bac.clean.nolog.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs

b.meta.ra <- sample_data(bac.clean.ss.RA)
b.meta.ra <- data.frame(b.meta.ra)

b.meta.ra$active_disease_2 <- b.meta.ra$active_disease
b.meta.ra$active_disease_2<-as.character(b.meta.ra$active_disease_2)
b.meta.ra$active_disease_2[which(b.meta.ra$active_disease_2 %in% c("0", "1"))] <- "low"
b.meta.ra$active_disease_2[which(b.meta.ra$active_disease_2 %in% c("2", "3"))] <- "high"

head(b.meta.ra)
class(b.meta.ra$Ds_Duration)
b.meta.ra$Ds_Duration <- as.numeric(as.character(b.meta.ra$Ds_Duration))

b.otu.ra <- otu_table(bac.clean.log.RA)
b.permanova.ra <- adonis(formula = t(b.otu.ra) ~ (RA_factor+anti_CCP+seropositive_RA+GC+TX+TNF_inhibitor+functionalfood_type+exercise_light+exercise_oderate+exercise_heavy+active_disease_2+Duration+Age+BMI+autoimmune_disease+diabetes), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova.ra

b.permanova.ra <- adonis(formula = t(b.otu.ra) ~ (RA_factor+anti_CCP+BMI+Age+Ds_Duration), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova.ra

b.permanova.ra <- adonis(formula = t(b.otu.ra) ~ (TNF_inhibitor+TX), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova.ra

b.meta$Gender
###Only women samples
bac.clean.log.W <- subset_samples(bac.clean.log, Gender == "2")
bac.clean.log.W <- phyloseq::filter_taxa(bac.clean.log.W, function(x) sum(x) != 0, TRUE)


bac.clean.nolog.W <- subset_samples(bac.clean.log, Gender == "2")
bac.clean.nolog.W <- phyloseq::filter_taxa(bac.clean.nolog.W, function(x) sum(x) != 0, TRUE)


bac.clean.ss.rel.W <- subset_samples(bac.clean.ss.rel, Gender == "2")
bac.clean.ss.rel.W <- phyloseq::filter_taxa(bac.clean.ss.rel.W, function(x) sum(x) != 0, TRUE)


b.meta.women$Ds_Duration <- as.numeric(as.character(b.meta.women$Ds_Duration))
b.meta.women$Ds_Duration[is.na(b.meta.women$Ds_Duration)] <- 0

b.otu <- otu_table(bac.clean.log.W)
b.permanova <- adonis(formula = t(b.otu) ~ (Diagnosis+RA_factor+anti_CCP+seropositive_RA+GC+TX+TNF_inhibitor+exercise_light+exercise_oderate+exercise_heavy+active_disease_2+Duration+Age+BMI+autoimmune_disease+diabetes), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (RA_factor+anti_CCP+BMI+Age+Ds_Duration+active_disease+active_disease_2), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (TNF_inhibitor+TX), data = b.meta.women, permutations=9999, method = "bray")
b.permanova


### RA and women
bac.clean.log.RA.W <- subset_samples(bac.clean.log.RA, Gender == "2")
bac.clean.log.RA.W <- phyloseq::filter_taxa(bac.clean.log.RA.W, function(x) sum(x) != 0, TRUE)


bac.clean.nolog.RA.W <- subset_samples(bac.clean.nolog.RA, Gender == "2")
bac.clean.nolog.RA.W <- phyloseq::filter_taxa(bac.clean.nolog.RA.W, function(x) sum(x) != 0, TRUE)


bac.clean.ss.RA.rel.W <- subset_samples(bac.clean.ss.RA.rel, Gender == "2")
bac.clean.ss.RA.rel.W <- phyloseq::filter_taxa(bac.clean.ss.RA.rel.W, function(x) sum(x) != 0, TRUE)

b.meta.ra.women$Ds_Duration <- as.numeric(as.character(b.meta.ra.women$Ds_Duration))
b.meta.ra.women$Ds_Duration[is.na(b.meta.ra.women$Ds_Duration)] <- 0

b.otu.ra <- otu_table(bac.clean.log.RA.W)
b.permanova <- adonis(formula = t(b.otu.ra) ~ (RA_factor+Ds_Duration+anti_CCP+seropositive_RA+GC+TX+TNF_inhibitor+exercise_light+exercise_oderate+exercise_heavy+active_disease_2+Duration+Age+BMI+autoimmune_disease+diabetes), data = b.meta.ra.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu.ra) ~ (RA_factor+anti_CCP+BMI+Age+Ds_Duration+active_disease+active_disease_2), data = b.meta.ra.women, permutations=9999, method = "bray")
b.permanova


b.permanova <- adonis(formula = t(b.otu.ra) ~ (TNF_inhibitor+TX), data = b.meta.ra.women, permutations=9999, method = "bray")
b.permanova


dev.off()




############### Fungal community #################


f.meta <- sample_data(fun.clean.log)
f.meta <- data.frame(f.meta)

f.meta$exercise_light <- as.character(f.meta$exercise_light)
f.meta$exercise_oderate <- as.character(f.meta$exercise_oderate)
f.meta$exercise_heavy <- as.character(f.meta$exercise_heavy)
f.meta$active_disease <- as.character(f.meta$active_disease)
f.meta$Duration <- as.character(f.meta$Duration)
f.meta$TNF_inhibitor <- as.character(f.meta$TNF_inhibitor)
f.meta$GC <- as.character(f.meta$GC)
f.meta$MTX <- as.character(f.meta$MTX)
f.meta$seropositive_RA <- as.character(f.meta$seropositive_RA)
f.meta$periodental_score <- as.character(f.meta$periodental_score)
f.meta$periodental <- as.character(f.meta$periodental)

### factor to numeric
f.meta$Age <- as.numeric(as.character(f.meta$Age))
f.meta$BMI <- as.numeric(as.character(f.meta$BMI))
f.meta$RA_factor <- as.numeric(as.character(f.meta$RA_factor))
f.meta$anti_CCP <- as.numeric(as.character(f.meta$anti_CCP))
f.meta$MTX_dose <- as.numeric(as.character(f.meta$MTX_dose))
f.meta$CRP <- as.numeric(as.character(f.meta$CRP))
f.meta$ESR <- as.numeric(as.character(f.meta$ESR))
f.meta$BUN <- as.numeric(as.character(f.meta$BUN))
f.meta$Cr <- as.numeric(as.character(f.meta$Cr))

f.meta$Total_cholesterol <- as.numeric(as.character(f.meta$Total_cholesterol))
f.meta$Triglyceride <- as.numeric(as.character(f.meta$Triglyceride))
f.meta$HDL <- as.numeric(as.character(f.meta$HDL))



##Ordering

f.meta$Diagnosis <- factor(f.meta$Diagnosis, levels = order.diagnosis)
f.meta$active_disease <- factor(f.meta$active_disease, levels = order.active_disease)

sample_data(fun.clean.log) <- sample_data(f.meta)
sample_data(fun.clean.ss) <- sample_data(f.meta)
sample_data(fun.clean.nolog) <-sample_data(f.meta)

fun.clean.log
fun.clean.nolog
fun.clean.ss.rel

bray1.fun <-  ordinate(fun.clean.log.2, 'PCoA', 'bray')
write.csv(bray1.fun$vectors, "Supplementary Fig. A_PCoA_fun.csv")
#cap_uni <-  ordinate(fun.clean.log, 'CAP', 'unifrac', weight = TRUE, ~ Diagnosis)


#bray3.fun <-  ordinate(fun.clean.log, 'PCoA', 'unifrac')

##Unconstrained PCoA
#Control vs RA
plot_ordination(fun.clean.log, bray1.fun, type = "samples", color='Diagnosis', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#Age
plot_ordination(fun.clean.log.2.its1, bray1.fun, type = "samples", color='Age', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


## MTX
plot_ordination(fun.clean.log.2.its1, bray1.fun, type = "samples", color='MTX', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

## TNF inhibitor
plot_ordination(fun.clean.log, bray1.fun, type = "samples", color='TNF_inhibitor', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


##Active disease
plot_ordination(fun.clean.log, bray1.fun, type = "samples", color='active_disease', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### CAP (Canonical analysis of principals)


##Relative abundance
fun.clean.ss.rel <- microbiome::transform(fun.clean.ss, "compositional")

f.cap.diagnosis <- ordinate(fun.clean.ss.rel, "CAP", "bray", ~ Diagnosis)

perm_anova.ord <- anova.cca(f.cap.diagnosis)
perm_anova.ord2 <- permutest(f.cap.diagnosis)

## Plotting
# Diagnosis
plot.f.cap.type     <- plot_ordination(fun.clean.ss.rel, f.cap.diagnosis, type = "samples", color='Diagnosis', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type


##Duration
f.cap.duration <- ordinate(fun.clean.ss.rel, "CAP", "bray", ~ Duration)

perm_anova.ord <- anova.cca(f.cap.duration)
perm_anova.ord2 <- permutest(f.cap.duration)

## Plotting
plot.f.cap.duration     <- plot_ordination(fun.clean.ss.rel, f.cap.duration, type = "samples", color='Duration', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.duration


##Active disease
f.cap.active <- ordinate(fun.clean.ss.rel, "CAP", "bray", ~ active_disease)

perm_anova.ord <- anova.cca(f.cap.active)
perm_anova.ord2 <- permutest(f.cap.active)

## Plotting
plot.f.cap.active     <- plot_ordination(fun.clean.ss.rel, f.cap.active, type = "samples", color='active_disease', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.active


##Age
f.cap.age <- ordinate(fun.clean.ss.rel, "CAP", "bray", ~ Age)

perm_anova.ord <- anova.cca(f.cap.age)
perm_anova.ord2 <- permutest(f.cap.age)

## Plotting
plot.f.cap.age     <- plot_ordination(fun.clean.ss.rel, f.cap.age, type = "samples", color='Age', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+ 
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.age

##Normalized abundance
f.cap.diagnosis <- ordinate(fun.clean.nolog, "CAP", "bray", ~ Diagnosis)

perm_anova.ord <- anova.cca(f.cap.diagnosis)
perm_anova.ord2 <- permutest(f.cap.diagnosis)

## Plotting
# Diagnosis
plot.f.cap.type     <- plot_ordination(fun.clean.nolog, f.cap.diagnosis, type = "samples", color='Diagnosis', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type


##Duration
f.cap.duration <- ordinate(fun.clean.nolog, "CAP", "bray", ~ Duration)

perm_anova.ord <- anova.cca(f.cap.duration)
perm_anova.ord2 <- permutest(f.cap.duration)

## Plotting
plot.f.cap.duration     <- plot_ordination(fun.clean.nolog, f.cap.duration, type = "samples", color='Duration', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.duration


##Active disease
f.cap.active <- ordinate(fun.clean.nolog, "CAP", "bray", ~ active_disease)

perm_anova.ord <- anova.cca(f.cap.active)
perm_anova.ord2 <- permutest(f.cap.active)

## Plotting
plot.f.cap.active     <- plot_ordination(fun.clean.nolog, f.cap.active, type = "samples", color='active_disease', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.active


##Age
f.cap.age <- ordinate(fun.clean.nolog, "CAP", "bray", ~ Age)

perm_anova.ord <- anova.cca(f.cap.age)
perm_anova.ord2 <- permutest(f.cap.age)

## Plotting
plot.f.cap.age     <- plot_ordination(fun.clean.nolog, f.cap.age, type = "samples", color='Age', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.age

##Normalized abundance and log transformed
f.cap.diagnosis <- ordinate(fun.clean.log.2.its1, "CAP", "bray", ~ Diagnosis)

perm_anova.ord <- anova.cca(f.cap.diagnosis)
perm_anova.ord2 <- permutest(f.cap.diagnosis)

## Plotting
# Diagnosis
plot.f.cap.type     <- plot_ordination(fun.clean.log.2.its1, f.cap.diagnosis, type = "samples", color='Diagnosis', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type

plot.f.cap.type     <- plot_ordination(fun.clean.log.2.its1, f.cap.diagnosis, type = "samples", color='Candida', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+scale_colour_gradient(low = "#CCCCFF",
                                             high = "#333366",
                                             space = "Lab",
                                             na.value = "grey50",
                                             guide = "colourbar",
                                             aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type

plot.f.cap.type     <- plot_ordination(fun.clean.log.2.its1, f.cap.diagnosis, type = "samples", color='Aspergillus', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+scale_colour_gradient(low = "#CCFFFF",
                                             high = "#003333",
                                             space = "Lab",
                                             na.value = "grey50",
                                             guide = "colourbar",
                                             aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type


##Duration
f.cap.duration <- ordinate(fun.clean.log, "CAP", "bray", ~ Duration)

perm_anova.ord <- anova.cca(f.cap.duration)
perm_anova.ord2 <- permutest(f.cap.duration)

## Plotting
plot.f.cap.duration     <- plot_ordination(fun.clean.log, f.cap.duration, type = "samples", color='Duration', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.duration


##Active disease
f.cap.active <- ordinate(fun.clean.log, "CAP", "bray", ~ active_disease)

perm_anova.ord <- anova.cca(f.cap.active)
perm_anova.ord2 <- permutest(f.cap.active)

## Plotting
plot.f.cap.active     <- plot_ordination(fun.clean.log, f.cap.active, type = "samples", color='active_disease', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.active


##Age
f.cap.age <- ordinate(fun.clean.log, "CAP", "bray", ~ Age)

perm_anova.ord <- anova.cca(f.cap.age)
perm_anova.ord2 <- permutest(f.cap.age)

## Plotting
plot.f.cap.age     <- plot_ordination(fun.clean.log, f.cap.age, type = "samples", color='Age', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.age


##MTX
f.cap.active <- ordinate(fun.clean.log, "CAP", "bray", ~ MTX)

perm_anova.ord <- anova.cca(f.cap.active)
perm_anova.ord2 <- permutest(f.cap.active)

## Plotting
plot.f.cap.active     <- plot_ordination(fun.clean.log, f.cap.active, type = "samples", color='MTX', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.active


## Only RA samples
fun.clean.log.RA

f.meta.ra$active_disease <- as.character(f.meta.ra$active_disease)
sample_data(fun.clean.log.RA) <- sample_data(f.meta.ra) 

#Active disease
f.cap.active <- ordinate(fun.clean.log.RA, "CAP", "bray", ~ active_disease)

perm_anova.ord <- anova.cca(f.cap.active)
perm_anova.ord2 <- permutest(f.cap.active)

## Plotting
# active disease
plot.f.cap.active     <- plot_ordination(fun.clean.log.RA, f.cap.active, type = "samples", color='active_disease', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.active


#Duration
f.cap.duration <- ordinate(fun.clean.log.RA, "CAP", "bray", ~ Duration)

perm_anova.ord <- anova.cca(f.cap.duration)
perm_anova.ord2 <- permutest(f.cap.duration)

## Plotting
# Diagnosis
plot.f.cap.duration     <- plot_ordination(fun.clean.log.RA, f.cap.duration, type = "samples", color='Duration', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.duration



#Duration (quantitative variable)
f.cap.duration <- ordinate(fun.clean.log.RA, "CAP", "bray", ~ Ds_Duration)

perm_anova.ord <- anova.cca(f.cap.duration)
perm_anova.ord2 <- permutest(f.cap.duration)

## Plotting
# Diagnosis
plot.f.cap.duration     <- plot_ordination(fun.clean.log.RA, f.cap.duration, type = "samples", color='Duration', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.duration



########## Subset samples ###########
######### Bacterial community ##########
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


######### fungal community ##########
## RA
fun.clean.ss.RA <- subset_samples(fun.clean.ss, Diagnosis == "RA") 
fun.clean.ss.RA <- phyloseq::filter_taxa(fun.clean.ss.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs

fun.clean.ss.RA.rel <- microbiome::transform(fun.clean.ss.RA, "compositional")

fun.clean.log.RA <- subset_samples(fun.clean.log, Diagnosis == "RA") 
fun.clean.log.RA <- phyloseq::filter_taxa(fun.clean.log.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs

fun.clean.nolog.RA <- subset_samples(fun.clean.nolog, Diagnosis == "RA") 
fun.clean.nolog.RA <- phyloseq::filter_taxa(fun.clean.nolog.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs


##Women
fun.clean.log.W <- subset_samples(fun.clean.log, Gender == "2")
fun.clean.log.W <- phyloseq::filter_taxa(fun.clean.log.W, function(x) sum(x) != 0, TRUE)


fun.clean.nolog.W <- subset_samples(fun.clean.log, Gender == "2")
fun.clean.nolog.W <- phyloseq::filter_taxa(fun.clean.nolog.W, function(x) sum(x) != 0, TRUE)


fun.clean.ss.rel.W <- subset_samples(fun.clean.ss.rel, Gender == "2")
fun.clean.ss.rel.W <- phyloseq::filter_taxa(fun.clean.ss.rel.W, function(x) sum(x) != 0, TRUE)


##RA and Women
fun.clean.log.RA.W <- subset_samples(fun.clean.log.RA, Gender == "2")
fun.clean.log.RA.W <- phyloseq::filter_taxa(fun.clean.log.RA.W, function(x) sum(x) != 0, TRUE)


fun.clean.nolog.RA.W <- subset_samples(fun.clean.nolog.RA, Gender == "2")
fun.clean.nolog.RA.W <- phyloseq::filter_taxa(fun.clean.nolog.RA.W, function(x) sum(x) != 0, TRUE)


fun.clean.ss.RA.rel.W <- subset_samples(fun.clean.ss.RA.rel, Gender == "2")
fun.clean.ss.RA.rel.W <- phyloseq::filter_taxa(fun.clean.ss.RA.rel.W, function(x) sum(x) != 0, TRUE)

f.meta
f.meta$active_disease_2 <- f.meta$active_disease
f.meta$active_disease_2<-as.character(f.meta$active_disease_2)
f.meta$active_disease_2[which(f.meta$active_disease_2 %in% c("0", "1"))] <- "low"
f.meta$active_disease_2[which(f.meta$active_disease_2 %in% c("2", "3"))] <- "high"

f.meta.ra <- subset(f.meta, Diagnosis == "RA")
f.meta.ra.women <- subset(f.meta.ra, Gender == "2")
f.meta.women <- subset(f.meta, Gender == "2")


### ITS1
## RA
fun.its1.clean.ss.RA <- subset_samples(fun.its1.clean.ss, Diagnosis == "RA") 
fun.its1.clean.ss.RA <- phyloseq::filter_taxa(fun.its1.clean.ss.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs

fun.its1.clean.ss.RA.rel <- microbiome::transform(fun.its1.clean.ss.RA, "compositional")

fun.its1.clean.log.RA <- subset_samples(fun.its1.clean.log, Diagnosis == "RA") 
fun.its1.clean.log.RA <- phyloseq::filter_taxa(fun.its1.clean.log.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs

fun.its1.clean.nolog.RA <- subset_samples(fun.its1.clean.nolog, Diagnosis == "RA") 
fun.its1.clean.nolog.RA <- phyloseq::filter_taxa(fun.its1.clean.nolog.RA, function(x) sum(x) != 0, TRUE) # 1339 OTUs


##Women
fun.its1.clean.log.W <- subset_samples(fun.its1.clean.log, Gender == "2")
fun.its1.clean.log.W <- phyloseq::filter_taxa(fun.its1.clean.log.W, function(x) sum(x) != 0, TRUE)


fun.its1.clean.nolog.W <- subset_samples(fun.its1.clean.log, Gender == "2")
fun.its1.clean.nolog.W <- phyloseq::filter_taxa(fun.its1.clean.nolog.W, function(x) sum(x) != 0, TRUE)


fun.its1.clean.ss.rel.W <- subset_samples(fun.its1.clean.ss.rel, Gender == "2")
fun.its1.clean.ss.rel.W <- phyloseq::filter_taxa(fun.its1.clean.ss.rel.W, function(x) sum(x) != 0, TRUE)


##RA and Women
fun.its1.clean.log.RA.W <- subset_samples(fun.its1.clean.log.RA, Gender == "2")
fun.its1.clean.log.RA.W <- phyloseq::filter_taxa(fun.its1.clean.log.RA.W, function(x) sum(x) != 0, TRUE)


fun.its1.clean.nolog.RA.W <- subset_samples(fun.its1.clean.nolog.RA, Gender == "2")
fun.its1.clean.nolog.RA.W <- phyloseq::filter_taxa(fun.its1.clean.nolog.RA.W, function(x) sum(x) != 0, TRUE)


fun.its1.clean.ss.RA.rel.W <- subset_samples(fun.its1.clean.ss.RA.rel, Gender == "2")
fun.its1.clean.ss.RA.rel.W <- phyloseq::filter_taxa(fun.its1.clean.ss.RA.rel.W, function(x) sum(x) != 0, TRUE)



##PERMANOVA ######
### Bacteria #####
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


#Women
b.otu <- otu_table(bac.clean.ss.rel.W)
b.otu <- otu_table(bac.clean.nolog.W)
b.otu <- otu_table(bac.clean.log.W)

b.permanova <- adonis(formula = t(b.otu) ~ (MTX+MTX_dose+TNF_inhibitor+enstruation), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (enstruation), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (periodental_disease+gu_edea_bleeding+sore_on_the_tongue_hollow), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (periodental+periodental_score), data = b.meta.women, permutations=9999, method = "bray")
b.permanova


## RA
b.otu.ra <- otu_table(bac.clean.ss.RA.rel)
b.otu.ra <- otu_table(bac.clean.nolog.RA)
b.otu.ra <- otu_table(bac.clean.log.RA)

b.permanova <- adonis(formula = t(b.otu.ra) ~ (MTX+MTX_dose), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu.ra) ~ (enstruation), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu.ra) ~ (periodental_disease+gu_edea_bleeding+sore_on_the_tongue_hollow), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu.ra) ~ (periodental+periodental_score), data = b.meta.ra, permutations=9999, method = "bray")
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
rm1<-b.meta.ra$SampleID[-which(b.meta.ra$Total_cholesterol == 0)]
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

b.permanova <- adonis(formula = b.otu.ra.t ~ (Total_cholesterol), data = b.meta.ra.edit, permutations=9999, method = "bray")
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

b.permanova <- adonis(formula = b.otu.ra.t ~ (Total_cholesterol+HDL+Triglyceride), data = b.meta.ra.edit, permutations=9999, method = "bray")
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


##new variable for active disease
b.meta$active_disease_2 <- b.meta$active_disease
b.meta$active_disease_2<-as.character(b.meta$active_disease_2)
b.meta$active_disease_2[which(b.meta$active_disease_2 %in% c("0", "1"))] <- "low"
b.meta$active_disease_2[which(b.meta$active_disease_2 %in% c("2", "3"))] <- "high"

b.otu <- otu_table(bac.clean.log)

b.permanova <- adonis(formula = t(b.otu) ~ (Diagnosis+MTX+RA_factor+anti_CCP+seropositive_RA+GC+TNF_inhibitor+exercise_light+exercise_oderate+exercise_heavy+active_disease_2+Duration+Age+BMI+autoimmune_disease+diabetes), data = b.meta, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (MTX+TNF_inhibitor), data = b.meta, permutations=9999, method = "bray")
b.permanova

##Only RA
b.meta.ra$active_disease_2 <- b.meta.ra$active_disease
b.meta.ra$active_disease_2<-as.character(b.meta.ra$active_disease_2)
b.meta.ra$active_disease_2[which(b.meta.ra$active_disease_2 %in% c("0", "1"))] <- "low"
b.meta.ra$active_disease_2[which(b.meta.ra$active_disease_2 %in% c("2", "3"))] <- "high"

head(b.meta.ra)
class(b.meta.ra$Ds_Duration)
b.meta.ra$Ds_Duration <- as.numeric(as.character(b.meta.ra$Ds_Duration))

b.otu.ra <- otu_table(bac.clean.log.RA)

b.permanova.ra <- adonis(formula = t(b.otu.ra) ~ (RA_factor+anti_CCP+seropositive_RA+GC+MTX+TNF_inhibitor+functionalfood_type+exercise_light+exercise_oderate+exercise_heavy+active_disease_2+Duration+Age+BMI+autoimmune_disease+diabetes), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova.ra

b.permanova.ra <- adonis(formula = t(b.otu.ra) ~ (RA_factor+anti_CCP+BMI+Age+Ds_Duration+active_disease+active_disease_2), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova.ra

b.permanova.ra <- adonis(formula = t(b.otu.ra) ~ (TNF_inhibitor+MTX+MTX_dose), data = b.meta.ra, permutations=9999, method = "bray")
b.permanova.ra


###Only women samples
b.meta.women$Ds_Duration <- as.numeric(as.character(b.meta.women$Ds_Duration))
b.meta.women$Ds_Duration[is.na(b.meta.women$Ds_Duration)] <- 0

b.meta.women$Total_cholesterol

b.otu <- otu_table(bac.clean.log.W)
b.permanova <- adonis(formula = t(b.otu) ~ (Diagnosis+RA_factor+anti_CCP+seropositive_RA+GC+MTX+TNF_inhibitor+exercise_light+exercise_oderate+exercise_heavy+Duration+Age+BMI+autoimmune_disease+diabetes), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (RA_factor+anti_CCP+BMI+Age+Ds_Duration+active_disease), data = b.meta.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (MTX+MTX_dose+TNF_inhibitor), data = b.meta.women, permutations=9999, method = "bray")
b.permanova


### RA and women
b.meta.ra.women$Ds_Duration <- as.numeric(as.character(b.meta.ra.women$Ds_Duration))
b.meta.ra.women$Ds_Duration[is.na(b.meta.ra.women$Ds_Duration)] <- 0

b.otu.ra <- otu_table(bac.clean.log.RA.W)

b.permanova <- adonis(formula = t(b.otu.ra) ~ (RA_factor+Ds_Duration+anti_CCP+seropositive_RA+GC+MTX+TNF_inhibitor+exercise_light+exercise_oderate+exercise_heavy+Duration+Age+BMI+autoimmune_disease+diabetes), data = b.meta.ra.women, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu.ra) ~ (RA_factor+anti_CCP+BMI+Age+Ds_Duration+active_disease), data = b.meta.ra.women, permutations=9999, method = "bray")
b.permanova


b.permanova <- adonis(formula = t(b.otu.ra) ~ (TNF_inhibitor+MTX+MTX_dose), data = b.meta.ra.women, permutations=9999, method = "bray")
b.permanova



### Fungi #####
## remove samples which do not have blood test results
## Lipid
rm1<-f.meta.ra$SampleID[-which(f.meta.ra$Total_cholesterol == 0)]
rm2<-f.meta.ra$SampleID[-which(f.meta.ra$Triglyceride == 0)]
rm3<-f.meta.ra$SampleID[-which(f.meta.ra$HDL == 0)]
rm4<-f.meta.ra$SampleID[-which(f.meta.ra$CRP == 0)]
rm5<-f.meta.ra$SampleID[-which(f.meta.ra$WBC == 0)]
rm6<-f.meta.ra$SampleID[-which(f.meta.ra$ESR == 0)]
rm7<-f.meta.ra$SampleID[-which(f.meta.ra$Hb == 0)]
rm8<-f.meta.ra$SampleID[-which(f.meta.ra$BUN == 0)]
rm9<-f.meta.ra$SampleID[-which(f.meta.ra$Cr == 0)]

f.otu.ra <- otu_table(fun.clean.ss.RA.rel.W)
f.otu.ra <- otu_table(fun.clean.nolog.RA.W)
f.otu.ra <- otu_table(fun.clean.log.RA.W)

### contribution of total cholesterol
f.otu.ra.t <- subset(t(f.otu.ra), rownames(t(f.otu.ra)) %in% as.character(rm1))
f.meta.ra.edit <- subset(f.meta.ra.women, rownames(f.meta.ra.women)%in%rm1)

f.permanova <- adonis(formula = f.otu.ra.t ~ (Total_cholesterol), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova

### contribution of Triglyceride
f.otu.ra.t <- subset(t(f.otu.ra), rownames(t(f.otu.ra)) %in% as.character(rm2))
f.meta.ra.edit <- subset(f.meta.ra.women, rownames(f.meta.ra.women)%in%rm2)

f.permanova <- adonis(formula = f.otu.ra.t ~ (Triglyceride), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova

### contribution of HDL
f.otu.ra.t <- subset(t(f.otu.ra), rownames(t(f.otu.ra)) %in% as.character(rm3))
f.meta.ra.edit <- subset(f.meta.ra.women, rownames(f.meta.ra.women)%in%rm3)

f.permanova <- adonis(formula = f.otu.ra.t ~ (HDL), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = f.otu.ra.t ~ (Total_cholesterol+HDL+Triglyceride), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova

### contribution of CRP
f.otu.ra.t <- subset(t(f.otu.ra), rownames(t(f.otu.ra)) %in% as.character(rm4))
f.meta.ra.edit <- subset(f.meta.ra.women, rownames(f.meta.ra.women)%in%rm4)

f.permanova <- adonis(formula = f.otu.ra.t ~ (CRP), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova

### contribution of wbc
f.otu.ra.t <- subset(t(f.otu.ra), rownames(t(f.otu.ra)) %in% as.character(rm5))
f.meta.ra.edit <- subset(f.meta.ra.women, rownames(f.meta.ra.women)%in%rm5)

f.permanova <- adonis(formula = f.otu.ra.t ~ (WBC), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova

### contribution of Hb
f.otu.ra.t <- subset(t(f.otu.ra), rownames(t(f.otu.ra)) %in% as.character(rm7))
f.meta.ra.edit <- subset(f.meta.ra.women, rownames(f.meta.ra.women)%in%rm7)

f.permanova <- adonis(formula = f.otu.ra.t ~ (Hb), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova


### contribution of ESR
f.otu.ra.t <- subset(t(f.otu.ra), rownames(t(f.otu.ra)) %in% as.character(rm6))
f.meta.ra.edit <- subset(f.meta.ra.women, rownames(f.meta.ra.women)%in%rm6)

f.permanova <- adonis(formula = f.otu.ra.t ~ (ESR), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova


### contribution of Plt
f.otu.ra.t <- subset(t(f.otu.ra), rownames(t(f.otu.ra)) %in% as.character(rm5))
f.meta.ra.edit <- subset(f.meta.ra.women, rownames(f.meta.ra.women)%in%rm5)

f.permanova <- adonis(formula = f.otu.ra.t ~ (Plt), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova

### contribution of BUN
f.otu.ra.t <- subset(t(f.otu.ra), rownames(t(f.otu.ra)) %in% as.character(rm8))
f.meta.ra.edit <- subset(f.meta.ra.women, rownames(f.meta.ra.women)%in%rm8)

f.permanova <- adonis(formula = f.otu.ra.t ~ (BUN), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova


### contribution of Cr
f.otu.ra.t <- subset(t(f.otu.ra), rownames(t(f.otu.ra)) %in% as.character(rm9))
f.meta.ra.edit <- subset(f.meta.ra.women, rownames(f.meta.ra.women)%in%rm9)

f.permanova <- adonis(formula = f.otu.ra.t ~ (Cr), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova


##All
f.otu <- otu_table(fun.clean.ss.rel)
f.otu <- otu_table(fun.clean.nolog)
f.otu <- otu_table(fun.clean.log)

f.otu <- otu_table(fun.its1.clean.log)

f.permanova <- adonis(formula = t(f.otu) ~ (Age+BMI+RA_factor+anti_CCP+Diagnosis+seropositive_RA+GC+treatment2+MTX+TNF_inhibitor+exercise_light+exercise_oderate+exercise_heavy+Duration+autoimmune_disease+diabetes), data = f.meta, permutations=9999, method = "bray")
f.permanova


##Effect of disease
f.permanova <- adonis(formula = t(f.otu) ~ (Diagnosis), data = f.meta, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu) ~ (Duration + Ds_Duration+active_disease), data = f.meta, permutations=9999, method = "bray")
f.permanova

## Quantitative index
f.permanova <- adonis(formula = t(f.otu) ~ (RA_factor+anti_CCP), data = f.meta, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu) ~ (Age+BMI), data = f.meta, permutations=9999, method = "bray")
f.permanova


### Effect of treatment
f.permanova <- adonis(formula = t(f.otu) ~ (MTX+TNF_inhibitor), data = f.meta, permutations=9999, method = "bray")
f.permanova

### Effect of excercise
f.permanova <- adonis(formula = t(f.otu) ~ (exercise_heavy+exercise_oderate+exercise_light), data = f.meta, permutations=9999, method = "bray")
f.permanova


## Effect of predental problem
f.permanova <- adonis(formula = t(f.otu) ~ (periodental_disease+gu_edea_bleeding+sore_on_the_tongue_hollow), data = f.meta, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu) ~ (periodental+periodental_score), data = f.meta, permutations=9999, method = "bray")
f.permanova




##Only RA

f.otu.ra <- otu_table(fun.clean.ss.RA.rel)
f.otu.ra <- otu_table(fun.clean.nolog.RA)
f.otu.ra <- otu_table(fun.clean.log.RA)


f.meta.ra$active_disease_2 <- f.meta.ra$active_disease
f.meta.ra$active_disease_2<-as.character(f.meta.ra$active_disease_2)
f.meta.ra$active_disease_2[which(f.meta.ra$active_disease_2 %in% c("0", "1"))] <- "low"
f.meta.ra$active_disease_2[which(f.meta.ra$active_disease_2 %in% c("2", "3"))] <- "high"

class(f.meta.ra$Ds_Duration)
f.meta.ra$Ds_Duration <- as.numeric(as.character(f.meta.ra$Ds_Duration))

##Effect of disease
f.permanova <- adonis(formula = t(f.otu.ra) ~ (Duration + Ds_Duration+active_disease+active_disease_2), data = f.meta.ra, permutations=9999, method = "bray")
f.permanova

## Quantitative index
f.permanova <- adonis(formula = t(f.otu.ra) ~ (RA_factor+anti_CCP), data = f.meta.ra, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu.ra) ~ (Age+BMI), data = f.meta.ra, permutations=9999, method = "bray")
f.permanova

### Effect of treatment
f.permanova <- adonis(formula = t(f.otu.ra) ~ (MTX+TNF_inhibitor+MTX_dose), data = f.meta.ra, permutations=9999, method = "bray")
f.permanova

### Effect of excercise
f.permanova <- adonis(formula = t(f.otu.ra) ~ (exercise_heavy+exercise_oderate+exercise_light), data = f.meta.ra, permutations=9999, method = "bray")
f.permanova

## Effect of predental problem
f.permanova <- adonis(formula = t(f.otu.ra) ~ (periodental_disease+gu_edea_bleeding+sore_on_the_tongue_hollow), data = f.meta.ra, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu.ra) ~ (periodental+periodental_score), data = f.meta.ra, permutations=9999, method = "bray")
f.permanova




###Only women samples

#Women
f.otu <- otu_table(fun.clean.ss.rel.W)
f.otu <- otu_table(fun.clean.nolog.W)
f.otu <- otu_table(fun.clean.log.W)

##Effect of disease
f.permanova <- adonis(formula = t(f.otu) ~ (Diagnosis), data = f.meta.women, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu) ~ (Duration+active_disease+active_disease_2), data = f.meta.women, permutations=9999, method = "bray")
f.permanova

## Quantitative index
f.permanova <- adonis(formula = t(f.otu) ~ (RA_factor+anti_CCP), data = f.meta.women, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu) ~ (Age+BMI), data = f.meta.women, permutations=9999, method = "bray")
f.permanova


### Effect of treatment
f.permanova <- adonis(formula = t(f.otu) ~ (MTX+TNF_inhibitor), data = f.meta.women, permutations=9999, method = "bray")
f.permanova

### Effect of excercise
f.permanova <- adonis(formula = t(f.otu) ~ (exercise_heavy+exercise_oderate+exercise_light), data = f.meta.women, permutations=9999, method = "bray")
f.permanova


## Effect of predental problem
f.permanova <- adonis(formula = t(f.otu) ~ (periodental_disease+gu_edea_bleeding+sore_on_the_tongue_hollow), data = f.meta.women, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu) ~ (periodental+periodental_score), data = f.meta.women, permutations=9999, method = "bray")
f.permanova


### RA and women

f.otu.ra <- otu_table(fun.clean.ss.RA.rel.W)
f.otu.ra <- otu_table(fun.clean.nolog.RA.W)
f.otu.ra <- otu_table(fun.clean.log.RA.W)


f.permanova <- adonis(formula = t(f.otu.ra) ~ (enstruation), data = f.meta.ra.women, permutations=9999, method = "bray")
f.permanova


##Effect of disease
f.permanova <- adonis(formula = t(f.otu.ra) ~ (Duration + Ds_Duration+active_disease+active_disease_2), data = f.meta.ra.women, permutations=9999, method = "bray")
f.permanova

## Quantitative index
f.permanova <- adonis(formula = t(f.otu.ra) ~ (RA_factor+anti_CCP), data = f.meta.ra.women, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu.ra) ~ (Age+BMI), data = f.meta.ra.women, permutations=9999, method = "bray")
f.permanova

### Effect of treatment
f.permanova <- adonis(formula = t(f.otu.ra) ~ (MTX+TNF_inhibitor+MTX_dose), data = f.meta.ra.women, permutations=9999, method = "bray")
f.permanova

### Effect of excercise
f.permanova <- adonis(formula = t(f.otu.ra) ~ (exercise_heavy+exercise_oderate+exercise_light), data = f.meta.ra.women, permutations=9999, method = "bray")
f.permanova

## Effect of predental problem
f.permanova <- adonis(formula = t(f.otu.ra) ~ (periodental_disease+gu_edea_bleeding+sore_on_the_tongue_hollow), data = f.meta.ra.women, permutations=9999, method = "bray")
f.permanova

f.permanova <- adonis(formula = t(f.otu.ra) ~ (periodental+periodental_score), data = f.meta.ra.women, permutations=9999, method = "bray")
f.permanova


f.otu <- otu_table(fun.clean.log)
f.meta
f.permanova <- adonis(formula = t(f.otu) ~ (Diagnosis+MTX+RA_factor+anti_CCP+seropositive_RA+GC+TNF_inhibitor+exercise_light+exercise_oderate+exercise_heavy+active_disease+Duration+Age+BMI+autoimmune_disease+diabetes), data = f.meta, permutations=9999, method = "bray")
f.permanova

f.otu <- otu_table(fun.clean.log.RA)
f.meta.ra
f.permanova <- adonis(formula = t(f.otu) ~ (MTX+RA_factor+anti_CCP+seropositive_RA+GC+TNF_inhibitor+exercise_light+exercise_oderate+exercise_heavy+active_disease+Duration+Age+BMI+autoimmune_disease+diabetes), data = f.meta.ra, permutations=9999, method = "bray")
f.permanova

f.otu <- otu_table(fun.clean.log.RA.W)
f.meta.ra.edit
f.permanova <- adonis(formula = f.otu.ra.t ~ (Age+BMI+Total_cholesterol+Duration+HDL+Triglyceride+RA_factor+anti_CCP+CRP+ESR), data = f.meta.ra.edit, permutations=9999, method = "bray")
f.permanova

###CAP

f.cap.active <- ordinate(fun.clean.log.RA.W, "CAP", "bray", ~ anti_CCP)

perm_anova.ord <- anova.cca(f.cap.active)
perm_anova.ord2 <- permutest(f.cap.active)

## Plotting
plot.f.cap.active     <- plot_ordination(fun.clean.log.RA.W, f.cap.active, type = "samples", color='anti_CCP', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.active



bray1.bac <-  ordinate(fun.clean.log.RA.W, 'PCoA', 'bray')

#cap_uni <-  ordinate(bac.clean.log, 'CAP', 'unifrac', weight = TRUE, ~ Diagnosis)


#bray3.bac <-  ordinate(bac.clean.log, 'PCoA', 'unifrac')

##Unconstrained PCoA
#Control vs RA
plot_ordination(fun.clean.log.RA.W, bray1.bac, type = "samples", color='anti_CCP', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
