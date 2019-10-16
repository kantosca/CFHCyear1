#########################################################
# Calculate Generalized UniFrac Distances and plot them.# 
#########################################################
library(phyloseq)
library(GUniFrac)
library(ggplot2)

#extracting ASV table
ASV <- as.data.frame(phylo@otu_table@.Data)
#check that ASV matched samples dataframe
rownames(ASV) == rownames(sampdf)

#calculate GUnifracDistances
unifracs <- GUniFrac(ASV, tree, alpha = c(0, 0.5, 1))$unifracs
#extract distance 0.5
d5 <- unifracs[,,"d_0.5"]

#PERMANOVA
adonis(as.dist(d5) ~ sampdf$TimePeriod + sampdf$Status)


#ordination for plotting
pcs <- pcoa(d5, correction = "cailliez")

rownames(sampdf) == rownames(pcs$vectors)
#modify  sample data for plotting
sampdf <- cbind(sampdf, pcs$vectors[,1:3])
pcs$values$Broken_stick[1:2]



#Extract CF data to calculate these distances
CFASV <- ASV[sampdf$Status=="CF",]
CFsamp <- sampdf[sampdf$Status == "CF",]
rownames(CFASV) == rownames(CFsamp)

#plotting 
tiff(file = "fig2.tiff", height = 5, width = 6.875, units = "in",res = 600)
ggplot(revsamp, aes(Axis.1,Axis.2)) + 
  geom_point(size = 1.5,aes( color = Status, shape = Status))  + 
  scale_shape_manual(values =c(17,19), name = "Status", labels = c("Cystic Fibrosis", "Control")) +
  scale_color_manual(values =c("cornflowerblue","darkgreen"),name = "Status", labels = c("Cystic Fibrosis", "Control"))+
  stat_ellipse(aes(color = Status))+ 
  labs(x = "PC 1 (3%)", y = "PC 2 (2%)")+
  theme_bw(base_size = 15) +
  facet_wrap(revsamp$TimePeriod,labeller = as_labeller(c("6W" = "6 Weeks","4M" = "4 Months", "6M" = "6 Months", "9M" = "9 Months", "12M" = "12 Months")) )+
  
  theme(legend.position = c(1.011, 0.25), legend.justification = c(1.1, 0), legend.key.height = unit(2, "mm"))
dev.off()



#running again for CF samples
CFunifracs <- GUniFrac(CFASV, tree, alpha = c(0, 0.5, 1))$unifracs
#extract distance 0.5
CFd5 <- CFunifracs[,,"d_0.5"]


adonis(as.dist(CFd5) ~ CFsamp$TimePeriod + CFsamp$PulmEx1yr+ 
         CFsamp$PancIns + CFsamp$dFeedingModeED + 
         CFsamp$DeliveryMode + CFsamp$B4Abx)

#Extracting PCs and percent explained for plotting
CFpcs <- pcoa(CFd5, correction = "cailliez")
CFpc1_2 <- CFpcs$vectors[,1:2]
colnames(CFpc1_2) <- c("PC1","PC2")
CFsamp <- cbind(CFsamp,CFpc1_2)
CFpcs$values$Broken_stick[1:2]

#plotting
tiff(file = "fig3.tiff", height = 5, width = 6.875, units = "in", res = 800)
ggplot(CFsamp, aes(PC1,PC2, shape = PulmEx1yr, color = PulmEx1yr)) + 
  geom_point(size = 3) +
  labs(x = "PC 1 (7%)", y = "PC 2 (5%)") +
  theme_bw(base_size = 15) + 
  scale_shape_manual(values =c(17,19), name = "Exacerbation", labels = c("No", "Yes")) +
  scale_color_manual(values =Blues[2:5], name = "Exacerbation", labels = c("No", "Yes"))+
  facet_wrap(CFsamp$TimePeriod, labeller = as_labeller(c("6W" = "6 Weeks","4M" = "4 Months", "6M" = "6 Months", "9M" = "9 Months", "12M" = "12 Months")))+
  theme(legend.position = c(1, 0.1), legend.justification = c(1.1, 0))
dev.off()

set.seed(33)
#repeating with Bray-Curtis distances
#calculatBray-curtis distance
bray <- vegdist(ASV, method = "bray")
CFbray <- vegdist(CFASV, method = "bray")
rownames(CFASV) == CFsamp$mblid
for (i in levels(CFsamp$TimePeriod)){
  CFbrayTime <- vegdist(CFASV[CFsamp$TimePeriod == i,], method = "bray")
  print(i)
  print(adonis(as.dist(CFbrayTime) ~ CFsamp$PulmEx1yr[CFsamp$TimePeriod == i]))
}

adonis(bray~sampdf$TimePeriod + sampdf$Status)
adonis(CFbray~CFsamp$TimePeriod + CFsamp$PulmEx1yr)
adonis(CFbray~CFsamp$TimePeriod)
for (i in levels(CFsamp$TimePeriod)){
  print(i)
  print(adonis(as.dist(CFbray[CFsamp$TimePeriod == i ,CFsamp$TimePeriod == i]) ~ CFsamp$PancIns[CFsamp$TimePeriod == i]))
}


#setting up to plot Bray-Curtis PCoA
CFbraypcs <- pcoa(CFbray, correction = "cailliez")
CFbraypc1_2 <- CFbraypcs$vectors[,1:2]
colnames(CFbraypc1_2) <- c("BrayPC1","BrayPC2")
CFsamp <- cbind(CFsamp,CFbraypc1_2)
colnames(CFsamp)<- c(colnames(CFsamp)[1:18],"BrayPC1","BrayPC2")
CFbraypcs$values$Broken_stick[1:2]





#plotting by feeding

tiff(file = "SupFig3A.tiff", height = 5, width = 6.875, units = "in",res = 600)
ggplot(CFsamp, aes(PC1,PC2, color = TimePeriod, shape = as.factor(dFeedingModeED))) + 
  scale_shape_manual(values = c(15,17,19),name = "Feeding Mode", labels = c("Breast", "Formula", "Mixed"))+ 
  geom_point(size = 5) +
  labs(x = "PC 1", y = "PC 2")+
  theme_bw(base_size = 12) +
  scale_color_manual(values = AgeCol, name = "Age")
dev.off()







