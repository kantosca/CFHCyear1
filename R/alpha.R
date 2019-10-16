#########################################
# Calculate alpha diversity and plot it.# 
#########################################

library(phyloseq)
library(ggplot2)

set.seed(33)
#calculate alpha diversity
alphadiv <- estimate_richness(phylo)

#remove X from row names
rownames(alphadiv)<-gsub("X","",rownames(alphadiv))

#extract shannon diversity and add it to sampdf the 
#clinical covariate dataframe for ease of plotting
Shannon <- alphadiv$Shannon
Simpson <- alphadiv$Simpson
Fisher <- alphadiv$Fisher
rownames(alphadiv)==rownames(sampdf)
sampdf <- cbind(sampdf, Shannon, Simpson)
sampdf <- cbind(sampdf, Fisher)


#Set order for plotting by age and Status
sampdf$TimePeriod <- factor(sampdf$TimePeriod, levels = c("6W","4M","6M","9M","12M"))
sampdf$Status <- factor(sampdf$Status, levels = c("HC","CF"))

#Generate Colors for plots 
Greens <- c("darkolivegreen1","chartreuse","springgreen3","chartreuse4","#003300")
Blues <- c("#00FFFF","deepskyblue", "royalblue1","slateblue3","#0000FF")
Pinks <- c("#FF99FF","#FF66CC","#FF3399","violetred","deeppink4")
Purples <- c("#CC99FF","mediumpurple2","mediumorchid","mediumorchid4","#660099")
Reds <- c("#FF3333","firebrick2","#CC0033","#990000","indianred4")
Oranges <- c("#FFFF00","#FFCC00","#FF9900","#FF6600","darkorange4")
GreenBlues <- as.vector(rbind(Greens, Blues))
brewer.pal(5,"Greens")

plot(1:10,1:10,col =Blues ,pch = 19,cex = 2)
plot(1:30,1:30,col = Rainbow,pch = 19,cex = 2)
Rainbow <- as.vector(cbind(Pinks[1:4], Reds[1:4], Oranges[1:4], Greens[1:4], Blues[1:4], Purples[1:4]))
Rainbow2 <- as.vector(cbind(Blues[1:4], Greens[1:4], Oranges[1:4], Pinks[1:4], Purples[1:4]))
Rainbow2[4] <- Blues[5]

#wilcoxon rank sum on Shannon diversity with multiple testing correction
pvals <- vector()
for (i in levels(sampdf$TimePeriod)){
  print(i)
  pvals[i] <- print(wilcox.test(data = sampdf[sampdf$TimePeriod==i,],Shannon~Status))$p.value
}

p.adjust(pvals, method = "fdr", n=5)

#plotting 
BlueGreen5 <- c(rep("darkgreen",5),rep("cornflowerblue",5)) 
tiff(file = "fig1.tiff",height = 4, width = 6.875,units = "in", res = 600)
ggplot(data = sampdf, aes(x = TimePeriod, y = Shannon))+
  facet_wrap(facets=sampdf$Status,labeller = as_labeller(c("CON" = "Control", "CF" = "Cystic Fibrosis")))+ 
  geom_boxplot(fill = c(BlueGreen5) ) +
  theme_bw(base_size = 15) +
  xlab("Age (Weeks/Months)")+
  ylab("Shannon Diversity Index")
dev.off()




### Fisher and rarefaction ###
set.seed(100)
#dropped 4 low count samples from control data and 
#1343 ASVs because they are no longer present in 
#any sample after random subsampling
alpharare <- rarefy_even_depth(phylo, sample.size = 9500)
rowSums(alpharare@otu_table@.Data)#confirm all samples have 9500 read counts
alphaR <- estimate_richness(alpharare)

RareDF <- as.data.frame(alpharare@sam_data@.Data)
colnames(RareDF)<- alpharare@sam_data@names
rownames(RareDF)<- alpharare@sam_data@row.names
rownames(alphaR) <- gsub("X","",rownames(alphaR))
rownames(alphaR) == rownames(RareDF)
RareDF <- cbind(RareDF,alphaR)

#wilcoxon rank sum on Shannon diversity post rarefaction with multiple testing correction
pvalsRare<- vector()
for (i in levels(RareDF$TimePeriod)){
  print(i)
  pvalsRare[i]<-print(wilcox.test(RareDF$Shannon[RareDF$TimePeriod == i]~RareDF$Status[RareDF$TimePeriod == i]))$p.value
}
p.adjust(pvalsRare, method = "fdr", n=5)

#plotting SDI post rarefaction
tiff(file = "ShannonRarefy.tiff",height = 4, width = 6.875,units = "in", res = 600)
ggplot(data = RareDF, aes(x = TimePeriod, y = Shannon))+
  facet_wrap(facets=RareDF$Status,labeller = as_labeller(c("CON" = "Control", "CF" = "Cystic Fibrosis")))+ 
  geom_boxplot(fill = c(BlueGreen5) ) +
  theme_bw(base_size = 15) +
  xlab("Age(Months/Years")+
  ylab("SDI on Rarefied data")
dev.off()


#wilcoxon rank sum on Fisher diversity with multiple testing correction
pvalsFisher<-vector()
for (i in levels(sampdf$TimePeriod)){
  print(i)
  print(wilcox.test(sampdf$Fisher[sampdf$TimePeriod == i]~sampdf$Status[sampdf$TimePeriod == i]))
}
p.adjust(pvalsFisher, method = "fdr", n=5)

#plotting Fisher's alpha 


tiff(file = "FisherAll.tiff",height = 4, width = 6.875,units = "in", res = 600)
ggplot(data = sampdf, aes(x = TimePeriod, y = Fisher))+
  facet_wrap(facets=sampdf$Status,labeller = as_labeller(c("CON" = "Control", "CF" = "Cystic Fibrosis")))+ 
  geom_boxplot(fill = c(BlueGreen5) ) +
  theme_bw(base_size = 15) +
  xlab("Age(Months/Years")+
  ylab("Fisher's alpha")
dev.off()


summary(glm(sampdf$Shannon~sampdf$TimePeriod+ sampdf$Status+
              sampdf$dFeedingModeED + PancAll + PulmExAll +B4Abx))
