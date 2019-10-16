####################################
# Identify differences by taxonomy #
####################################
library(phyloseq)
library(ggplot2)
library(dplyr)



#subset data for binary taxa comparisons. 
#Starting with broader criteria here so that we can 
#identify taxa for further analyses (controlling for covariates) 
#and wetlab followup. 
phylo6W <-phylo %>% subset_samples(TimePeriod %in% "6W") 
phylo4M <-phylo %>% subset_samples(TimePeriod %in% "4M") 
phylo6M <-phylo %>% subset_samples(TimePeriod %in% "6M") 
phylo9M <-phylo %>% subset_samples(TimePeriod %in% "9M") 
phylo12M <-phylo %>% subset_samples(TimePeriod %in% "12M") 
phyloCF <-phylo %>% subset_samples(Status %in% "CF") 
multcomtaxa6W <- mt(phylo6W, "Status")
multcomtaxa4M <- mt(phylo4M, "Status")
multcomtaxa6M <- mt(phylo6M, "Status")
multcomtaxa9M <- mt(phylo9M, "Status")
multcomtaxa12M <- mt(phylo12M, "Status")
multcomtaxaPulmEx <- mt(phyloCFgenus, "PulmEx1yr")
multcomtaxaPanc <- mt(phyloCF, "PancIns")


#collapse phyloseq object at the genus and phylum level
phyloGenus <- tax_glom(phylo,taxrank = "Genus", NArm = FALSE)
phyloPh <- tax_glom(phylo,taxrank = "Phylum", NArm = FALSE)
phyloCFgenus <- tax_glom(phyloCF, taxrank = "Genus", NArm = FALSE)

#testing for differences in taxa level
multcompGenus <- mt(phyloGenus, "Status")
multcompPhylum <- mt(phyloPh, "Status")
mtCFPulmExGenus <- mt(phyloCFgenus, "PulmEx1yr")
write.csv(mtCFPulmExGenus, file = "CFPulmExGenus.csv")
mtCFPulmExStrep <- mt(CFstrep, "PulmEx1yr")
write.csv(mtCFPulmExStrep, file = "CFPulmExStrep.csv")

#transform for plotting relative abundance
RAphylo <- transform_sample_counts(phyloGenus, function(x) x / sum(x))
rowSums(RAphylo@otu_table@.Data)
mergeStat <- merge_samples(RAphylo, sampdf$Status)
mergeStat<- transform_sample_counts(mergeStat, function(x) x / sum(x))
mergeStat20 <- prune_taxa(Top20,mergeStat)
rowSums(mergeStat20@otu_table@.Data)

genNames <- (A20mergeRA@tax_table@.Data[,6])
vals <- A20mergeRA@otu_table@.Data
colnames(vals)<-genNames

rowSums(phylo20@otu_table@.Data)
#top 20 for bar plot
Top20 <- names(sort(taxa_sums(phyloGenus), TRUE)[1:20])
phylo20 <- prune_taxa(Top20, phyloGenus)

phylo20StatAge <- merge_samples(phylo20, sampdf$StatusAge)
phylo20RA <- transform_sample_counts(phylo20, function(x) x / sum(x))
phylo20StatAgeRA <- transform_sample_counts(phylo20StatAge, function(x) x / sum(x))
rowSums(phylo20StatAgeRA@otu_table@.Data)
rowSums(phylo20RA@otu_table@.Data)


colors = c( "#D896FF", "#E30074", "#FF4E50", "#FC913A", "#FFAAA5", "#FFD3B6", "#EAE374", "#F9D62E",  
            "#8AE429", "#E2F4C7", "#9AFE2E", "#FE00F6", "#77AAFF", "#99CCFF", "#FF9900",
            "#A8E6CF", "#00F9FF", "#005EFF", "#009FFF", "#0900FF", "#00308F", "#B8D000", "#4D7F17", "#6BB120",
            "#BBEEFF", "#800080", "#BE29EC", "#EFBBFF", "#FF8B94", "#0392CF", "#FDF498", "#FFD4E5", "#D4FFEA",
            "#EECBFF", "#DBDCFF", "#FE0000", "#E08822", "#5E8EB7", "#BD5757", "#EE4035", 
            "#73BA9F", "#FDFE02", "#9BA2FF", "#D4CAC5")


#colors palette to make figures accessible for color blind readers 
#adapted from: https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
#color differentiation tested vis color oracle: https://colororacle.org

cbcols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
            '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', 
            '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', 
            '#000075', '#808080', "#ff1b00", "#ff929b", "#0092ff", "#00ff00", 
            "#d0165e", "#b7cb48", "#af4ddb", "#00dfe3", "#008200", "#ec8efb", 
            "#beff64", "#0000f0", "#00b65b", "#a1027e")
cbcols2 <- c('#f9194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
             '#46f0f0', '#f032e6', "#d0165e", '#bcf60c', '#fabebe', '#008080', 
             '#e6beff',  '#800000', '#aaffc3', '#ffd8b1', '#000075',  "#0092ff", 
             "#00ff00", "#b7cb48")

#my attempt at at colorblind friendly palette
cbtest<- c("#ff0000", "#ffce00", "#0be567", "#005eff", "#cb06ff",
           "#fd8c8c","#e7f695", "#a5f69d", "#90e1f9", "#d99efb", 
           "#e6194b", "#ffe119", "#3cb44b","#4363d8", "#911eb4",
           "#f58231", "#fffac8","#526b2d", "#0d034d","#9600ff")


#Plot relative abundance barplot
tiff("fig4.tiff", height = 4, width = 6.875, units= "in", res = 600)
ggplot(BPtest2, aes(x = Status, y = Abundance, fill = Genus, group = Abundance)) + geom_bar(stat = 'identity') +
  theme_bw(base_size =12)+ 
  theme(legend.position = "bottom",legend.text = element_text(face = "italic",size = 7.5),legend.key.size = unit(3,"mm"),legend.title = element_text(size=8)) + 
  scale_fill_manual(values = cbtest) +xlab("Age and Health Status") +ylab("Relative Abundance %")+
  scale_y_continuous(labels = c("0","25","50","75","100"))+
  annotate(geom = "text" , data = labs,x = 1:10, y=ypos ,label= "a", inherit.aes = FALSE)#+
scale_x_discrete(labels = c(" CF" = "CF", "Con" = "Control"))+
  facet_wrap(BPtest2$Age, labeller = as_labeller(c("6W" = "6 Weeks","4M" = "4 Months", "6M" = "6 Months", "9M" = "9 Months", "12M" = "12 Months")),ncol = 5)#+
annotate(geom = "text" , data = labs,x = 1:10, y=ypos ,label= "a", inherit.aes = FALSE)
dev.off()


###############################################################################
# Dot plot showing taxa spread

#prep data
A20merge <- merge_samples(phylo20RA,sampdf$Status, fun = mean)
rowSums(A20mergeRA@otu_table@.Data)
colSums(A20mergeRA@otu_table@.Data)
A20mergeRA <- transform_sample_counts(A20merge, function(x) x / sum(x))


#Find difference in means
meandif <- A20mergeRA@otu_table@.Data[2,] - A20mergeRA@otu_table@.Data[1,]

Gen20 <- phylo20RA@tax_table[,6]
#replace NA with actual taxa
Gen20[4]<-paste("o_",phylo20RA@tax_table[4,5], "_g_NA",sep="")

Gen20<- Gen20[order(meandif)]
Gen20

#adding significance marks to labels
Gen20 <- gsub("Bacteroides","***Bacteroides", Gen20)
Gen20 <- gsub("Roseburia","***Roseburia", Gen20)
Gen20 <- gsub("Veillonella","*Veillonella", Gen20)

redlabs = c("***Bacteroides","***Roseburia","*Veillonella")
colorlist = c("black","red")


axiscolor = colorlist[Gen20 %in% redlabs +1]

factor1 <- rownames(sampdf[sampdf$Status == "CF",])
factor2 <- rownames(sampdf[sampdf$Status == "CON",])



#Relative abundance of all samples plot
tiff(file = "fig5.tiff", width = 6.875, height = 6, res = 600, units = "in")
par(mar=c(5,10,2,1))
plot(0,0, xlim = c(-1,1), ylim = c(1,20), xlab = "Relative abundance %", ylab = "", yaxt = 'n', xaxt = 'n',cex.lab = 1, cex = 0.01)
title(main = "Top taxa differences", cex.main = 1)
grid(ny = 20)
par(xpd = TRUE)
legend(-1.98,-0.8, legend = c ("Sample Relative Abundances", "Mean Relative Abundance","Difference in Means"), col = c("black","steelblue","red"), pch = 16, title = "Legend", pt.cex=1, cex=0.75)

for (i in 1:20){
  #plot relative abundances
  j = names(sort(meandif))[i]
  points(x = -phylo20RA@otu_table@.Data[factor = factor1,j],rep(i,length(factor1)), type = 'p', pch = 16, col = "black",cex = 0.75)
  points(x = phylo20RA@otu_table@.Data[factor = factor2,j],rep(i,length(factor2)), type = 'p', pch = 16, col = "black",cex = 0.75)
  
}

#plot means and difference in means
points(x = -A20mergeRA@otu_table@.Data[1,order(meandif)], y = c(1:20), xlim = c(-1,1),xlab = "", ylab = "", pch = 16, col = 'steelblue1',
       yaxt = 'n', xaxt = 'n',cex = 1)
points(x = A20mergeRA@otu_table@.Data[2,order(meandif)], y = c(1:20), xlim = c(-1,1),xlab = "", ylab = "", pch = 16, col = 'steelblue1',
       yaxt = 'n', xaxt = 'n',cex = 1)
points(x = sort(meandif), y = c(1:20), xlim = c(-1,1),xlab = "Relative abundance%",ylab = "", pch = 16, col = 'red',
       yaxt = 'n', xaxt = 'n', cex = 1.25)
#add axis labels
axis(2, at = c(1:20), labels = FALSE, las = 2,cex.axis =1, font = 3)
text(labels = Gen20, col = axiscolor,  x=rep(-1.1,length(Gen20)), y=1:length(Gen20),pos = 2, font = 3, xpd = TRUE, cex = 0.875)
axis(1,at = c(-1,-0.5,0,0.5,1), labels = c("100%","CF","0%","Control","100%"),cex.axis = 1)
#lines(x = rep(0,27), y = c(0:26))

lines(x = c(0,0),y = c(0.25,20.75))

dev.off()


#########################################################################


#extracting Bacteroides relative abundances for ease of plotting
BacteroidesRA <- phylo20RA@otu_table@.Data
colnames(BacteroidesRA) <- phylo20RA@tax_table@.Data[,6] 
#Double checked that these rows match prior to this step
BacteroidesRAdat <- cbind(BacteroidesRA[,"Bacteroides"], sampdf)
colnames(BacteroidesRAdat) <- c("Bacteroides", colnames(sampdf))


#Plot Bacteroides relative abundance


tiff("fig6.tiff", height = 3.67, width = 6.875, units= "in", res = 600)
ggplot(BacteroidesRAdat, aes(x = TimePeriod, y = Bacteroides, fill = StatusAge))+
  facet_wrap(facets=sampdf$Status,labeller = as_labeller(c("CON" = "Control", "CF" = "Cystic Fibrosis")))+ 
  geom_boxplot(fill = c(BlueGreen5) ) +
  theme_bw(base_size = 15) +
  xlab("Age (Weeks/Months)")+
  scale_y_continuous(labels = c("0","25","50","75","100"),limits = c(0,1))+
  theme_bw(base_size =12)+
  ylab("Bacteroides Relative Abundance %")
dev.off()

 
