#loading required packages

library(dada2)
library(ShortRead)
library(ggplot2)
library(magrittr)
library(dplyr)
library(vegan)
library(phyloseq)
library(ade4)
library(phangorn)
library(DECIPHER)

#The following code is a workflow modified from the tutorial is this link: 
#https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html 


#Get list of fastq files, 
#enter path for where you kept the files between the quotes
path <- ""
list.files(path)
fnFs <- sort(list.files(path,pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path,pattern = "_R2.fastq.gz", full.names = TRUE))

#extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_R"),'[',1)

#quality Filtering and trimming performed by Erika Dade #
#visualize quality profile
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
#filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#filter with standard parameters except trunQ = 5 for consistency across studies ***
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,170),
#             maxN=0, maxEE=c(2,2), truncQ=5, rm.phix=TRUE,
#            compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
#head(out)

#learn error rates 
#changed from filtFs and filtRs for this run instance for all subsequent 
#instances because Erika ran previous steps
errF <- learnErrors(fnFs, multithread=TRUE)
errR <- learnErrors(fnRs, multithread=TRUE)
plotErrors(errF, nominalQ = TRUE)

#dereplication
derepFs <- derepFastq(fnFs, verbose=TRUE)
derepRs <- derepFastq(fnRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inference
#note: dada(...,pool = "pseudo")performs pseudo-pooling, in which samples are processed 
#independently after sharing information between samples, approximating pooled 
#sample inference in linear time.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#inspect returned dada-class object
dadaFs[[1]]

#merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#remove the deredFs and Rs(118.3Gb) to create space in the R environment
rm(derepFs)
rm(derepRs)




#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#percent non-chimeras
sum(seqtab.nochim)/sum(seqtab)

#Track reads to make sure no step seems faulty
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#assign taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "~/taxa/silva_nr_v132_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "~/taxa/silva_species_assignment_v132.fa")
#check taxa assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#multiple alignments with DECIPHER then create the tree with phangorn RaxML via the ips package within R
seqs <- getSequences(year1seqtab)
# This propagates to the tip labels of the tree
names(seqs) <- seqs 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

#From https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html 
#The phangorn R package is then used to construct a phylogenetic tree. 
#Here we first construct a neighbor-joining tree, and then fit a GTR+G+I 
#(Generalized time-reversible with Gamma rate variation) 
#maximum likelihood tree using the neighbor-joining tree as a starting point.

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))


#extract tree
tree <- fitGTR$tree
tree <- midpoint(tree)

#read in sample data
sampdf <- read.csv("ClinData.csv")
rownames(sampdf) <- sampdf$mblid


#Check row order
sampdf <- sampdf[match(rownames(year1seqtab), sampdf$mblid),]
rownames(sampdf) == rownames(year1seqtab)


#convert to phyloseq object
phylo <- phyloseq(otu_table(year1seqtab, taxa_are_rows=FALSE), 
                  sample_data(sampdf), 
                  tax_table(taxa), phy_tree(tree)) 


#Checking read counts 

rownames(phylo@otu_table@.Data) == rownames(sampdf)
CFSampReads <- rowSums(phylo@otu_table@.Data)[sampdf$Status == "CF"]
HCSampReads <- rowSums(phylo@otu_table@.Data)[sampdf$Status == "CON"]
mean(CFSampReads)#116993.8
sd(CFSampReads)#50311.98
mean(HCSampReads)#109920
sd(HCSampReads)#62345.95

t.test(rowSums(phylo@otu_table@.Data)~sampdf$Status)



