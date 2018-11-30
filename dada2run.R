

#Get list of fastq files
path <- "~//Google Drive/Kathys Schoolwork/Dartmouth/CFHCyear1/samples"
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





