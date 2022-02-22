## Dada2 and phyloseq demo for Ganda lab
# Emily Bean, 11/2020

# ---- install and read data ----

### dada2 and phyloseq are both Bioconductor packages

# test if Bioconductor is already installed
install.packages("BiocManager", version = 3.14)

# load package
require(BiocManager)

# install dada2
BiocManager::install("dada2")

# load package
require(dada2)


# path to directory that contains fastq files
# Github doesn't handle file storage well so these are on OneDrive
# PATH = "~/The Pennsylvania State University/Ganda, Erika - Shared-Trello-Projects/Mouse_Opioid/Mouse_Opioid_Reads/"
PATH <- "/Users/anafonseca/bioinfo_trainning/"

# test this variable
PATH

# Let's see if anything's in this folder
ourlist <- list.files(PATH)

# How many fastq files are in my folder
length(ourlist)

#  set pathway to the database you prefer - we will use Greengenes
## NOTE:Github doesn't handle fastq.gz files well so download these locally
# dada2 maintains a list of databases: https://benjjneb.github.io/dada2/training.html
DB = "~/Desktop/gg_13_8_train_set_97.fa.gz"


# paired end characterization; most Illumina files are sample names + "_R1_001.fastq" for forward reads
# however, sequences downloaded from NCBI have patterns: "_1.fastq" for forward and "_2.fastq" for reverse
PATTERNF = "_R1_001.fastq.gz"
PATTERNR = "_R2_001.fastq.gz"


## ---- getFiles----

# get forward and reverse files
## for this demo, subset to the first 2 forward and 2 reverse files
fwdFiles <- list.files(PATH, pattern = PATTERNF, full.names = TRUE)[1:2]
revFiles <- list.files(PATH, pattern = PATTERNR, full.names = TRUE)[1:2]

# check to make sure that the lengths of both files are the same
if(length(fwdFiles) != length(revFiles)) {
  
  stop("There is an unequal number of forward and reverse files")
}

# get sample names
fwdNames <- sapply(strsplit(basename(fwdFiles), PATTERNF), `[`, 1)
revNames <- sapply(strsplit(basename(revFiles), PATTERNR), `[`, 1)

## NOTE: DEFAULT CODE ASSUMES FWD AND REV FILES ARE ORDERED
# error catch if unordered
if(any(!fwdNames %in% revNames)) {
  
  stop("forward and reverse files are out of order")
  
}

### ---- filterAndTrim ----

# create subdirectory for filtered files
filtForward <- file.path("./data/filtered", paste0(fwdNames, "_F_filt.fastq.gz"))
filtReverse <- file.path("./data/filtered", paste0(revNames, "_R_filt.fastq.gz"))

## Dada2 can plot Phred qualities but FastQC is much better
plotQualityProfile(fwdFiles[1])

# filter and trim
cleaned <- filterAndTrim(
  # set forward and reverse paths
  fwd = fwdFiles, rev = revFiles,
  # set paths for the filtered files that will be created
  filt = filtForward, filt.rev = filtReverse,
  # add any necessary filtering parameters
  maxN = 0,
  
  # MAC ONLY: multithread ability
  multithread = TRUE
)

# path to filtered and cleaned reads
CLEANEDPATH = "./data/filtered/"

# visualize the quality of filtered data on a forward & reverse read
# in reality: do this step in FastQC as well
plotQualityProfile(CLEANEDPATH[1])



## ---- dada2 algorithm----

### NOTE: some patterns are re-done because this script can also be run
# after filtering with Trimmomatic or another program like QIIME

# pattern that specifies which reads are forward or reverse
# if single-read pairs, only specifyforward 
FILTEREDF = "_F_filt.fastq.gz"
FILTEREDR = "_R_filt.fastq.gz"

# get forward and reverse reads
forward <- sort(list.files(CLEANEDPATH, pattern = FILTEREDF, full.names = TRUE))
reverse <- sort(list.files(CLEANEDPATH, pattern = FILTEREDR, full.names = TRUE))

# check to make sure that the lengths of both files are the same and that they match
fwdNames <- sapply(strsplit(basename(forward), FILTEREDF), `[`, 1)
revNames <- sapply(strsplit(basename(reverse), FILTEREDR), `[`, 1)
if(length(fwdNames) != length(revNames)) {
  stop("The number of forward and reverse files do not match.")
} else {
  
  if(any(!fwdNames%in% revNames)) {
    
    stop("Forward and reverse reads are out of order.")
  }
}

# perform error learning
errF <- learnErrors(forward, 
                    multithread = TRUE)
errR <- learnErrors(reverse, 
                    multithread = TRUE)

# perform denoising on forward and reverse reads
dadaForward <- dada(derep = forward, 
                    err = errF, 
                    multithread = TRUE)
dadaReverse <- dada(derep = reverse, 
                    err = errR, 
                    multithread = TRUE)

# merge paired reads
mergers <- mergePairs(dadaF = dadaForward,
                      derepF = forward,
                      dadaR = dadaReverse,
                      derepR = reverse,
                      verbose = TRUE)

# construct sequence table of ASVs
seqtab <- makeSequenceTable(samples = mergers)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(unqs = seqtab, 
                                    method = "consensus",
                                    multithread = TRUE,
                                    verbose = TRUE)

# assign taxonomy using the Silva database
tax <- assignTaxonomy(seqs = seqtab.nochim, 
                      refFasta = DB, 
                      multithread = TRUE,
                      verbose = TRUE)

## This is the end of the dada2 algorithm
## OPTIONAL: export ASV table and taxonomy table to text files
# to import into QIIME or other program


### ---- phyloseq ----

# create demo sample data - in reality, you would have metadata
# associated with your dataset
sampleDF <- data.frame(sample = sapply(strsplit(rownames(seqtab.nochim), "_F"), `[`, 1),
                       treatment = c("Treatment", "Control"),
                       sex = c("Male", "Female"),
                       row.names = rownames(seqtab.nochim))

# make phyloseq object
## NOTE: Phyloseq has wrappers for import directly from QIIME among others
# check phyloseq Github for more info

## Phyloseq objects need 3 things: ASV table, taxonomy table, and sample data
# OPTIONAL: add a phylogenetic tree
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(sampleDF), 
               tax_table(tax))

# default settings are messy; replace strings of DNA with an ASV sequence 
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# examine the phyloseq object
ps

# to examine parts of phyloseq object, use "@"
ps@sam_data

# to get columns within the sample data, use "@" followed by "$"
ps@sam_data$treatment

# what does the ASV table look like?
head(ps@otu_table)

# the taxonomy table serves as a lookup table for ASVs
head(ps@tax_table)

# how many Phyla were classified?
get_taxa_unique(physeq = ps,
                taxonomic.rank = "Phylum")

## NOTE: phyloseq plots are wrappers for ggplot and take gg commands
#Abundance barplot by Sample
plot_bar(physeq = ps, 
         x = "sample", 
         y = "Abundance", 
         fill = "Phylum", 
         title = "Abundance Barplot by Sample")

# Alpha diversity plot
plot_richness(ps, x="treatment", measures=c("Shannon", "Simpson"), color="treatment")



