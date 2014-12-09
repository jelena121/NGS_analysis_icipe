#loading the library
library("DESeq") 

# reading in counts data from a file
rawData <- read.table("counts_table.txt", header=T, row.names=1, sep="\t")

head(rawData)
names(rawData) <- c("WT_a", "WT_b","WT_c","WT_d","patient_a",
                 "patient_b","patient_c","patient_d")

# specifying conditions
conds <- c(rep("WT",4), rep("patient",4))

# creating counts object
cds <- newCountDataSet(rawData, conds)

# estimating size factors
cds <- estimateSizeFactors( cds )

# viewing size factors
sizeFactors( cds )

# saving normalized counts in a variable
nCounts <- counts(cds, normalized=TRUE)

# writing out the normalized counts to a table
write.table(nCounts, file="normalized_counts.txt", sep="\t", col.names=NA, quote=F)

# estimating dispersions
cds = estimateDispersions( cds )

# looking up model information
str( fitInfo(cds) )

# plotting dispersions
plotDispEsts( cds )

# calculating results
res = nbinomTest( cds, "WT", "patient" )

head(res)

# filtering the results table

# setting cutoffs
pvalCutoff <- 0.01
readsCutoff <- 10
foldCutoff <- 1

# filtering infinities and incomplete results
resFiltered <- res[complete.cases(res) & res$baseMeanA != 0 & res$baseMeanB != 0 & 
                       res$log2FoldChange!="-Inf" 
                     & res$log2FoldChange!="Inf",]

# only including results above a certain read count
resFiltered10Reads <-resFiltered[resFiltered$baseMeanA > readsCutoff & 
                                   resFiltered$baseMeanB > readsCutoff,]

# applying p-value cutoff

resSig <- resFiltered10Reads[ resFiltered10Reads$padj < pvalCutoff, ]; 
head(resSig)

# filtering for magnitude of expression change
resFinal <- resSig[ abs(resSig$log2FoldChange) > foldCutoff, ]

# renaming individual columns
names(resFinal)[3] <- "meanWT"
names(resFinal)[4] <- "meanPatient"
names(resFinal)[6] <- "log2FoldChange (patient/WT)"

head(resFinal)

# writing out results to file

write.table(resFinal, file="patient_vs_WT_p0_01.txt", sep="\t", 
            row.names=TRUE, col.names=NA)
