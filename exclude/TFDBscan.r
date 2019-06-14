# TFDBscan
# http://bioconductor.org/packages/release/bioc/vignettes/TFBSTools/inst/doc/TFBSTools.html


# retrieve 2018 JASPAR database
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("JASPAR2018", version = "3.8")
# install TFBSTools
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("TFBSTools", version = "3.8")

library(JASPAR2018)
library(TFBSTools)

# get plant PFM
opts <- list()
opts[["tax_group"]] <- "plants"
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
PFMatrixList
# PFMatrixList of length 489
# names(489): MA0020.1 MA0021.1 MA0034.1 MA0044.1 ... MA1085.2 MA1416.1 MA1417.1

# matrix conversion
pfm = PFMatrixList[[1]]
pfm
pwm <- toPWM(pfm, pseudocounts=0.8) # default pseudocounts=0.8, necessary for correcting the small number of counts or eliminating the zero values before log transformation. In TFBS perl module, the square root of the number of sequences contributing to each column. However, it has been shown to too harsh (Nishida, Frith, and Nakai 2009). Hence, a default value of 0.8 is used. Of course, it can be changed to other customized value or even different values for each column.
pwm
icm <- toICM(pfm, pseudocounts=0.8, schneider=TRUE)
icm
seqLogo(icm)

# scans a nucleotide sequence with the pattern represented in the PWM
library(Biostrings)
pwmList <- do.call(PWMatrixList, lapply(PFMatrixList,toPWM))
subject <- DNAString("GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG")
sitesetList <- searchSeq(pwmList, subject, seqname="seq1", min.score="60%", strand="*")
# write gff
head(writeGFF3(sitesetList))
## get the relative scores
relScore(sitesetList)
#> $MA0003.2
#> [1] 0.6196884
#>
#> $MA0004.1
#> [1] 0.6652185 0.6652185 0.6141340 0.6633668 0.6141340 0.6633668
## calculate the empirical p-values of the scores
pvalues(sitesetList, type="TFMPvalue")
#> [1] 0.02734375 0.02734375 0.04638672 0.04052734 0.04638672 0.04052734
pvalues(siteset, type="sampling")
#> [1] 0.0258 0.0258 0.0581 0.0396 0.0581 0.0396
## generate gff2 or gff3 style output
head(writeGFF3(siteset))

