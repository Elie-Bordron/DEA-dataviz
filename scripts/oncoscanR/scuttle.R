setwd("C:\\Users\\e.bordron\\Desktop\\CGH-scoring\\M2_internship_Bergonie") #ctrl maj H
library(devtools)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
BiocManager::install("GenomicRanges")
library(oncoscanR)
segs.filename <- system.file("extdata", "chas_example.txt", package = "oncoscanR")
workflow_oncoscan.run(segs.filename, "M")
