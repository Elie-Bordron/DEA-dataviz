#!/usr/bin/Rscript
# csv_formatting.R
# 
# create column "Full Location" from column "Microarray.Nomenclature..ISCN.2013."
# 
#
# Author: Elie Bordron
# Date: 2022/02/22

pathDataDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/data/working_data/2-AD"
celFile = "2-ADREC.RC.OSCHP.segments.csv"
pathSample = file.path(pathDataDir, celFile)
sample2AD = read.csv(pathSample, header = TRUE, sep=';')

chromosomeCol = sample2AD$Chromosome
microarrNomenclCol = sample2AD[3]

newColAsVec = c()
for (i in 1:nrow(microarrNomenclCol)) {
    row = microarrNomenclCol[i, ]
    rowSliced = stringr::str_split(row,'\\(')
    vecstr = unlist(rowSliced)
    secondSplit = unlist(stringr::str_split(vecstr[2], '\\)'))
    full_location = secondSplit[1]
    chrm = chromosomeCol[i]
    rowOfNewCol=paste("chr", chrm, ':', full_location, sep='')
    newColAsVec = append(newColAsVec, rowOfNewCol)
}
sample2AD["Full Location"] = newColAsVec
# sample2AD["Full location "]
celFile_csv = paste(gsub('.{4}$', '', celFile), "_FULL_LOCATION", ".csv", sep='')
celFile_txt = paste(gsub('.{4}$', '', celFile), "_FULL_LOCATION", ".txt", sep='')
path_out = file.path(pathDataDir, celFile_txt)
print(c("path_out: ", path_out))
# write.csv2(sample2AD, path_out)
colCN = which(colnames(sample2AD)=="CN.State")
colnames(sample2AD)[colCN] = "CN State"
sample2AD <- sample2AD[, c(
"Full Location", "File", "Chromosome", "Microarray.Nomenclature..ISCN.2013.",
"Cytoband.Start", "CN State", "Median.Log2Ratio",  
"Type", "Size..kbp.", "Genes", 
"Marker.Count", "Gene.Count", "OMIM...Genes.Count",  
"OMIM...Genes", "CytoRegions", "Call",  
"Segment.Interpretation", "Curation.By", "Materially.Modified.Segment", 
"Mean.Weighted.Log2Ratio"
)]
write.table(sample2AD, path_out, sep='\t', quote=FALSE, row.names=FALSE)
colnames(sample2AD) 


