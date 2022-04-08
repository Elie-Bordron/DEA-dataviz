require(EaCoN)
library(EaCoN)
# to load function custom_OS.Process()
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/EaCoN_functions.R")








## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)



pathToData = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_cel_files"
pathToATCelFile = file.path(pathToData, "5-LD_AT_(OncoScan_CNV).CEL")
pathToGCCelFile = file.path(pathToData, "5-LD_GC_(OncoScan_CNV).CEL")
# pathToGCCelFile = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/data/working_data/2-AD/2-AD_GC_(OncoScan_CNV).CEL"
print(c("pathToATCelFile: ", pathToATCelFile))
print(c("pathToGCCelFile: ", pathToGCCelFile))
outputFolder = "2-AD"

custom_OS.Process(ATChannelCel = pathToATCelFile, GCChannelCel = pathToGCCelFile, samplename = outputFolder)


library(affxparser, aroma.light, ASCAT, Biostrings, BSgenome,
        changepoint, data.table, doParallel, dplyr, DT, facets,
        foreach, GenomeInfoDb, GenomicRanges, iotools, IRanges, limma,
        mclust, R.utils, rhdf5, seqinr, sequenza)
