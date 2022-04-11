library(EaCoN)
library("ASCAT")
packagesImportedByEaCoN = c("affxparser", "aroma.light", "ASCAT", "Biostrings", "BSgenome",
                            "changepoint", "data.table", "doParallel", "dplyr", "DT",
                            "foreach", "GenomeInfoDb", "GenomicRanges", "iotools", "IRanges", "limma",
                            "mclust", "R.utils", "rhdf5", "seqinr", "sequenza")
lapply(packagesImportedByEaCoN, require, character.only = TRUE)

# to load functions from EaCoN_functions.R

## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)


pathToOSCHP = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_OSCHP/5-LD.OSCHP"
pathToCel = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_CEL"
pathToATCelFile = file.path(pathToCel, "5-LD_AT_(OncoScan_CNV).CEL")
pathToGCCelFile = file.path(pathToCel, "5-LD_GC_(OncoScan_CNV).CEL")
# print(c("pathToATCelFile: ", pathToATCelFile))
# print(c("pathToGCCelFile: ", pathToGCCelFile))
outputFolder = "5-LD"

source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/EaCoN_functions.R")
ASCAT_obj = custom_OS.Process(ATChannelCel = pathToATCelFile, GCChannelCel = pathToGCCelFile, samplename = outputFolder, oschp_file=pathToOSCHP, force=T, oschp.keep=T, return.data=TRUE) 
#arguments:
    #force=T is to discard results files before recreating them.
    #oschp.keep=T is to not discard oschp file after processing.

### to check that readRDS() returns the ASCAT object. this means we can do the pipeline without saving to .RDS if needed.
# x = data.frame(c(9,8,5,6,3),c("a","a","z","e","a"))
# saveRDS(x, file.path(working_dir,"testDf.RDS"))
# yvaldftest = readRDS(file.path(working_dir,"testDf.RDS"))
# x!=yvaldftest
RDSInitPath = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/5-LD_OncoScan_CNV_hg19_processed.RDS"
Segment.ff(RDS.file = RDS_path, segmenter = "ASCAT", force=T)
RDSSegPath = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/ASCAT/L2R/5-LD.SEG.ASCAT.RDS"
ASCN.ff(RDS.file = RDSSegPath)
RDSCallPath = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/ASCAT/ASCN/gamma0.60/5-LD.ASCN.ASCAT.RDS"
Annotate.ff(RDS.file = RDSCallPath, author.name = "E.Bordron")

