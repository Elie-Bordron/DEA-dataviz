library(EaCoN)
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


loaded_OSCHP <- oschp.load(file = pathToOSCHP)
class(loaded_OSCHP)

### this works
# x = data.frame(c(98,55,68,5,21),c("p", "a","p","p","a"))
# saveRDS(x,"C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/5-LD/a.RDS", compress = "bzip2")
