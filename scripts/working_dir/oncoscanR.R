"
2022/03/09
Elie Bordron
OncoscanR.R

This script loads OSCHP-derived segmentation text files, processes them individually through oncoscanR, calculates Genomic Index from these results, then writes the results as a text file.
dataDir is the directory contaning all segments.txt files.
resultsDir is the output directory.
genders.txt is a table of two columns: sample(e.g. '2-AD') and gender ('M' or 'F')
"


## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
# working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/OncoscanR"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)

## importing libraries
library(oncoscanR)
library(dplyr)

## defining constants
dataDir = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/premiers_E_recus/all_segmentFiles_from_ChAS" 
resultsDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/OncoscanR"
gendersTable = read.table("C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/genders.tsv", h=T)

# loading custom Oncoscan Run functions 
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/oncoscanR_functions.R")
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/oncoscanR_scores.R")
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/onsoscanR_utils.R")

sampleNames = c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC")
# sampleNames = c("12-BC")
## running oncoscanR
main = function(dataDir, sampleNames, gendersTable) {
    source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/OncoscanR_functions.R")
    ######################################################################################
    resList = list()
    for(i in 1:length(sampleNames))  {
        currSample = sampleNames[i]
        segments_txtPattern = paste0(".*",currSample,".*","\\.segments.txt")
        nameSeg_txt = list.files(dataDir, pattern=segments_txtPattern)
        message(paste0("Processing ",nameSeg_txt))
        filepath = file.path(dataDir,nameSeg_txt)
        gender_row = gendersTable %>% dplyr::filter(sample==currSample)
        curr_gender = gender_row$gender
        before = Sys.time()
        curr_res = computeOneFileWithOncoscanR(filepath,curr_gender)
        after = Sys.time(); runTime = round(difftime(after, before, unit="secs"), 2)
        curr_res$runTime = runTime
        resList = append(resList, list(curr_res)) # caution: using `append(resList, curr_res)` would concatenate instead of append
    }
    ######################################################################################
    source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/OncoscanR_functions.R")
    GIResultsAllMethods = oncoscanR_GIs_to_dataframe(resList, sampleNames)
    ############## save GI data to file
    source(file.path(working_dir, "crossPackagesFunctions.R"))
    saveGI_ResToFile(GIResultsAllMethods, "oncoscanR")
    
    # GI_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/GI_all_methods"
    # write.table(GIResultsAllMethods,file.path(GI_dir, "gi_oncoscanR.txt"),sep="\t",row.names=FALSE, quote=F)
    
   return(GIResultsAllMethods)
}

   

res13VT = resList[[4]]
ALA13VT = res13VT[["armlevel"]]
nbAlter = 0
for (i in ALA13VT) {
    nbAlter = nbAlter + length(i)
}
nbAlter

GItable = main(dataDir, gendersTable)
GItable















