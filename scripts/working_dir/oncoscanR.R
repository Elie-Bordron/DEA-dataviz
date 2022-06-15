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
# sampleNames = c("11-BG")
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

##### build segTables
lengthOfChrs = c(0, 247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
posOfCentromeres = c(123.4, 93.9, 90.9, 50.0, 48.8, 59.8, 60.1, 45.2, 43.0, 39.8, 53.4, 35.5, 17.7, 17.2, 19.0, 36.8, 25.1, 18.5, 26.2, 28.1, 12.0, 15.0, 61.0, 10.4) * 10**6

CN_per_alter = list("GAIN"=3, "LOSS"=1)

segTables = list()
for (s in 1:length(sampleNames)) {
    print(c("treating sample", sampleNames[s]))
    ## initialize by-segments table 
    segTable = data.frame(matrix(ncol = 5, nrow = 0))
    colnames(segTable) = c("sample", "Chromosome", "Start", "End", "CN")
    segId = 1
    currSampleRes = resList[[s]]
    ALA = currSampleRes$armlevel
    for (currAlter in c("GAIN", "LOSS")) {
        if(length(ALA[currAlter][[1]])!=0) {
            for (arm in ALA[currAlter][[1]]) {
                currRow = list()
                ## find if arm is q or p
                print(c("arm: ", arm))
                if(grepl("p", arm)) {
                    armID = "p"
                } else if(grepl("q", arm)) {
                    armID = "q"
                } else {print("arm name doesn't contain p nor q")}
                splitted = stringr::str_split(arm, armID)
                chr = splitted[[1]][1]
                currRow$sample = sampleNames[s]
                currRow$chr = chr
                ## get arm start and stop positions
                chrNum = as.numeric(chr)
                if(armID=="p") {
                    currRow$start = sum(lengthOfChrs[1:chrNum])
                    currRow$end = sum(lengthOfChrs[1:chrNum]) + posOfCentromeres[chrNum]
                } else if(armID=="q") {
                    currRow$start = sum(lengthOfChrs[1:chrNum]) + (lengthOfChrs[chrNum+1] - posOfCentromeres[chrNum])
                    currRow$end =  sum(lengthOfChrs[1:chrNum+1])
                } else {print("armID is not p nor q")}
                currRow$CN = CN_per_alter[currAlter][[1]]
                print(c("currRow$end - currRow$start: ", currRow$end - currRow$start) )
                print(c("currRow: ", currRow))
                segTable[segId,] = currRow
                segId = segId+1
            }
        }
        
    }
    # segTables[sampleNames[s]] = list(segTable)
    segTables = append(segTables, list(segTable))
}




############## save segments table to file
saveSegTables = function(segTables, outputDir, sampleNames=c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )) {
    segTablesDir = file.path(outputDir, "segTables")
    if(!dir.exists(segTablesDir))dir.create(segTablesDir)
    for (i in 1:length(segTables)) {
        currSegTable = segTables[i]
        currFilePath = file.path(segTablesDir, paste0(sampleNames[i], ".tsv"))
        write.table(currSegTable, currFilePath, sep="\t", row.names=FALSE, quote=F)
    }
}
saveSegTables(segTables, resultsDir)




GItable = main(dataDir, gendersTable)
GItable















