"
2022/03/09
Elie Bordron
OncoscanR.R

This script loads OSCHP-derived segmentation text files, processes them individually through oncoscanR, calculates Genomic Index from these results, then writes the results as a text file.
dataDir is the directory contaning all segments.txt files.
resultsDir is the output directory.
genders.txt is a table of two columns: sample(e.G. '2-AD') and gender ('M' or 'F')
"


## importing libraries
library(oncoscanR)
library(dplyr)


## defining constants
setwd( "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results")
dataDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/data/working_data/from_CEL"
resultsDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results"
gendersTable = read.table("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/data/working_data/genders.txt", h=T)


######## functions for iterating over files

computeOneFileWithOncoscanR = function(curr_file,dataDir,gendersTable) {
    fileFullPath = file.path(dataDir, curr_file)
    print(c("fileFullPath: ", fileFullPath))
    gender_row = gendersTable %>% dplyr::filter(stringr::str_detect(curr_file, sample))
    curr_gender = gender_row$gender
    print(c("curr_gender: ", curr_gender))
    curr_res = oncoscanR::workflow_oncoscan.run(fileFullPath, curr_gender)
    return(curr_res)
}


computeAllFilesWithOncoscanR = function(dataDir, gendersTable) {
    resList = list()
    files <- list.files(dataDir)
    for(i in 1:length(files))  {
        curr_file = files[i]
        isSegmentsFile = grepl("segments.txt", curr_file)
        if(isSegmentsFile) {
            curr_res = computeOneFileWithOncoscanR(curr_file,dataDir,gendersTable)
            resList = append(resList, list(curr_res)) # caution: using `append(resList, curr_res)` would concatenate instead of append
        }
    }
    return(resList)
}



######## functions for calculating GI

getChrIdsFromAlterVec = function(alterVec) {
    chrIds = c()
    if (length(alterVec)!=0) {
        for(i in 1:length(alterVec)) {
            currentAlter = alterVec[i]
            nbChar = nchar(currentAlter)
            chr_id = substr(currentAlter, 0, nbChar-1) 
            chrIds = append(chrIds, chr_id)
        }
    }
    return(chrIds)
}

getChrIds_oncoscanR = function(armLevelAlter) {
    losses = armLevelAlter$LOSS
    gains = armLevelAlter$GAIN
    chrIds_loss = getChrIdsFromAlterVec(losses)
    chrIds_gain = getChrIdsFromAlterVec(gains)
    chrIds = c(chrIds_loss, chrIds_gain)
    return(chrIds)
}

getNbChrs_oncoscanR = function(chrIds) {
    alteredChrs = unique(chrIds)
    nbChrs = length(alteredChrs)
    return(nbChrs)
}

calcGI = function(nbAlter, nbChrs) {
    if (nbAlter==nbChrs) {
        GI=nbAlter ## to avoid dividing by zero
    }else {
        GI=(nbAlter**2)/nbChrs
    }
    return(GI)
        
}

getGIParams_oncoscanR = function(armLevelAlter) {
    losses = armLevelAlter$LOSS
    gains = armLevelAlter$GAIN
    nbAlter = length(losses) + length(gains)
    chrIds = getChrIds_oncoscanR(armLevelAlter)
    nbChrs = getNbChrs_oncoscanR(chrIds)
    return(list(nbAlter, nbChrs))
}

calcGI_oncoscanR = function(armLevelAlter) {
    params = getGIParams_oncoscanR(armLevelAlter)
    nbAlter = params[[1]]
    nbChrs = params[[2]]
    GI = calcGI(nbAlter, nbChrs)
    return(GI)    
}

getGIFromOncoscanR_result = function(currRes) {
    armlevelCN = currRes$armlevel
    GI = calcGI_oncoscanR(armlevelCN)
    return(GI)
}

######## functions for exporting results to text file
getSampleNameFromOncoscanR_result = function(currRes) {
    filename = currRes$file
    splittedFilename = stringr::str_split(filename, '\\.') # need to escape the period because this is a regex
    print(c("splittedFilename[1]: ", splittedFilename[[1]][1]))
    sampleName = splittedFilename[[1]][1]
    return(sampleName)
}

oncoscanR_results_to_dataframe = function(resultsOncoscanR) {
    GIsList = list()
    samplesList = list()
    for (i in 1:length(resultsOncoscanR)) {
        currRes = resultsOncoscanR[[i]]
        ### get GI
        curr_GI = getGIFromOncoscanR_result(currRes)
        GIsList = append(GIsList, curr_GI)
        ### get sample name
        currSampleName = getSampleNameFromOncoscanR_result(currRes)
        samplesList = append(samplesList, currSampleName)
    }
    GIResultsAllMethods = cbind(samplesList, GIsList)
    return(GIResultsAllMethods)
}


main = function(dataDir, gendersTable) {
    resList = computeAllFilesWithOncoscanR(dataDir, gendersTable)
    GIResultsAllMethods = oncoscanR_results_to_dataframe(resList)
    colnames(GIResultsAllMethods) = c("sample", "GI_oncoscanR")
    write.table(GIResultsAllMethods,file.path(resultsDir, "gi_results_all_methods.txt"),sep="\t",row.names=FALSE)
    return(GIResultsAllMethods)
}



val = main(dataDir, gendersTable)
plot()