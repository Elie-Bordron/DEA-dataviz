"
2022/03/09
Elie Bordron
OncoscanR.R

This script loads OSCHP-derived segmentation text files, processes them individually through oncoscanR, calculates Genomic Index from these results, then writes the results as a text file.
dataDir is the directory contaning all segments.txt files.
resultsDir is the output directory.
genders.txt is a table of two columns: sample(e.g. '2-AD') and gender ('M' or 'F')
"



# loading custom Oncoscan Run function 
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/OncoscanR_functions.R")


## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)

## importing libraries
library(oncoscanR)
library(dplyr)


## defining constants
dataDir = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_segmentFiles_from_ChAS"
resultsDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/OncoscanR"
gendersTable = read.table("C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/genders.tsv", h=T)

### viewing data ------------------------------------------#
# paths to probeset.txt files
s5_LD_path = file.path(dataDir, "5-LD.OSCHP.segments.txt")
s6_VJ_path = file.path(dataDir, "6-VJ.OSCHP.segments.txt")
s8_MM_path = file.path(dataDir, "8-MM.OSCHP.segments.txt")
# loading files as tables
s5_LD_segments = read.table(s5_LD_path, sep='\t', h=T)
s6_VJ_segments = read.table(s6_VJ_path, sep='\t', h=T)
s8_MM_segments = read.table(s8_MM_path, sep='\t', h=T)
### -------------------------------------------------------#


######## functions for iterating over files

computeOneFileWithOncoscanR = function(filepath,gender) {
    print(c("filepath: ", filepath))
    print(c("gender: ", gender))
    #### using custom run fx
    custom_workflow_oncoscan.run(filepath, gender)
    # curr_res = oncoscanR::workflow_oncoscan.run(filepath, gender)
    return(curr_res)
}


computeAllFilesWithOncoscanR = function(dataDir, gendersTable) {
    resList = list()
    files <- list.files(dataDir)
    for(i in 1:length(files))  {
        curr_file = files[i]
        isSegmentsFile = grepl("segments.txt", curr_file)
        # print(c("files: ", files))
        # print(c("isSegmentsFile: ", isSegmentsFile))
        if(isSegmentsFile) {
            filepath = file.path(dataDir, curr_file)
            # print(c("filepath: ", filepath))
            # print(c("curr_file: ", curr_file))
            # print(c("sample: ", sample))
            gender_row = gendersTable %>% dplyr::filter(stringr::str_detect(curr_file, sample))
            curr_gender = gender_row$gender
            curr_res = computeOneFileWithOncoscanR(filepath,curr_gender)
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
    splittedFilename = stringr::str_split(filename, '\\.') # need to escape the period (the dot) because this is a regex
    # print(c("splittedFilename[1]: ", splittedFilename[[1]][1]))
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



GItable = main(dataDir, gendersTable)
GItable








################################ oncoscanR::workflow_oncoscan.run(filepath, gender) :
chas.fn = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_segmentFiles_from_ChAS/6-VJ.OSCHP.segments.txt"
gender = "M"

# Remove the 21p arm from the Oncoscan coverage as it is only partly covered and we don't
# want to return results on this arm.
oncoscan.cov <- oncoscanR::oncoscan_na33.cov[seqnames(oncoscanR::oncoscan_na33.cov) != '21p']

#### load ####
# Load the ChAS file and assign subtypes.
segments <- load_chas(chas.fn, oncoscan.cov)
segments$cn.subtype <- get_cn_subtype(segments, gender)
segsDf = as.data.frame(segments)

#### clean ####
# Clean the segments: resctricted to Oncoscan coverage, LOH not overlapping with copy loss
# segments, smooth&merge segments within 300kb and prune segments smaller than 300kb.
segs.clean <- trim_to_coverage(segments, oncoscan.cov) %>% # All segments that are not entirely contained within the kit coverage will be trimmed to the coverage's limits.
    adjust_loh() %>% # LOH segments completely contained within (or equal to) a copy loss segment are deleted. LOH segments partially overlapping (on one end only) with a copy loss segment are trimmed to remove the overlap. If a copy loss segment is completely contained within (but not equal to) a LOH segment, then nothing is done; the overlap remains.
    merge_segments() %>% # merge segments that are at a distance smaller than the resolution (300Kb)  (only occurs if the segments have the same copy number).
    prune_by_size() # remove segments smaller than 300kb. this value = Oncoscan assay resolution outsode of cancer genes
segsCleanDf = as.data.frame(segs.clean) # for 6VJ sample, in this step, the two 22q segs are fused together.



#### get ALA ####
## split segs in 4 objects, one per type. ##
# Split segments by type: Loss, LOH, gain or amplification and get the arm-level alterations.
# Note that the segments with copy gains include all amplified segments.
armlevel.loss <- segs.clean[segs.clean$cn.type == cntype.loss] %>%
    armlevel_alt(kit.coverage = oncoscan.cov)
armlevel.loh <- segs.clean[segs.clean$cn.type == cntype.loh] %>%
    armlevel_alt(kit.coverage = oncoscan.cov)
armlevel.gain <- segs.clean[segs.clean$cn.type == cntype.gain] %>%
    armlevel_alt(kit.coverage = oncoscan.cov)
armlevel.amp <- segs.clean[segs.clean$cn.subtype %in% c(cntype.strongamp, cntype.weakamp)] %>%
    armlevel_alt(kit.coverage = oncoscan.cov)
## armlevel.loss contains all chromosome arms characterized with a loss, and the % of coverage of each chromosome arm by loss-altered segments.
## same for the others respectively.


# Remove amplified segments from armlevel.gain
armlevel.gain <- armlevel.gain[!(names(armlevel.gain) %in% names(armlevel.amp))]


#### get scores ####
# Get the number of nLST and TDplus
### wgd: the NUMBER of WGD events.
wgd <- score_estwgd(segs.clean, oncoscanR::oncoscan_na33.cov) # Get the avg CN, including 21p
# if average copy number is > 3.4, wgd=2.
# else, if average copy number is > 2.2, wgd=1.
# else, wgd=0.

hrd <- score_nlst(segs.clean, wgd['WGD'], oncoscan.cov)

n.td <- score_td(segs.clean)

# Get the alterations into a single list and print it in a JSON format.
armlevel_alt.list <- list(AMP=sort(names(armlevel.amp)),
                          LOSS=sort(names(armlevel.loss)),
                          LOH=sort(names(armlevel.loh)),
                          GAIN=sort(names(armlevel.gain)))
scores.list <- list(HRD=paste0(hrd['HRD'], ', nLST=', hrd['nLST']), TDplus=n.td$TDplus,
                    avgCN=substr(as.character(wgd['avgCN']), 1, 4))

return(list(armlevel=armlevel_alt.list,
            scores=scores.list,
            gender=gender,
            file=basename(chas.fn)))