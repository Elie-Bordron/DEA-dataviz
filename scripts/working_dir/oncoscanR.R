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
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/OncoscanR_functions.R")
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/oncoscanR_scores.R")
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/onsoscanR_utils.R")

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
dataDir = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_segmentFiles_from_ChAS"
resultsDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/OncoscanR"
gendersTable = read.table("C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/genders.tsv", h=T)

### viewing data ------------------------------------------#
# paths to segments files
s5_LD_segs_path = file.path(dataDir, "5-LD.OSCHP.segments.txt")
s6_VJ_segs_path = file.path(dataDir, "6-VJ.OSCHP.segments.txt")
s8_MM_segs_path = file.path(dataDir, "8-MM.OSCHP.segments.txt")
# loading files as tables
s5_LD_segments = read.table(s5_LD_segs_path, sep='\t', h=T)[c("Type", "CN.State", "Full.Location")]
s6_VJ_segments = read.table(s6_VJ_segs_path, sep='\t', h=T)[c("Type", "CN.State", "Full.Location")]
s8_MM_segments = read.table(s8_MM_segs_path, sep='\t', h=T)[c("Type", "CN.State", "Full.Location")]
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


######## functions for exporting results to text file
getSampleNameFromOncoscanR_result = function(currRes) {
    filename = currRes$file
    splittedFilename = stringr::str_split(filename, '\\.') # need to escape the period (the dot) because this is a regex
    # print(c("splittedFilename[1]: ", splittedFilename[[1]][1]))
    sampleName = splittedFilename[[1]][1]
    return(sampleName)
}

oncoscanR_GIs_to_dataframe = function(resultsOncoscanR) {
    GIsList = list()
    samplesList = list()
    for (i in 1:length(resultsOncoscanR)) {
        currRes = resultsOncoscanR[[i]]
        ### get GI
        armlevelCN = currRes$armlevel
        curr_GI = calcGI_oncoscanR(armlevelCN)
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
    GIResultsAllMethods = oncoscanR_GIs_to_dataframe(resList)
    colnames(GIResultsAllMethods) = c("sample", "GI_oncoscanR")
    write.table(GIResultsAllMethods,file.path(resultsDir, "gi_results_all_methods.txt"),sep="\t",row.names=FALSE)
    return(GIResultsAllMethods)
}



GItable = main(dataDir, gendersTable)
GItable

















################################ oncoscanR::workflow_oncoscan.run(filepath, gender) :

#### functions for plotting seg data
plotSeg = function(seg_df, lengthOfChrs){
    ## get nb cols
    nbcols = length(seg_df)
    ## get endpos of last segment of previous chromosome
    if(grepl("X", seg_df["seqnames"]) || grepl("Y", seg_df["seqnames"]) || grepl("24", seg_df["seqnames"])){
        pos0CurrChr = sum(lengthOfChrs[1:22])
    }
    else {
        currChrArm = seg_df["seqnames"] 
        currChr = substring(currChrArm, 1, nchar(currChrArm)-1) # removing "p" or "q" char at the end of arm name, because we want the chromosome number
        if (currChr!=1){
            pos0CurrChr = sum(lengthOfChrs[1:currChr-1])
        } else {
            pos0CurrChr=0
        }
    }
    # print(c("pos0CurrChr: ", pos0CurrChr))
    ## drawing a segment on plot for each segment of the genome, using estimated values
    segStartPos = as.numeric(seg_df[["start"]])
    segEndPos = as.numeric(seg_df[["end"]])
    estimCN = as.numeric(seg_df[["cn"]])
    # print(c("pos0CurrChr + segStartPos: ", pos0CurrChr+segStartPos))
    segments(pos0CurrChr+segStartPos, estimCN, pos0CurrChr+segEndPos, estimCN, col="black", lwd=2)
}

generateGrid = function(graph_title, addAblines=F) {
    #create empty plot to add things in
    plot(1, ylim=c(0,3), xlim=c(0,3*10**9),col="white", xaxt="n", yaxt="n", xlab="nombre de bases", ylab="nombre de copies", main=graph_title)
    if (addAblines) {
        #  add horizontal grid
        for (CN in c(-7:8)) {
            abline(h=CN, col="dark grey")
        }
    }
    # X-axis
    axis(1, at = c(0, 5*10**8, 1*10**9, 1.5*10**9, 2*10**9, 2.5*10**9, 3*10**9))
    # Y-axis
    axis(2, at = c(-7:8))
}




# sampleName = "5-LD"
sampleName = "6-VJ"
# sampleName = "8-MM"
chas.fn = paste0("C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_segmentFiles_from_ChAS/",sampleName,".OSCHP.segments.txt")
# chas.fn = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_segmentFiles_from_ChAS/6-VJ_from_analysis_setup.segment.txt"
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

if(F){
    ### plotting before/after each filter/smooth step
    plotBeforeAfterCleanStep = function(segs_df, cleaned_df, action){
        lengthOfChrs = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
        generateGrid(paste0(sampleName," before ", action))
        apply(segs_df, 1, plotSeg, lengthOfChrs)
        generateGrid(paste0(sampleName," after ", action))
        apply(cleaned_df, 1, plotSeg, lengthOfChrs)
    }
    ## trimming
    trimmed = as.data.frame(trim_to_coverage(segments, oncoscan.cov))
    plotBeforeAfterCleanStep(segsDf, trimmed, "trimming segments to coverage")
    ## adjusting LOH
    lohAdjusted = as.data.frame(adjust_loh(segments))
    plotBeforeAfterCleanStep(segsDf, lohAdjusted, "adjusting LOH")
    ## merging segments
    merged = as.data.frame(merge_segments(segments))
    plotBeforeAfterCleanStep(segsDf, merged, "merging segments")
    ## removing segments
    pruned = as.data.frame(prune_by_size(segments))
    plotBeforeAfterCleanStep(segsDf, pruned, "removing segments")
}

segs.clean <- trim_to_coverage(segments, oncoscan.cov) %>% # All segments that are not entirely contained within the kit coverage will be trimmed to the coverage's limits.
    adjust_loh() %>% # LOH segments completely contained within (or equal to) a copy loss segment are deleted.   LOH segments partially overlapping (on one end only) with a copy loss segment are trimmed to remove the overlap. If a copy loss segment is completely contained within (but not equal to) a LOH segment, then nothing is done; the overlap remains.
    merge_segments() %>% # merge segments that are at a distance smaller than the resolution (300Kb)  (only occurs if the segments have the same copy number).
    prune_by_size() # remove segments smaller than 300kb. this value = Oncoscan assay resolution outside of cancer genes
## all steps
if (F){
    segsCleanDf = as.data.frame(segs.clean) # for 6VJ sample, in this step, the two 22q segs are fused together.
    row.names(segsCleanDf) = segsCleanDf$seqnames
    plotBeforeAfterCleanStep(segsDf, segsCleanDf, "4-steps cleaning")
}


#### get ALA ####
## split segs in 4 objects, one per type. ##
# Split segments by type: Loss, LOH, gain or amplification and get the arm-level alterations.
# Note that the segments with copy gains include all amplified segments.
custom_threshold = 0.1
armlevel.loss <- segs.clean[segs.clean$cn.type == cntype.loss] %>%
    armlevel_alt(kit.coverage = oncoscan.cov, custom_threshold)
armlevel.loh <- segs.clean[segs.clean$cn.type == cntype.loh] %>%
    armlevel_alt(kit.coverage = oncoscan.cov, custom_threshold)
armlevel.gain <- segs.clean[segs.clean$cn.type == cntype.gain] %>%
    armlevel_alt(kit.coverage = oncoscan.cov, custom_threshold)
armlevel.amp <- segs.clean[segs.clean$cn.subtype %in% c(cntype.strongamp, cntype.weakamp)] %>%
    armlevel_alt(kit.coverage = oncoscan.cov, custom_threshold)

### plotting ALA
OScovDf = as.data.frame(oncoscan.cov)
getPAAs = function(oncoscan_cov,armlevelAlter){
    currChrArm = oncoscan_cov["seqnames"]
    if (currChrArm %in% names(armlevelAlter)) {
        return(armlevelAlter[currChrArm])
    }
    return(0)
}
plotPAAWithThreshold = function(armLevelAlt_1type, plotTitle){
    OScovDf$PAA_1type = apply(OScovDf, 1, getPAAs, armLevelAlt_1type)
    OScovDf = mutate(OScovDf, isTrueAlter=ifelse(OScovDf$PAA_1type > 0.8, "#b75869", "#5869b7")) # red, blue
    barplot(OScovDf$PAA_1type, names.arg = OScovDf$seqnames, ylim=c(0,1), main=plotTitle, col = "#5869b7")
    barplot(OScovDf$PAA_1type, names.arg = OScovDf$seqnames, ylim=c(0,1), main=plotTitle, col = "#5869b7")
    abline(h=0.8)
    barplot(OScovDf$PAA_1type, names.arg = OScovDf$seqnames, ylim=c(0,1), main=plotTitle, col = OScovDf$isTrueAlter)
    abline(h=0.8)
}

## displaying chromosomes arms in the normal order
row.names(OScovDf) = OScovDf$seqnames
rowsLayout = c("1p","1q", "2p","2q", "3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q","Xp", "Xq","Yp", "Yq")
OScovDf = OScovDf[rowsLayout,]
# creating color group
# save_plots_here = file.path(resultsDir, sampleName)
# if(!dir.exists(save_plots_here)) dir.create(save_plots_here)
# barplot(OScovDf$PAA_loss, names.arg = OScovDf$seqnames, ylim=c(0,1), main=paste0(sampleName, " loss"), col = "#5869b7")
# abline(h=0.8)
if (F) {
    plotPAAWithThreshold(armlevel.loss, paste0(sampleName, " loss"))
    plotPAAWithThreshold(armlevel.loh, paste0(sampleName, " loh"))
    plotPAAWithThreshold(armlevel.gain, paste0(sampleName, " gain"))
    plotPAAWithThreshold(armlevel.amp, paste0(sampleName, " amp"))
}

## armlevel.loss contains all chromosome arms characterized with a loss,and the % of coverage of each chromosome arm by loss-altered segments.
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

## known. see  logbook.md or cahier. number of LST

n.td <- score_td(segs.clean)
## known. see  logbook.md or cahier. tdplus score .

# Get the alterations into a single list and print it in a JSON format.
armlevel_alt.list <- list(AMP=sort(names(armlevel.amp)),
                          LOSS=sort(names(armlevel.loss)),
                          LOH=sort(names(armlevel.loh)),
                          GAIN=sort(names(armlevel.gain)))

scores.list <- list(HRD=paste0(hrd['HRD'], ', nLST=', hrd['nLST']), TDplus=n.td$TDplus,
                    avgCN=substr(as.character(wgd['avgCN']), 1, 4))


obj_to_return = list(armlevel=armlevel_alt.list,
            scores=scores.list,
            gender=gender,
            file=basename(chas.fn))
