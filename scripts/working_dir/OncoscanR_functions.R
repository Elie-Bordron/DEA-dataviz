####### 26/04/2022 this code was obtained from https://rdrr.io/github/yannchristinat/oncoscanR/src/R/workflow_oncoscan.R
# workflow_oncoscan.R
#
# Functions to run the complete workflow from input files to scores and arm-level alterations.
#
# Author: Yann Christinat
# Date: 23.06.2020

#' Run the standard workflow for Oncoscan ChAS files.
#'
#' @details Identifies the globally altered arms (\>=80\% of arm altered), computes the HRD and
#' TD+ scores. The amplification is defined as a CN subtype \code{cntype.weakamp} or
#' \code{cntype.strongamp}. An arm is gained if of CN type \code{cntype.gain} unless the arm is
#' amplified.
#'
#' @param chas.fn Path to the text-export ChAS file
#' @param gender Gender of the sample (M or F)
#'
#' @return A list of lists with the following elements:
#' \code{armlevel = list(AMP= list of arms, GAIN= list of arms, LOSS= list of arms, LOH= list of arms),
#' scores = list(LST= number, LOH= number, TDplus= number, TD= number),
#' gender = gender as given by the parameter,
#' file = path of the ChAS file as given by the parameter)}
#'
#' @export
#'
#' @import magrittr
#'
#' @examples
#' segs.filename <- system.file("extdata", "chas_example.txt", package = "oncoscanR")
#' workflow_oncoscan.run(segs.filename, "M")
custom_workflow_oncoscan.run <- function(chas.fn, gender){
    if(!(gender %in% c('M', 'F'))){
        stop("The gender (second argument) has to be F or M.")
    }

    # Remove the 21p arm from the Oncoscan coverage as it is only partly covered and we don't
    # want to return results on this arm.
    oncoscan.cov <- oncoscanR::oncoscan_na33.cov[seqnames(oncoscanR::oncoscan_na33.cov) != '21p']

    # Load the ChAS file and assign subtypes.
    segments <- load_chas(chas.fn, oncoscan.cov)
    segments$cn.subtype <- get_cn_subtype(segments, gender)

    # Clean the segments: resctricted to Oncoscan coverage, LOH not overlapping with copy loss
    # segments, smooth&merge segments within 300kb and prune segments smaller than 300kb.
    segs.clean <- trim_to_coverage(segments, oncoscan.cov) %>%
    adjust_loh() %>%
    merge_segments() %>%
    prune_by_size()

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

    # Remove amplified segments from armlevel.gain
    armlevel.gain <- armlevel.gain[!(names(armlevel.gain) %in% names(armlevel.amp))]

    # Get the number of nLST and TDplus
    wgd <- score_estwgd(segs.clean, oncoscanR::oncoscan_na33.cov) # Get the avg CN, including 21p
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
}





score_estwgd <- function(segments, kit.coverage){
  # Get the average copy number
  avgcn <- score_avgcn(segments, kit.coverage)

  # Estimate the number of WGD events
  wgd.est <- ifelse(avgcn > 3.4, 2, ifelse(avgcn > 2.2, 1, 0))

  return(c(WGD=wgd.est, avgCN=avgcn))
}










############################################# functions for automating GI calculation

######## functions for iterating over files

computeOneFileWithOncoscanR = function(filepath,gender) {
    print(c("filepath: ", filepath))
    print(c("gender: ", gender))
    #### using custom run fx
    custom_workflow_oncoscan.run(filepath, gender)
    curr_res = oncoscanR::workflow_oncoscan.run(filepath, gender)
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
