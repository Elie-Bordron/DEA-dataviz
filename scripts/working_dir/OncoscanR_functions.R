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
custom_workflow_oncoscan.run <- function(chas.fn, gender, cleanSegments=TRUE, plotPAA=F){
    if(!(gender %in% c('M', 'F'))){
        stop("The gender has to be 'F' or 'M'.")
    }

    # Remove the 21p arm from the Oncoscan coverage as it is only partly covered and we don't
    # want to return results on this arm.
    oncoscan.cov <- oncoscanR::oncoscan_na33.cov[seqnames(oncoscanR::oncoscan_na33.cov) != '21p']

    # Load the ChAS file and assign subtypes.
    
    
    if (FALSE) { ## to make internship defense plotSeg
        ############ pour slides soutenance
        soutenance_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/res_soutenance"
        # ylim = c(-3,3)
        ylim = c(1,4)
        PercentRowsToRemove = 30
        xlab = "position genomique"
        pkgName="CGHcall"
        ###### plot dimensions:
        pngWidth = 500
        pngHeight = 430
        sampleName = params$sampleNames[1]


        chas.fn = filepath
        gender = curr_gender
        sampleName = "11-BG"
        working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
        source(file.path(working_dir, "CGHcall_functions.R"))
        source(file.path(working_dir, "oncoscanR_functions.R"))
        O_segs = read.table(chas.fn, sep='\t', h=TRUE)
        
        ### get start and end positions
        for (i in 1:length(O_segs$Full.Location)) {
            segLoc = O_segs$Full.Location[i]
            val = stringr::str_split(segLoc, ":")[[1]][2]
            val = stringr::str_split(val, "-")[[1]]
            O_segs[i, "loc.start"] = val[1]
            O_segs[i, "loc.end"] = val[2]
        }

        O_segsForPlot = O_segs[c("Chromosome", "loc.start", "loc.end", "CN.State", "Marker.Count")] ## Marker.Count is for replacing "nbProbes column".
        pathToSaveFile = paste0(working_dir, "/plots/oncosR_", sampleName, "_raw.png")
        png(pathToSaveFile, width = pngWidth, height = pngHeight)
        plotSegTable(O_segsForPlot, sampleName, savePlot=FALSE, ylim = ylim, mainTitlePrefix = " raw ")
        dev.off()
    }
    
    segments <- load_chas(chas.fn, oncoscan.cov)
    segments$cn.subtype <- get_cn_subtype(segments, gender)
    # print(c("segs with subtypes: ", segments$cn.subtype))

    # Clean the segments: resctricted to Oncoscan coverage, LOH not overlapping with copy loss segments, smooth&merge segments within 300kb and prune segments smaller than 300kb.
    if(!cleanSegments) {
        message("not cleaning segments.")
        segs.clean = trim_to_coverage(segments, oncoscan.cov)
        # print(c("segs after minimal trimming: ", segs.clean))
    } else{
        message("cleaning segments...")
        segs.clean <- trim_to_coverage(segments, oncoscan.cov) %>%
        adjust_loh() %>%
        merge_segments() %>%
        prune_by_size()
    }
    
    if(FALSE) {  ## to make internship defense plots
        O_segs_cleaned = as.data.frame(segs.clean)
        O_segs_cleaned_for_plot = O_segs_cleaned[c("seqnames", "start", "end", "cn", "width")] ## width is for replacing "nbProbes column".
        ## get chrs from arms
        O_segs_cleaned_for_plot$seqnames = as.character(O_segs_cleaned_for_plot$seqnames)
        for(i in 1:length(O_segs_cleaned_for_plot$seqnames)) {
            currSeg = O_segs_cleaned_for_plot$seqnames[[i]]
            noq = stringr::str_split(currSeg, "q")[[1]][1]
            nop = stringr::str_split(noq, "p")[[1]][1]
            print(nop)
            O_segs_cleaned_for_plot[i, "chrom"] = nop
        }
        O_segs_cleaned_for_plot = O_segs_cleaned_for_plot[c("chrom", "start", "end", "cn", "width")]
        
        pathToSaveFile = paste0(working_dir, "/plots/oncosR_", sampleName, "_postClean.png")
        png(pathToSaveFile, width = pngWidth, height = pngHeight)
        plotSegTable(O_segs_cleaned_for_plot, sampleName, savePlot=FALSE, ylim = c(1,4), mainTitlePrefix = " after cleaning")
        dev.off()
    }
        
    # Split segments by type: Loss, LOH, gain or amplification and get the arm-level alterations.
    # Note that the segments with copy gains include all amplified segments.
    customThreshold = 0.8
    armlevel.loss <- segs.clean[segs.clean$cn.type == cntype.loss] %>%
    armlevel_alt(kit.coverage = oncoscan.cov)
    armlevel.loh <- segs.clean[segs.clean$cn.type == cntype.loh] %>%
    armlevel_alt(kit.coverage = oncoscan.cov)
    armlevel.gain <- segs.clean[segs.clean$cn.type == cntype.gain] %>%
    armlevel_alt(kit.coverage = oncoscan.cov)
    armlevel.amp <- segs.clean[segs.clean$cn.subtype %in% c(cntype.strongamp, cntype.weakamp)] %>%
    armlevel_alt(kit.coverage = oncoscan.cov)

    ## retrieve sample Name
    fileName = stringr::str_split(chas.fn,'/')[[1]]
    lastElmt = fileName[[length(fileName)]]
    sampleName = stringr::str_split(lastElmt,'\\.')[[1]]
    sampleName = sampleName[[1]]
    if(plotPAA) {
        message("plotting proportions of altered arms...")
        plotPAAWithThreshold(armlevel.loss, paste0(sampleName, " loss"),oncoscan.cov)
        plotPAAWithThreshold(armlevel.loh, paste0(sampleName, " loh"),oncoscan.cov)
        plotPAAWithThreshold(armlevel.gain, paste0(sampleName, " gain"),oncoscan.cov)
        plotPAAWithThreshold(armlevel.amp, paste0(sampleName, " amp"),oncoscan.cov)
    }



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





# score_estwgd <- function(segments, kit.coverage){
#     # Get the average copy number
#     avgcn <- score_avgcn(segments, kit.coverage)

#     # Estimate the number of WGD events
#     wgd.est <- ifelse(avgcn > 3.4, 2, ifelse(avgcn > 2.2, 1, 0))

#     return(c(WGD=wgd.est, avgCN=avgcn))
# }










############################################# functions for automating GI calculation

######## functions for iterating over files

computeOneFileWithOncoscanR = function(filepath,gender, plotPAA=F) {
    # print(c("filepath: ", filepath))
    # print(c("gender: ", gender))
    #### using custom run fx
    curr_res = custom_workflow_oncoscan.run(filepath, gender,plotPAA)
    # curr_res = oncoscanR::workflow_oncoscan.run(filepath, gender)
    # print(c("computed result: ", curr_res))
    return(curr_res)
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
    chrIds = c()
    for (AlterType in armLevelAlter) {
        chrIds = c(chrIds,getChrIdsFromAlterVec(AlterType))
    }
    # losses = armLevelAlter$LOSS
    # gains = armLevelAlter$GAIN
    # amps = armLevelAlter$AMP
    # lohs = armLevelAlter$LOH
    # chrIds_loss = getChrIdsFromAlterVec(losses)
    # chrIds_gain = getChrIdsFromAlterVec(gains)
    # chrIds_amp = getChrIdsFromAlterVec(amps)
    # chrIds_loh = getChrIdsFromAlterVec(lohs)
    # chrIds = c(chrIds_loss, chrIds_gain, chrIds_amp, chrIds_loh)
    return(chrIds)
}

getNbChrs_oncoscanR = function(chrIds) {
    alteredChrs = unique(chrIds)
    nbChrs = length(alteredChrs)
    return(nbChrs)
}

calcGI = function(nbAlter, nbChrs) {
    if(nbChrs > nbAlter) {
        warning("nbChrs > nbAlter. Did you invert them?")
        return(NULL)
    if(nbChrs>24) {
        warning("nbChrs > 24. Did you invert nbChrs and nbAlter?")
        return(NULL)
    }
    if(nbChrs>22) {
        warning("nbChrs > 22. Please make sure you did not invert nbChrs and nbAlter, and that you did not include sex chromosomes data")
    }
    }else{
        if (nbAlter==nbChrs) {
            GI=nbAlter ## to avoid dividing by zero
        }else {
            print(c("nbAlter: ", nbAlter))
            print(c("nbAlter^2: ", nbAlter**2))
            print(c("nbChrs: ", nbChrs))
            GI=(nbAlter**2)/nbChrs
            print(c("GI in calcGI() : ", GI))
        }
    }
    return(GI)   
}

getGIParams_oncoscanR = function(armLevelAlter) {
    # print(c("armLevelAlter: ", armLevelAlter))
    nbAlter = 0
    # for (AlterType in armLevelAlter) {
    #     print(c("AlterType: ", AlterType))
    #     nbAlter = nbAlter + length(AlterType)
    #     print(c("nbAlter: ", nbAlter))
    # }
    losses = armLevelAlter$LOSS
    gains = armLevelAlter$GAIN
    amps = armLevelAlter$AMP
    lohs = armLevelAlter$LOH
    alterGain = c(gains, amps)
    alterLoss = c(losses, lohs)
    # print("==================================")
    alterGain = unique(alterGain)
    alterLoss = unique(alterLoss)
    # print(c("alterGain: ", alterGain))
    # print(c("alterLoss: ", alterLoss))
    chrIds = getChrIds_oncoscanR(armLevelAlter)
    # print(c("chrIds: ", chrIds))
    nbChrs = length(unique(chrIds))
    nbAlter = length(c(alterGain, alterLoss))
    return(list(nbAlter, nbChrs))
}

calcGI_oncoscanR = function(armLevelAlter) {
    # print(c("armLevelAlter passed: ", armLevelAlter))
    params = getGIParams_oncoscanR(armLevelAlter)
    nbAlter = params[[1]]
    nbChrs = params[[2]]
    GI = calcGI(nbAlter, nbChrs)
    return(list(GI,nbAlter, nbChrs))    
}


######## functions for exporting results to text file
getSampleNameFromOncoscanR_result = function(currRes) {
    filename = currRes$file
    splittedFilename = stringr::str_split(filename, '\\.') # need to escape the period (the dot) because this is a regex
    sampleName = splittedFilename[[1]][1]
    return(sampleName)
}

oncoscanR_GIs_to_dataframe = function(resultsOncoscanR, sampleNames) {
    GI_oncoscanR_df = data.frame(matrix(ncol = 4, nrow = length(sampleNames)))
    colnames(GI_oncoscanR_df) = c("GI", "nbAlter", "nbChr", "runTime")
    rownames(GI_oncoscanR_df) = sampleNames
    for (i in 1:length(resultsOncoscanR)) {
        currRes = resultsOncoscanR[[i]]
        ### get GI
        # print(c("currRes: ", currRes))
        # print(c("currRes$armlevel: ", currRes$armlevel))
        armlevelAlters = currRes$armlevel
        # print(c("armlevelAlters: ", armlevelAlters))
        curr_GI_res = calcGI_oncoscanR(armlevelAlters)
        # print(c("curr_GI_res: ", curr_GI_res))
        # print( list(sampleNames[i],curr_GI_res))
        GI_oncoscanR_df[i,] = append(curr_GI_res, currRes$runTime)
        # print(GI_oncoscanR_df)
    }
    return(GI_oncoscanR_df)
}












#### functions for plotting seg data
plotSeg = function(seg_df, lengthOfChrs){
    ## get nb cols
    nbcols = length(seg_df)
    ## get endpos of last segment of previous chromosome
    if(grepl("X", seg_df["seqnames"]) || grepl("Y", seg_df["seqnames"]) || grepl("24", seg_df["seqnames"])){
        pos0CurrChr = sum(lengthOfChrs[1:22]) ## 22 because we exclude sex chromosomes
    }
    else {
        currChrArm = seg_df["seqnames"] 
        currChr = substring(currChrArm, 1, nchar(currChrArm)-1) # removing "p" or "q" char at the end of arm name, because we want the chromosome number
        # print(c("currChr in oncoscanR functions: ", currChr))
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

# generateGrid = function(graph_title, addAblines=F) {
#     #create empty plot to add things in
#     plot(1, ylim=c(0,3), xlim=c(0,3*10**9),col="white", xaxt="n", yaxt="n", xlab="nombre de bases", ylab="nombre de copies", main=graph_title)
#     if (addAblines) {
#         #  add horizontal grid
#         for (CN in c(-7:8)) {
#             abline(h=CN, col="dark grey")
#         }
#     }
#     # X-axis
#     axis(1, at = c(0, 5*10**8, 1*10**9, 1.5*10**9, 2*10**9, 2.5*10**9, 3*10**9))
#     # Y-axis
#     axis(2, at = c(-7:8))
# }

generateGrid = function(graph_title, mode="CN", addAblines=TRUE, addChrGrid=TRUE, ylim=c(-2,3)) {
    print("in oncosanR's generateGrid()'")
    if (mode=="CN") {
        #create empty plot to add things in
        plot(1, ylim=ylim, xlim=c(0,3*10**9),col="white", xaxt="n", yaxt="n", xlab="nombre de bases", ylab="nombre de copies", main=graph_title)
        #  add horizontal grid
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
    } else{ 
        if(mode=="LRR") {
            plot(1, ylim=c(-6,2), xlim=c(0,3*10**9),col="white", xlab="nombre de bases", ylab="Log Ratio", main=graph_title)
        }
    }
    if(addChrGrid) {
        lengthOfChrs = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
        abline(v=0, col="light grey", lwd=0.5, lty=20)
        for (chr in 1:length(lengthOfChrs)) {
                chrStart = sum(lengthOfChrs[1:chr])
                abline(v=chrStart, col="light grey", lwd=0.5, lty=20)
            }
    }
}



getPAAs = function(oncoscan_cov,armlevelAlter){
    currChrArm = oncoscan_cov["seqnames"]
    if (currChrArm %in% names(armlevelAlter)) {
        return(armlevelAlter[currChrArm])
    }
    return(0)
}
plotPAAWithThreshold = function(armLevelAlt_1type, plotTitle, oncoscanCoverage){
    OScovDf = as.data.frame(oncoscanCoverage)
    OScovDf$PAA_1type = apply(OScovDf, 1, getPAAs, armLevelAlt_1type)
    OScovDf = mutate(OScovDf, isTrueAlter=ifelse(OScovDf$PAA_1type > 0.8, "#b75869", "#5869b7")) # red, blue
    barplot(OScovDf$PAA_1type, names.arg = OScovDf$seqnames, ylim=c(0,1), main=plotTitle, col = "#5869b7")
    # barplot(OScovDf$PAA_1type, names.arg = OScovDf$seqnames, ylim=c(0,1), main=plotTitle, col = "#5869b7")
    # abline(h=0.8)
    # barplot(OScovDf$PAA_1type, names.arg = OScovDf$seqnames, ylim=c(0,1), main=plotTitle, col = OScovDf$isTrueAlter)
    # abline(h=0.8)
}

if (F) {
    ############################### custom run
    # oncoscanR::workflow_oncoscan.run(filepath, gender) :
    # sampleName = "5-LD"
    sampleName = "6-VJ"
    # sampleName = "8-MM"
    chas.fn = paste0("C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/premiers_E_recus/all_segmentFiles_from_ChAS/",sampleName,".OSCHP.segments.txt")
    # chas.fn = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_segmentFiles_from_ChAS/6-VJ_from_analysis_setup.segment.txt"
    gender = "M"

    # Remove the 21p arm from the Oncoscan coverage as it is only partly covered and we don't
    # want to return results on this arm.
    oncoscan.cov <- oncoscanR::oncoscan_na33.cov[seqnames(oncoscanR::oncoscan_na33.cov) != '21p']
    chrArmsPos = as.data.frame(oncoscan.cov)
    #### load ####
    # Load the ChAS file and assign subtypes.
    segments <- load_chas(chas.fn, oncoscan.cov)
    segments$cn.subtype <- get_cn_subtype(segments, gender)
    segsDf = as.data.frame(segments)

    #### clean ####
    # Clean the segments: resctricted to Oncoscan coverage, LOH not overlapping with copy loss
    # segments, smooth&merge segments within 300kb and prune segments smaller than 300kb.

    if(T){
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

    ### plotting ALA for slides
    OScovDf = as.data.frame(oncoscan.cov)



    ## displaying chromosomes arms in the normal order
    row.names(OScovDf) = OScovDf$seqnames
    rowsLayout = c("1p","1q", "2p","2q", "3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q","Xp", "Xq","Yp", "Yq")
    OScovDf_ordered = OScovDf[rowsLayout,]
    # creating color group
    # save_plots_here = file.path(resultsDir, sampleName)
    # if(!dir.exists(save_plots_here)) dir.create(save_plots_here)
    # barplot(OScovDf$PAA_loss, names.arg = OScovDf$seqnames, ylim=c(0,1), main=paste0(sampleName, " loss"), col = "#5869b7")
    # abline(h=0.8)
    # if (F) {
    plotPAAWithThreshold(armlevel.loss, paste0(sampleName, " loss"),oncoscan.cov)
    plotPAAWithThreshold(armlevel.loh, paste0(sampleName, " loh"),oncoscan.cov)
    plotPAAWithThreshold(armlevel.gain, paste0(sampleName, " gain"),oncoscan.cov)
    plotPAAWithThreshold(armlevel.amp, paste0(sampleName, " amp"),oncoscan.cov)
    # }

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

}

