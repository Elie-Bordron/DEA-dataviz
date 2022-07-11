## ----Functions -----------------------------------------------------
removePointsForQuickPlotting = function(cghDf, pointsToRemove=20) {
    # if pointsToRemove=5, we keep 1 value out of 5.
    cghDfFiltered = dplyr::slice(cghDf, seq(1,length(cghDf[,1]),pointsToRemove))
    return(cghDfFiltered)
}


###########These functions are now in CGHcall.R
# getAbspos_probe = function(probeSet, lengthOfChrs) {
#     # Check that it doesn't match any non-number
#     numbers_only <- function(x) !grepl("\\D", x)
#     print(c("probeSet: ", probeSet))
#     # print(c("probeSet[\"ChrNum\"]: ", probeSet["ChrNum"]))
#     if(numbers_only(probeSet["ChrNum"])) {
#         currentChr = as.numeric(probeSet["ChrNum"])
#     } else {
#         ## this works is string is like "chr13"; doesn't work for i.e. "Chr13"
#         currentChr = stringr::str_split(probeSet["ChrNum"], "chr")[[1]]
#         currentChr = as.numeric(currentChr[[2]])
#     }
#     currentChr
#     if (currentChr==1){
#         pos0CurrChr=0
#     } else {
#         pos0CurrChr = sum(lengthOfChrs[1:currentChr-1])
#     }
#     newPos = pos0CurrChr + as.numeric(probeSet["ChrStart"])
#     return(newPos)
# }

# getAbspos_probeset = function(probeSet) {
#     lengthOfChrs = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
#     probeSet$absPos = apply(probeSet, 1, getAbspos_probe, lengthOfChrs)
#     return(probeSet)
# }

plotSeg_rCGH = function(seg_df, value_col, indivSeg=FALSE, segColor="dark blue", alreadyGoodPos=FALSE){
    ########### Used by rCGH.R and CGHcall.R ########### 
    # print(c("value_col: ", value_col))    
    lengthOfChrs = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
    ## get endpos of last segment of previous chromosome
    if((seg_df["chrom"]=="X" ) || (seg_df["chrom"]=="Y") || (seg_df["chrom"]==24) || (seg_df["chrom"]==25)){
        pos0CurrChr = sum(lengthOfChrs[1:22]) ## 22 because we exclude sex chromosomes
    }else {
        # print(c("seg_df['chrom']: ", seg_df["chrom"]))
        # print(c("class(seg_df['chrom']): ", class(seg_df["chrom"])))
        currChromosome = seg_df["chrom"]
        currChromosome = as.numeric(currChromosome)
        # print(c("currChromosome: ", currChromosome))
        if (currChromosome!=1){
            pos0CurrChr = sum(lengthOfChrs[1:currChromosome-1])
        } else {
            pos0CurrChr=0
        }
    }
    ## drawing a segment on plot for each segment of the genome, using estimated values
    segStartPos = pos0CurrChr + as.numeric(seg_df[["loc.start"]])
    segEndPos = pos0CurrChr + as.numeric(seg_df[["loc.end"]])
    if(alreadyGoodPos) {
        segStartPos = as.numeric(seg_df[["loc.start"]])
        segEndPos = as.numeric(seg_df[["loc.end"]])
    }
    # estimCN = as.numeric(seg_df[["Log2Ratio"]])
    # estimCN = as.numeric(seg_df[["estimCopy"]])
    # print(c("seg_df[[value_col]]: ", seg_df[[value_col]]))
    estimCN = as.numeric(seg_df[[value_col]])
    # print(c("pos0CurrChr, segStartPos: ", pos0CurrChr+segStartPos))
    segments(segStartPos, estimCN, segEndPos, estimCN, col=segColor, lwd=2)
    if(indivSeg) {
        height = 0.05
        segments(segStartPos, estimCN+height, segStartPos, estimCN-height, col=segColor, lwd=0.1)
        segments(segEndPos, estimCN+height, segEndPos, estimCN-height, col=segColor, lwd=0.1)
    }
}

getPrbLvSegments_rCGH = function(cghNorm) {
    ### extract LRR & CN
    print(c("cghNorm@cnSet: ", cghNorm@cnSet))
    # print(c("class(cghNorm@cnSet): ", class(cghNorm@cnSet)))
    cghNorm@cnSet = cghNorm@cnSet[c("ProbeName", "ChrNum", "ChrStart", "Log2Ratio", "Segm", "estimCopy")] # absPos not needed for get_seg_table()
    print(c("cghNorm@cnSet after selecting columns: ", cghNorm@cnSet))
    colnames(cghNorm@cnSet) = c("probeID", "Chromosome", "Start", "absPos", "rawLRR", "Log2Ratio", "CN")
    print(c("cghNorm@cnSet after renaming columns: ", cghNorm@cnSet))
    cghNorm@cnSet
}


hush=function(code){ ## function to silence another function's prints while still returning its output
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
}


#################################### calculate GI

getNbChrs = function(segmentsTable) { # function from ASCAT.R
    chrs = as.vector(segmentsTable$chrom)
    nbChr = length(unique(chrs))
    return(nbChr)
}

calcGI_rCGH = function(segmentsTable) {
    ## removing segments of copy number 2 as they are not aberrations (actually,they could be Copy Number-neutral events, but rCGH can't detect these)
    segmentsTable = dplyr::filter(segmentsTable, estimCopy!=2)
    nbChr = getNbChrs(segmentsTable)
    nbAlter = dim(segmentsTable)[1]
    print(c("nbAlter: ", nbAlter))
    print(c("nbChr: ", nbChr))
    print(c("nbAlter**2: ", nbAlter**2))
    print(c("(nbAlter**2)/nbChr: ", (nbAlter**2)/nbChr))
    GI = calcGI(nbAlter, nbChr)
    return(list(GI,nbAlter,nbChr))
}

