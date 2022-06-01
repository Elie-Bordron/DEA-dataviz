## ----Functions -----------------------------------------------------
removePointsForQuickPlotting = function(cghDf, pointsToRemove=20) {
    # if pointsToRemove=5, we keep 1 value out of 5.
    cghDfFiltered = dplyr::slice(cghDf, seq(1,length(cghDf[,1]),pointsToRemove))
    return(cghDfFiltered)
}

getNewPos_iterative = function(cghDf, lengthOfChrs) {
    currentChr = as.numeric(cghDf["ChrNum"])
    # print(c("currentChr: ", currentChr))
    if (currentChr!=1){
        pos0CurrChr = sum(lengthOfChrs[1:currentChr-1])
    } else {
        pos0CurrChr=0
    }
    newPos = pos0CurrChr + as.numeric(cghDf["ChrStart"])
    return(newPos)
}

getNewPos = function(cghDf) {
    lengthOfChrs = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
    cghDf$absPos = apply(cghDf, 1, getNewPos_iterative, lengthOfChrs)
    return(cghDf)
}

plotSeg_rCGH = function(seg_df, value_col, indivSeg=FALSE){
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
    segStartPos = as.numeric(seg_df[["loc.start"]])
    segEndPos = as.numeric(seg_df[["loc.end"]])
    # estimCN = as.numeric(seg_df[["Log2Ratio"]])
    # estimCN = as.numeric(seg_df[["estimCopy"]])
    # print(c("seg_df[[value_col]]: ", seg_df[[value_col]]))
    estimCN = as.numeric(seg_df[[value_col]])
    # print(c("pos0CurrChr, segStartPos: ", pos0CurrChr+segStartPos))
    segments(pos0CurrChr+segStartPos, estimCN, pos0CurrChr+segEndPos, estimCN, col="black", lwd=2)
    if(indivSeg) {
        segments(pos0CurrChr+segStartPos, estimCN, pos0CurrChr+segEndPos, estimCN, col="black", lwd=2)

    }
}


hush=function(code){ ## function to silence another function's prints while still returning its output
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
}

# foo=function(){
#     print("BAR!")
#     return(42)
# }

# x = hush(foo())

#################################### calculate GI

getNbChrs = function(segmentsTable) { # function from ASCAT.R
    chrs = as.vector(segmentsTable$chrom)
    nbChr = length(unique(chrs))
    return(nbChr)
}

calcGI_rCGH = function(segmentsTable) {
    ## removing segments of copy number 2 as they are not aberrations (actually,they could be Copy Number-neutral events, but rCGH can't detect these')
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

