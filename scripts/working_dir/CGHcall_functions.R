cleanProbeset = function(pathToMultiProbeset="allSamples_2_3_centeredProbeset.txt"){
    ######################## to clean all-samples probeset.txt file (it comes from ChAS analysis workflow):
    sampleNames_ChasOrder = c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )
    # sampleNames_ChasOrder = c("1-RV", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "2-AD", "20-CJ", "21-DC", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA")
    warning("check samplenames order before processing this block of code.")
    ## for loading data from a single probeset.txt file that contains all samples data
    allSamplesPath = file.path(dataDir, pathToMultiProbeset)
    ProbeData = read.table(allSamplesPath, sep='\t', h=T)
    ## remove weightedlog2Ratio, BAF and allelicDifference columns
    ProbeData_filtered = dplyr::select(ProbeData, -contains("Weighted"))
    ProbeData_filtered = dplyr::select(ProbeData_filtered, -contains("Backup")) #removing 5-LD_backup.OSCHP
    ProbeData_filtered = dplyr::select(ProbeData_filtered, 1:3, contains("Log2Ratio") )
    # sampleNamesByChas = colnames(ProbeData_filtered[4:length(ProbeData_filtered)])
    ## add END_POS column
    ProbeData_filtered$END_POS = ProbeData_filtered$Position ## previously ProbeData_filtered$END_POS = ProbeData_filtered$Position+20; 20 being the length of a probe (not so meaningful because probes provide information about SNPs). Tony said it is actually around 100.
    colnames(ProbeData_filtered)= c("probeID",  "CHROMOSOME", "START_POS", sampleNames_ChasOrder, "END_POS")
    ## order columns
    ProbeData_ordered = dplyr::select(ProbeData_filtered, c("probeID", "CHROMOSOME", "START_POS", "END_POS", all_of(sampleNames)))
    ## write table to file so we don't have to do this at every run
    allSamplesClean_path = file.path(dataDir, "allSamplesCleanProbeset_2_3Rec.txt")
    write.table(ProbeData_ordered, allSamplesClean_path, sep='\t', quote=F)
}



getSeg = function(currSampleSegs, s){
    ### name of value column must be "CN"
    i=s
    # print(c("s: ", s))
    # print(c("i: ", i))
    segToReturn = list()
    # segToReturn$currSegChr = chrVec[s]
    # print(c("currSegChr: ", currSegChr))
    segToReturn$currSegChr = currSampleSegs[s,]$Chromosome
    segToReturn$currSegVal = currSampleSegs$CN[s]
    # print(c("length(currSampleSegs[,1]): ", length(currSampleSegs[,1])))
    # print(c("currSampleSegs$CN[i+1]: ", currSampleSegs$CN[i+1]))
    # print(c("chrVec[i+1]: ", chrVec[i+1]))
    while((i<length(currSampleSegs[,1])) && (currSampleSegs$CN[i+1]==segToReturn$currSegVal) && (currSampleSegs[i+1,]$Chromosome==segToReturn$currSegChr)) {i=i+1}
    # print(c("chrVec[i+1]: ", chrVec[i+1]))
    # print(c("chrVec[i]: ", chrVec[i]))
    # print(c("chrVec[i-1]: ", chrVec[i-1]))
    segToReturn$currSegStart = currSampleSegs[s,]$Start
    segToReturn$currSegEnd = currSampleSegs[i,]$End
    segToReturn$currSegNbProbes = i-(s-1)
    segToReturn$i = i
    # print(c("class(currSegVal)", class(currSegVal)))
    # segToReturn = dplyr::select(segToReturn, c("currSegChr", "currSegStart", "currSegEnd", "currSegVal", "currSegNbProbes", "i"))
    segToReturn = segToReturn[c("currSegChr", "currSegStart", "currSegEnd", "currSegVal", "currSegNbProbes", "i")]
    # print(c("segToReturn: ", segToReturn))
    # print(c("one segment: ", currSegChr, currSegStart, currSegEnd, currSegVal, currSegNbProbes, i))
    # return(list(currSegChr, currSegStart, currSegEnd, currSegVal, currSegNbProbes, i))
    return(segToReturn)
}

getSegTable = function(currSampleSegs) {
    ### name of value column must be "CN"
    ## initialize by-segments table 
    segTableBySegment = data.frame(matrix(ncol = 5, nrow = 0))
    colnames(segTableBySegment) = c("Chromosome", "Start", "End", "Value", "nbProbes")
    segId = 1
    s=1 # start of segment
    i=1 # end of segment
    # print(c("currSampleSegs: ", currSampleSegs)) 
    # print(c("length(currSampleSegs[,1]): ", length(currSampleSegs[,1])))
    while (s < length(currSampleSegs[,1])) {
        print(paste0("------------- segId = ", segId," -------------"))
        resSeg = getSeg(currSampleSegs, s)
        # print(c("resSeg[1:5]: ", resSeg[1:5]))
        i = resSeg$i
        segTableBySegment[segId,] = resSeg[1:5]
        # print(c("resSeg: ", resSeg))
        s = resSeg[[length(resSeg)]]+1
        segId = segId+1
        # print(c("segId: ", segId))
    }
    return(segTableBySegment)
}



getSegTables = function(segTableByProbe, sampleNames) {
    # get segTables for all samples; concatenate them in a list
    segTablesList = list()
    for(sample in sampleNames) {
        print(paste0("==================================== sample = ", sample," ===================================="))
        currSample_SegTableByProbe = cbind(rowsInfo, segTableByProbe[[sample]])
        colnames(currSample_SegTableByProbe) = c(colnames(rowsInfo), "CN")
        
        print(c("currSample_SegTableByProbe: ", currSample_SegTableByProbe))
        currSegTable = getSegTable(currSample_SegTableByProbe)
        segTablesList = append(segTablesList, list(currSegTable))
    }
    return(segTablesList)
}

plotSegTable = function(currSegTable,currSampleName,resultsDir) {
    print(c("plotting sample ", currSampleName))
    imgName = paste0(resultsDir,"/",currSampleName,"segsUsedForGI.png")
    png(imgName, width=700, height=484)
    generateGrid(paste0(currSampleName, " Copy number"), mode="CN")
    colnames(currSegTable) = c("chrom", "loc.start", "loc.end", "callVal", "nbProbes")
    apply(currSegTable, 1, plotSeg_rCGH, "callVal", indivSeg=TRUE)
    dev.off()
}

plotSegTables = function(segTablesList, sampleNames, resultsDir) {
    if(!dir.exists(resultsDir)) dir.create(resultsDir)
    for (i in 1:length(sampleNames)) {
        currSegTable = segTablesList[[i]]
        currSampleName = sampleNames[[i]]
        plotSegTable(currSegTable,currSampleName,resultsDir)
    }    
}


getPrbLvSegmentsFromCallObj = function(callRes) {
    ### callRes must be a cghCall object containing results of one or more samples, but can't be a list of cghCall objects
    sampleNames = colnames(callRes@assayData[["calls"]])
    rowsInfo = fData(callRes)
    CGHcall_segments = as.data.frame(calls(callRes))
    CGHcall_segments = cbind(rowsInfo, CGHcall_segments)
    colnames(CGHcall_segments) = c(colnames(rowsInfo), sampleNames)
    return(CGHcall_segments)
}


getPrbLvSegments = function(pipelineRes) {
    ### this works whatever the input is, as long as the pipeline ran correctly.
    if(is.list(pipelineRes)) {
        print("result is a list of individual cghCall objects")
        currRes = pipelineRes[[1]]
        probeLevelSegments = getPrbLvSegmentsFromCallObj(currRes)
        for(i in 2:length(pipelineRes)) { ## I must make it so a list containing only 1 cghCall result can't exist
            currRes = pipelineRes[[i]]
            curr_prbLevSegs = getPrbLvSegmentsFromCallObj(currRes)
            probeLevelSegments = cbind(probeLevelSegments, curr_prbLevSegs[,length(curr_prbLevSegs)])
            colnames(probeLevelSegments)[length(probeLevelSegments)] = colnames(curr_prbLevSegs)[length(colnames(curr_prbLevSegs))]
            # print(c("probeLevelSegments: ", probeLevelSegments))
        }
    } else if(class(pipelineRes)=="cghCall") {
            print("result is a single cghCall object.")
            probeLevelSegments = getPrbLvSegmentsFromCallObj(pipelineRes)
    } else {
        print("invalid input")
    }
    return(probeLevelSegments)
}




#################################### calculate GI
getNbChrsCGHcall = function(segmentsTable) { # function from ASCAT.R
    chrs = as.vector(segmentsTable$Chromosome)
    print(c("chrs: ", chrs))
    print(c("chrs: ", unique(chrs)))
    nbChr = length(unique(chrs))
    print(c("nbChr: ", nbChr))
    return(nbChr)
}

calcGI_CGHcall = function(segmentsTable) {
    ## removing segments of log ratio 0 as they are not aberrations
    segmentsTableClean = dplyr::filter(segmentsTable, Value!=0)
    nbOrigin = dim(segmentsTable)[1]
    nbClean = dim(segmentsTableClean)[1]
    nbRemoved = nbOrigin-nbClean
    print(paste0(nbRemoved, " segments removed out of ", nbOrigin))

    nbChr = getNbChrsCGHcall(segmentsTableClean)
    nbAlter = dim(segmentsTableClean)[1]
    print(c("nbAlter: ", nbAlter))
    print(c("nbChr: ", nbChr))
    GI = calcGI(nbAlter, nbChr)
    return(list(GI, nbAlter, nbChr))
}




## pour visualiser les sondes du groupe -1 du plot precedent (je parle de calling) au sein du plot logratios_for_hist
# allProbes_niveauPlusUn = c()
getProbesOfCallLevelfromOtherDf = function(lvl, df_calls, df_logCN) {
    nb_cols = dim(df_calls)[2]
    allProbes_GivenLevel = c()
    for(i in 1:nb_cols) {
        print("~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*new iteration*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~")
        #1 obtenir la liste des noms de ces sondes (et le sample auquel elles appartiennent)
            #1.1 extraire les 5 cols de df_calls.
            call_currSample = df_calls[,i]
            samplename = colnames(df_calls)[i]
            #1.2 pour chacune, donner la valeur TRUE aux cases qui ont la valeur 1. les autres cases recoivent la valeur FALSE.
            summaryAnyVec(call_currSample)
            currSample__call_equals_lvl = call_currSample==lvl
            summaryAnyVec(currSample__call_equals_lvl)
        #2 extraire les valeurs de ces sondes A partir du dataframe logratios CN.
            #2.1 ajouter la colonne call_equals_lvl au df logratios CN.
            df_logCN = cbind(df_logCN, currSample__call_equals_lvl)
            colnames(df_logCN)[i+nb_cols] = paste(samplename, "__call_equals_lvl")
            #2.2 filtrer les lignes pour lesquelles la colonne call_equals_lvl vaut TRUE
            df_logCN = as.data.frame(df_logCN)
            filteredOnCurrSample_callEqualsLvl = dplyr::filter(df_logCN, df_logCN[i+5]==TRUE)
            #2.3 extraire la colonne de notre sample de ce sous-tableau
            CurrSample__callEquals1 = filteredOnCurrSample_callEqualsLvl[,1]
            cl = class(CurrSample__callEquals1)
            # View(CurrSample__callEquals1)
            print(c("class of appended result: ", cl))
        #3 ajouter ces sondes a la liste
        allProbes_GivenLevel = append(allProbes_GivenLevel, CurrSample__callEquals1)
        # allProbes_niveauPlusUn = append(allProbes_niveauPlusUn, CurrSample__callEquals1)
    }
    return(allProbes_GivenLevel)
}

getProbesOfCallLevel = function(lvl, df_calls) {
    nb_cols = dim(df_calls)[2]
    allProbes_GivenLevel = c()
    for(i in 1:nb_cols) {
        print("~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*new iteration*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~")
        call_currSample = df_calls[,i]
        # summaryAnyVec(call_currSample)
        call_currSample = as.data.frame(call_currSample)
        currSample__call_equals_lvl = dplyr::filter(call_currSample, call_currSample==lvl)
        currSample__call_equals_lvl__asVec = unlist(currSample__call_equals_lvl)
        # summaryAnyVec(currSample__call_equals_lvl__asVec)
        allProbes_GivenLevel = append(allProbes_GivenLevel, currSample__call_equals_lvl__asVec)
    }
    return(allProbes_GivenLevel)
}










## to understand postsegnormalize

customPostsegnormalize = function(segmentData, inter = c(-0.1, 0.1)) 
{
    seg <- segmented(segmentData)
    values <- c()
    for (i in 1:ncol(seg)) {
        values <- c(values, median(seg[, i]))
    }
    matrixValues <- matrix(rep(values, nrow(seg)), ncol = ncol(seg), 
                           byrow = TRUE)
    seg <- seg - matrixValues
    countlevall <- apply(seg, 2, function(x) {
        as.data.frame(table(x))
    })
    print(c("countlevall aka segvec: ", countlevall))
    intcount <- function(int, sv) {
        print("OO======================intcount=======================OO")
        print(c("int: ", int))
        # print(c("sv: ", sv)) #table with two columns: the probe log ratio value, and the number of times it appears (== length of 1 segment)
        sv1 <- as.numeric(as.vector(sv[, 1])) #sv1: log ratio column
        wh <- which(sv1 <= int[2] & sv1 >= int[1]) # position of probes which value is within the interval
        # print(c("sv1: ", sv1))
        # print(c("class(sv1): ", class(sv1)))
        sv1 = as.data.frame(sv1)
        print(c("wh: ", wh)) 
        print(c("count of probes within interval: ", (sv[wh, 2])))
        print(c("values of probes within interval: ", (sv1[wh,])))
        print(c("returned_sum: ", sum(sv[wh, 2])))
        return(sum(sv[wh, 2])) # the number of probes(not their value) within the interval is returned
    }
    postsegnorm <- function(segvec, int = inter, intnr = 3) { #int = c(-0.5, 0.3)
        print("OO=============================================postsegnorm================================================OO")
        # print(c("interval received: ", int))
        intlength <- (int[2] - int[1])/2
        gri <- intlength/intnr
        intst <- int[1] + (0:intnr) * gri # e.g. (-0.50000000 -0.23333333  0.03333333  0.30000000)
        intend <- intst + intlength # e.g. (0.3000000 0.5666667 0.8333333 1.1000000)
        ints <- cbind(intst, intend) # intervals
        # print(c("intst: ", intst))
        # print(c("intend: ", intend))
        # print(c("intervals list 1/2: ", ints))
        intct <- apply(ints, 1, intcount, sv = segvec) # for each interval, finds the probes contained in it then counts them. we often end up with around 3 logR values repeated each dozens/hundreds of times. The total number of probes is then counted for this interval.
        # intct contains then one value per interval, which represents how much segmented the data is.
        whmax <- which.max(intct) # finds highest value, representing the best interval.
        # print(c("segvec: ", segvec))
        print(c("interval count: ", intct))# of all
        print(c("whmax: ", whmax))
        print(c("ints: ", ints))
        print(c("best interval found (ints[whmax, ]): ", ints[whmax, ]))
        return(ints[whmax, ]) #returns said best interval.
    }
    postsegnorm_rec <- function(segvec, int, intnr = 3) {
        print("============== postsegnorm_rec =============")
        newint <- postsegnorm(segvec, int, intnr)
        newint <- postsegnorm(segvec, newint, intnr)
        newint <- postsegnorm(segvec, newint, intnr)
        newint <- postsegnorm(segvec, newint, intnr)
        newint <- postsegnorm(segvec, newint, intnr)
        return(newint[1] + (newint[2] - newint[1])/2) # we return the middle point of the interval.
    }
    listres <- lapply(countlevall, postsegnorm_rec, int = inter) # this runs postsegnorm_rec once for every sample in our seg dataset.
    vecres <- c()
    for (i in 1:length(listres)) {
        vecres <- c(vecres, listres[[i]])
    }
    print(c("listres: ", listres))
    print(c("vecres: ", vecres))
    segmented(segmentData) <- t(t(seg) - vecres) # this substracts vecres[i](the middle point of the interval found) to all segmentation data, hence normalizing.
    copynumber(segmentData) <- t(t(copynumber(segmentData) - 
                                       matrixValues) - vecres)
    return(segmentData)
}


#function call
if (F) {
        
    postseg.cghdata <- customPostsegnormalize(seg.cghdata) # argument: objet cghSeg
    # plot before/after postsegnormalize
    png("./plots/before_postsegnorm.png")
    sampleNames(seg.cghdata) = "segmented data"
    plot(seg.cghdata, ylim = c(-2,2), ylab="log ratio", xlab="position genomique", main="segmented data")
    dev.off()
    png("./plots/after_postsegnorm.png")
    sampleNames(postseg.cghdata) = "segmented data after normalization"
    plot(postseg.cghdata, ylim = c(-2,2), ylab="log ratio", xlab="position genomique", main="segmented data after normalization")
    dev.off()

    segTable = segmented(seg.cghdata)
    # function for plotting 
    plotInterval = function(data, plot_ylim, interval) {
        png(paste("./plots/postsegnorm_recherche_intervalle_", toString(interval), ".png", sep = ""))
        try(plot(data, ylim=plot_ylim, main="recherche du meilleur intervalle", ylab="log ratio", xlab="position genomique")) #,pch=1, cex=0.2
        abline(h=interval[1], col="#3b2dbc"); abline(h=interval[2], col="#3b2dbc")
        abline(h=plot_ylim[1], col="#3b2dbc"); abline(h=plot_ylim[2], col="#3b2dbc") #59, 45, 188
        nb_points = length(data)
        x = -10000:nb_points+10000
        polygon(x=c(1, nb_points, nb_points, 1), y=c(interval[1],interval[1],interval[2],interval[2]), col = "#65BFFF", density = 10)
        dev.off()
    }
    prevInter = c(-2.5,1) # c(-2.5,1) is default ylim for this data
    currInter = c(-0.1,0.1)
    plotInterval(segTable,prevInter,currInter) 
    prevInter = currInter
    currInter = c(-0.0333,0.0667) 
    plotInterval(segTable,prevInter,currInter) 
    prevInter = currInter
    # currInter = c(-0.0167,0.0333) # this is the actual interval but it does not appear to have more segments in it in the plot.
    currInter = c(0.0167,0.0667) # this is visually the best interval, I use it for visualisation purposes
    plotInterval(segTable,prevInter,currInter) 
    prevInter = currInter
    currInter = c(-0.0083,0.0167)
    plotInterval(segTable,prevInter,currInter) 
    prevInter = currInter
    currInter = c(-0.0083,0.0042)
    plotInterval(segTable,prevInter,currInter) 
    prevInter = currInter
    currInter = c(-0.0042,0.0021)
    plotInterval(segTable,prevInter,currInter) 

    intnr=3
    int = c(-0.5, 0.3)
    intlength = 0.8
    gri= intlength/intnr
    intst = -0.5 + (0:intnr) * gri
    intend = intst + intlength
    print(intst)
    print(intend)


    x1 = -0.033333
    x2 = 0.666667

    x1 + (x2-x1)/2
}

