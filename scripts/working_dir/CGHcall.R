## set working directory from shiny app
if(exists("working_dir_shiny")) {
    working_dir = working_dir_shiny
} else {
    ## bergo
    # working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
    ## C
    working_dir = "C:/Users/warew/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
    
}


if (FALSE) {
    setwd(working_dir)
    options("max.print"=100)
    ## open working directory in Files tab
    rstudioapi::filesPaneNavigate(working_dir)
}

## import libraries
library(dplyr)
library(CGHcall)

pipelineCGHcall = function(osData, params) {
    # osData = rawProbesData
    # print(c("osData: ", os    Data))

    before=Sys.time()
    # osData = s2Probes # to run on one sample
    if(is.null(params$tumor_prop)) {params$tumor_prop=1}
    # ACGH_data <- make_cghRaw(Wilting)
    ACGH_data <- make_cghRaw(osData)
    # we want to apply fewest changes possible to data, so we want to do our own preprocess if we have time
    # cghdata = removedNaNProbes = dplyr::filter(ACGH_data, !is.na(ACGH_data[5]))
    print("- preprocess -")
    cghdata <- preprocess(ACGH_data, maxmiss=params$Maxmiss, nchrom=22) # because we don't need sex chromosomes data for GI.
    # plot(cghdata)
    print("- normalize -")
    norm.cghdata <- normalize(cghdata, method=params$NormMethod, smoothOutliers=params$SmoothOutliers)
    print("- segment -")
    seg.cghdata <- segmentData(norm.cghdata, method="DNAcopy", undo.splits="sdundo",undo.SD=params$UndoSD, clen=params$Clen, relSDlong=params$RelSDlong)
    # segTable = segmented(seg.cghdata)
    print("- postsegnormalize -")
    postseg.cghdata <- postsegnormalize(seg.cghdata, inter=params$Inter)
    # plot(postseg.cghdata, ylimit=c(-2,2))
    # print(c("params$tumor_prop,params$Prior, params$Robustsig, params$Minlsforfit: ", params$tumor_prop,params$Prior, params$Robustsig, params$Minlsforfit))
    print("- call -")
    rawCghResult <- CGHcall(postseg.cghdata,nclass=5,cellularity=params$tumor_prop, prior=params$Prior, robustsig = params$Robustsig, minlsforfit = params$Minlsforfit)
    print("- expand call -")
    CghResult <- ExpandCGHcall(rawCghResult,postseg.cghdata,CellularityCorrectSeg=params$CellularityCorrectSeg) # use CellularityCorrectSeg=TRUE to correct using cellularity
    after = Sys.time()
    CghResult$processingTime = round(as.numeric(difftime(after, before, units = "secs"))/dim(CghResult)[2], 2) ## If more than 1 sample is processed for this run, the average value is given to each sample. This value is stored in CghResult@phenoData@data[["processingTime"]].
    
    return(CghResult)
    # return(list(ACGH_data, cghdata, norm.cghdata, seg.cghdata, postseg.cghdata, rawCghResult, CghResult))
}

getDefParams = function() {
    params = list()
    ## run parameters
    params$date = Sys.time()
    params$runName = "Default run"
    params$runAsCohort = FALSE
    # params$sampleNames = c("11-BG") # for single sample run
    # params$sampleNames = c("1-RV") # for single sample run
    params$saveAsValidGI = FALSE
    params$saveSamePkg = FALSE
    ## pipeline parameters
    params$Maxmiss=0.95
    params$NormMethod = "median"
    params$SmoothOutliers = TRUE
    params$UndoSD=3
    params$Clen=10
    params$RelSDlong=5
    params$Inter = c(-0.1,0.1)
    params$Prior = "not all"
    params$Robustsig = "no"
    params$Minlsforfit = 0.01 # Minimum length of the segment (in Mb) to be used for fitting the model 
    params$CellularityCorrectSeg = T
    params$tumor_prop = 1
    return(params)
}








main = function() {
    ## setting paths
    # dataDirProbesets = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/premiers_E_recus/all_probeset"
    ## only true on Cass
    dataDirProbesets = working_dir_shiny
    CGHcallDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/CGHcall"
    resultsDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/CGHcall"
    ######################## Define parameters
    params = getDefParams()
    params$sampleNames = "1-RV"
    # params$sampleNames = c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )
    

    ######################## start analysis
    if(FALSE){ ## set this to TRUE to use one file for all samples
        ## load all-samples probeset.txt file
        inputFile = "allSamplesCleanProbeset_2_3Rec.txt"
        CleanInputFile_path = file.path(dataDirProbesets, inputFile)
        params$tumor_prop = c(0.9,0.9,0.9,0.9,0.9,0.8,0.9,0.8,0.9,0.8,0.8,0.9,0.8,0.8,0.8,0.9,0.8,1,0.95,0.8,0.8) # samples 1 to 21
        params$saveSamePkg = FALSE
        rawProbesData = read.table(CleanInputFile_path, h=T)
        colnames(rawProbesData) = c("probeID", "CHROMOSOME", "START_POS", "END_POS", c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" ))
    } else {
        source(file.path(working_dir, "CGHcall_functions.R"))
        rawProbesData = loadProbeset(dataDirProbesets)
        rawProbesData = nameprobeset(rawProbesData)
    }
    
    ### Run pipeline
    if(params$runAsCohort) {
        callAllSamples = pipelineCGHcall(rawProbesData, params)
    } else {
        callAllSamples = list()
        tumor_prop = params$tumor_prop
        for (s in 1:length(params$sampleNames)) {
            currSampleName = params$sampleNames[s]
            print(c("currSampleName: ", currSampleName))
            currProbesData = dplyr::select(rawProbesData, c("probeID", "CHROMOSOME", "START_POS", "END_POS", all_of(currSampleName)))
            params$tumor_prop = tumor_prop[s]
            currCallRes = pipelineCGHcall(currProbesData, params)
            callAllSamples = append(callAllSamples, currCallRes)
        }
        print(c("callAllSamples: ", callAllSamples))
        if(length(callAllSamples)==1) {
            callAllSamples = callAllSamples[[1]] ## if only one sample was treated, we don't want it to be in a list.
        }
    }

    ############### convert probes table to segments tables
    ## retrieve call segments
    CGHcall_segments = getPrbLvSegments(callAllSamples, segsType="both")
    
    # ## get segments tables
    # # source(file.path(working_dir, "rCGH_functions.R"))
    # # source(file.path(working_dir, "oncoscanR_functions.R"))
    # source(file.path(working_dir, "CGHcall_functions.R"))
    # allSegTables = getSegTables(CGHcall_segments,params$sampleNames)
    # probeData = removePointsForQuickPlotting(rawProbesData)
    # colnames(probeData)[c(2:3)] = c("ChrNum", "ChrStart")
    # probeData = getAbspos_probeset(probeData)
    # colnames(probeData)[c(2:3)] = c("CHROMOSOME", "START_POS")
    # currSegTable = allSegTables[[1]]
    # # colnames(currSegTable)[6] <- "Log2Ratio"
    # currSampleName = params$sampleNames
    # colnames(probeData)[which(colnames(probeData)==currSampleName)] <- "Log2Ratio"
    # plotSegTableForWGV_GG(currSegTable, probeData)
    # # colnames(currSegTable)[6] <- "Value"
    # colnames(probeData)[which(colnames(probeData)=="Log2Ratio")] <- currSampleName
    # ## plot called data on all profiles
    # # plotSegTables(allSegTables,params$sampleNames,resultsDir)
    # 
    # 
    # 
    # 
    # 
    
    #--------------------------- test here, delete this section asap
    source(file.path(working_dir, "CGHcall_functions.R"))
    CGHcall_segments = getPrbLvSegments(callAllSamples, segsType="both")
    segTableByProbe = prepareSegtableByProbe(CGHcall_segments)
    segTable = get_seg_table(segTableByProbe)
    GI = calcGI_CGHcall(segTable)
    GI
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############### compute GI
    ## initialize GI df
    GI_CGHcall_df = data.frame(matrix(ncol = 4, nrow = length(params$sampleNames)))
    colnames(GI_CGHcall_df) = c("GI", "nbAlter", "nbChr", "runTime")
    rownames(GI_CGHcall_df) = params$sampleNames
    for (s in 1:length(allSegTables)) {
        print(paste0("======= sample ", params$sampleNames[s], " ======="))
        if("Value" %in% colnames(allSegTables[[s]])){ print(" a Value col is present. perfect")}else{warning(" need column named Value. Try again, or GI will be miscalculated")}
        GI_res = calcGI_CGHcall(allSegTables[[s]])
        print(c("GI_res: ", GI_res))
        if(params$runAsCohort) {
            GI_CGHcall_df[params$sampleNames[s],] = append(GI_res, callAllSamples[,s]$processingTime)
        } else {
            GI_CGHcall_df[params$sampleNames[s],] = append(GI_res, callAllSamples[[s]]$processingTime)
        }
    }
    
    ############### saving segmentation tables
    saveSegTables = function(allSegTables, CGHcallDir, sampleNames = c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC")) {
        for (s in 1:length(allSegTables)) {
            currSample = sampleNames[s]
            currSegTablePath = paste0(CGHcallDir,"/",currSample,"_segtable.txt")
            write.table(allSegTables[s], currSegTablePath, sep="\t", row.names=FALSE, quote=F)
        }
    }
    segTablesDir = file.path(CGHcallDir, "segTables")
    if(!dir.exists(segTablesDir)) dir.create(segTablesDir)
    saveSegTables(allSegTables, segTablesDir)
    


    
    ############### saving GI table
    ### along GIs of other packages
    if(saveAsValidGI) {
        source(file.path(working_dir, "crossPackagesFunctions.R"))
        saveGI_ResToFile(GI_CGHcall_df, "CGHcall")
    }

    ### or along other runs of the same package
    if(saveSamePkg) {
        CGHcallOutputFile = file.path(CGHcallDir, "GI_res.txt")
        CGHcallOutputFile
        saveGI_ResToFile(GI_CGHcall_df, runName, CGHcallOutputFile)
    }
    ### write log file
    x = capture.output(params)
    x[1] = paste0("\n==================================================\n", x[1])
    logText = x
    logfile = file.path(CGHcallDir, "log.txt")
    write(logText, logfile, append=T)
}











############ pour slides soutenance
if (FALSE) {
    soutenance_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/res_soutenance"
    # ylim = c(-3,3)
    ylim = c(-1,1)
    PercentRowsToRemove = 30
    xlab = "position genomique"
    pkgName="CGHcall"
    ###### plot dimensions:
    pngWidth = 350
    pngHeight = 300
    sampleName = params$sampleNames[1]

    ## result of pipeline
    ACGH_data = currCallRes[[1]]
    cghdata = currCallRes[[2]]
    norm.cghdata = currCallRes[[3]]
    seg.cghdata = currCallRes[[4]]
    postseg.cghdata = currCallRes[[5]]
    rawCghResult = currCallRes[[6]]
    CghResult = currCallRes[[7]]


    ### plot raw value per probe 
    state = "raw"
    vals = as.data.frame(ACGH_data@assayData[["copynumber"]])
    vals = removePointsForQuickPlotting(vals, PercentRowsToRemove)
    savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
    png(savePath, width = pngWidth, height = pngHeight)
    plot(vals[[sampleName]], pch = 20, cex = 0.01, ylab = 'log Ratio', xlab = xlab, ylim = ylim)
    dev.off()

    ### After preprocess
    state="preprocess"
    vals = as.data.frame(cghdata@assayData[["copynumber"]])
    vals = removePointsForQuickPlotting(vals, PercentRowsToRemove)
    savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
    png(savePath, width = pngWidth, height = pngHeight)
    plot(vals[[sampleName]], pch = 20, cex = 0.01, ylab = 'log Ratio', xlab = xlab, ylim = ylim)
    dev.off()

    ### After norm 1
    state="norm1"
    vals = as.data.frame(norm.cghdata@assayData[["copynumber"]])
    vals = removePointsForQuickPlotting(vals, PercentRowsToRemove)
    savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
    png(savePath, width = pngWidth, height = pngHeight)
    plot(vals[[sampleName]], pch = 20, cex = 0.01, ylab = 'log Ratio', xlab = xlab, ylim = ylim)
    dev.off()

    ### After seg
    state="seg"
    vals = as.data.frame(seg.cghdata@assayData[["segmented"]])
    vals = removePointsForQuickPlotting(vals, PercentRowsToRemove)
    savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
    png(savePath, width = pngWidth, height = pngHeight)
    plot(vals[[sampleName]], pch = 20, cex = 0.01, ylab = 'log Ratio', xlab = xlab, ylim = ylim)
    dev.off()

    ### After postsegnorm
    state="postsegnorm"
    vals = as.data.frame(postseg.cghdata@assayData[["segmented"]])
    vals = removePointsForQuickPlotting(vals, PercentRowsToRemove)
    savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
    png(savePath, width = pngWidth, height = pngHeight)
    plot(vals[[sampleName]], pch = 20, cex = 0.01, ylab = 'log Ratio', xlab = xlab, ylim = ylim)
    dev.off()

    ### plot Call value per probe
    state = "call"
    vals = as.data.frame(CghResult@assayData[["calls"]])
    vals = removePointsForQuickPlotting(vals, PercentRowsToRemove)[[1]]
    vals = vals + 2 # values are given as a normal CN == 0 , a single loss == -1, etc. For consistency, we make it like this: normal==2, single loss==1, etc.
    savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
    png(savePath, width = pngWidth, height = pngHeight)
    plot(y = vals, x=1:length(vals), pch=20, cex = 0.1, ylab = 'Nombre de copies', xlab = xlab)
    dev.off()

}





if (sys.nframe() == 0){
    source(file.path(working_dir, "rCGH.R"))
    source(file.path(working_dir, "rCGH_functions.R"))
    source(file.path(working_dir, "CGHcall_functions.R"))
    source(file.path(working_dir, "crossPackagesFunctions.R"))
    source(file.path(working_dir, "oncoscanR_functions.R"))
    main()
}
























