##### this is the guide to rCGH package. use this with the pdf in docs

### run this when starting R session 
## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)
## loading libraries
library(rCGH)
library(dplyr)

#### load rCGH objects from .RDS files
if(F) {
    sampleNames = c("1-RV", "3-ES", "9-LA")
    i=1
    for (sampleName in sampleNames) {
        probesetTxtFolder = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/premiers_E_recus/all_probeset"
        pathToProbesetTxt = paste0(probesetTxtFolder,"/",sampleName,".probeset.txt")
        print(pathToProbesetTxt)
        cgh = rCGH::readAffyOncoScan(pathToProbesetTxt, sampleName=i)
        cgh@cnSet = removePointsForQuickPlotting(cgh@cnSet, 100)
        saveRDS(cgh, paste0(working_dir, "/rCGH_", i, ".RDS"))
        i=i+1
    }
    s1 = readRDS(paste0("rCGH_", sampleName, ".RDS"))
}


pipeline_rCGH = function(sampleName, silent = FALSE) {
    # sampleName="3-ES"
    # sampleName="2-AD"
    # sampleName="1-RV"
    sampleName="11-BG"
    # sampleName="9-LA" # has LOH
    ## ----readFiles----------------------------------
    probesetTxtFolder = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/premiers_E_recus/all_probeset"
    pathToProbesetTxt = paste0(probesetTxtFolder,"/",sampleName,".probeset.txt")
    print(pathToProbesetTxt)
    before = Sys.time()
    cgh = rCGH::readAffyOncoScan(pathToProbesetTxt, sampleName=sampleName)
    ## remove sex chromosomes data
    cgh@cnSet = dplyr::filter(cgh@cnSet, ChrNum<23)
    ##-- create a column "absolute position" for better plots
    cgh@cnSet = getAbspos_probeset(cgh@cnSet)
    # Organize cghSet columns in the same layout as test rCGH file
    cgh@cnSet = cgh@cnSet %>% mutate(SmoothSignal=rep(NA, length(cgh@cnSet[,1])), .before=Allele.Difference)
    colnames(cgh@cnSet) = c("ProbeName","ChrNum","ChrStart","Log2Ratio","WeightedLog2Ratio","SmoothSignal","Allele.Difference","NormalDiploid","BAF","absPos")
    cgh@cnSet = cgh@cnSet %>% relocate(Allele.Difference, .after = NormalDiploid)

    ## ----adjustSignal-------------------------------------------------------------
    ## removing probes with Log2Ratio=NaN 
    LRRData = cgh@cnSet
    cgh@cnSet = dplyr::filter(LRRData, !is.na(LRRData["Log2Ratio"]))
    if (silent) {
        cghAdj <- hush(rCGH::adjustSignal(cgh, nCores=1, suppOutliers=T, verbose=F, Scale=T))
    } else {
        cghAdj <- rCGH::adjustSignal(cgh, nCores=1, suppOutliers=T, verbose=F, Scale=T)
    }

    if (F) {
        ### see difference between and after adjusting data
        adj = removePointsForQuickPlotting(cghAdj@cnSet)
        raw = removePointsForQuickPlotting(cgh@cnSet)
        plot(y=raw$Log2Ratio, x=raw$absPos, pch=20, main=paste0(sampleName," Raw data"), cex=0.001, xlab="Genomic position (bp)", ylab="Log Ratio", ylim=c(-6,2))
        plot(y=adj$Log2Ratio, x=adj$absPos, pch=20, main=paste0(sampleName," Adjusted data"), cex=0.001, xlab="Genomic position (bp)", ylab="Log Ratio", ylim=c(-6,2))
        # plot(y=adj$adjMan, x=adj$absPos, pch=20, main=paste0(sampleName," adjusted data using scale()"), cex=0.001, xlab="Genomic position (bp)", ylab="Log Ratio", ylim=c(-6,2))
        # plot(y=adj$adjManMan, x=adj$absPos, pch=20, main=paste0(sampleName," adjusted data using manual sd calculation then dividing by it"), cex=0.001, xlab="Genomic position (bp)", ylab="Log Ratio", ylim=c(-6,2))
    }


    ## ----SegmentCGH---------------------------------------------------------------
    cghSeg <- rCGH::segmentCGH(cghAdj, Smooth=TRUE, nCores=1, minLen=10, verbose=TRUE)
    
    # source(file.path(working_dir, "rCGH_dev.R"))
    # cghSeg <- segmentCGH_custom(cghAdj, Smooth=TRUE, nCores=1, minLen=10, verbose=TRUE)
    

     ## ----segTable-----------------------------------------------------------------
    # head(segTable_rCGH)

    ## ----EMnormalize--------------------------------------------------------------
    cghNorm <- rCGH::EMnormalize(cghSeg)
    after = Sys.time()
    setInfo(cghNorm, "runTime") = difftime(after, before, units="secs")
    # Retrieve this using currCallRes@info[["runTime"]] or  getInfo(currCallRes, "runTime")
    
    if(F) {
        ### plot to compare with seg before this 2nd normalisation
        segDfNorm = cghNorm@cnSet
        segDfNorm = removePointsForQuickPlotting(segDfNorm)
        plot(y=segDfNorm$Segm , x=segDfNorm$absPos, pch=20, main=paste0(sampleName," Normalized segmented data"), cex=0.01, xlab="Genomic position (bp)", ylab="Log Ratio", ylim=c(-6,2))
    }

    ## ----plotDensity, fig.width=7, fig.height=5, fig.show='hide'------------------
    # plotDensity(cghNorm)

    
    if(F) {
        ## plot gene positions on genome
        segTable_rCGH <- getSegTable(cghSeg)
        geneTable <- byGeneTable(segTable_rCGH)
        head(geneTable, n=3)
        colTransitoire = colnames(geneTable)
        index_chr = which(colTransitoire=="chr")
        index_start = which(colTransitoire=="chrStart")
        index_end = which(colTransitoire=="chrEnd")
        colTransitoire[index_chr] = "chrom"
        colTransitoire[index_start] = "loc.start"
        colTransitoire[index_end] = "loc.end"
        colnames(geneTable) = colTransitoire
        generateGrid("gene table", mode="LRR")
        apply(geneTable, 1, plotSeg, "Log2Ratio")
        }

    if (F) {
        ############ possibilities of the package
        ## ----byGeneTable2-------------------------------------------------------------
        byGeneTable(segTable_rCGH, "erbb2", genome = "hg19")[,1:6]
        byGeneTable(segTable_rCGH, "erbb2", genome = "hg18")[,1:6]
    
        ## ----getParams----------------------------------------------------------------
        paramsAfterRecenter = getParam(cghNorm)
        # getParam(cghNorm)[1:3] ## K params
    
        ## ----getProfile, fig.width=7.7, fig.height=9.5, fig.show='hide'---------------
        multiplot(cghNorm, symbol = c("egfr", "erbb2"))
    
        ## ----recenter, fig.width=7.5, fig.height=4, fig.show='hide'-------------------
        # Recentering on peak #2
        recenter(cghNorm) <- 3 ## this does not change the parameters, although the R help on this function says different
        
        plotProfile(cghNorm, symbol = c("egfr", "erbb2"))
    
        ### also plotting LOH
        plotLOH(cghNorm)
    
        ## ----view, eval=FALSE, echo=TRUE----------------------------------------------
        view(cghNorm)
    
        ## ----exampleFiles-------------------------------------------------------------
        list.files(system.file("extdata", package = "rCGH"))
    
        ## ----session------------------------------------------------------------------
        sessionInfo()
    }
    
    return(cghNorm)
}



main = function() {
    source(file.path(working_dir, "oncoscanR_functions.R"))
    source(file.path(working_dir, "rCGH_functions.R"))
    source(file.path(working_dir, "CGHcall_functions.R"))
    ## define paths
    # sampleNames = c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT", "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )
    sampleNames = c("2-AD")
    ## initialize list to contain all rCGH objects, and list of segTables
    rCGHResults = list()
    segTables = list()
    ## initialize df to contain all GIs
    GIdf_rCGH = data.frame(matrix(ncol=4, nrow=length(sampleNames)))
    colnames(GIdf_rCGH) = c("GI", "nbAlter", "nbChr", "runTime")
    rownames(GIdf_rCGH) = sampleNames
    rCGHdir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/rCGH/"
    
    segsType = "raw"
    if(!dir.exists(rCGHdir)) dir.create(rCGHdir)
    for (s in 1:length(sampleNames)) {
        currSampleName = sampleNames[s]
        print(paste0("processing pipeline for sample ", currSampleName))
        currCallRes = pipeline_rCGH(currSampleName, silent=FALSE)
        rCGHResults = append(rCGHResults,list(currCallRes))
        # extract seg table
        if(segsType=="raw") {
            segTable_rCGH = currCallRes@cnSet[["Segm"]]
        } else if(segsType=="CN") {
            segTable_rCGH = getSegTable(currCallRes)
            # extract GI result
            GI_res = calcGI_rCGH(segTable_rCGH)
            cleanRunTime = round(as.numeric(getInfo(currCallRes, "runTime")), 2)
            GIdf_rCGH[s,] = c(GI_res, cleanRunTime)
        }
        segTables = append(segTables,list(segTable_rCGH))
    }
    
    if(segsType == "raw") {
        source(file.path(working_dir, "CGHcall_functions.R"))
        segTables_rCGH = data.frame(segTables[[1]])
        for (sample in 1:length(segTables)) {
            segTables_rCGH[sample] = segTables[[sample]]
            colnames(segTables_rCGH) = sampleNames
        }
        rowsInfo = currCallRes@cnSet[c("ChrNum", "ChrStart")]
        rowsInfo$ChrEnd = rowsInfo$ChrStart+20
        segTables_rCGH = cbind(rowsInfo, segTables_rCGH)
        colnames(segTables_rCGH) = c("Chromosome", "Start", "End", colnames(segTables_rCGH)[4:length(segTables_rCGH)])
        segTables_rCGH = getSegTables(segTables_rCGH, sampleNames)
        segTables_rCGH
    }
    res_rCGH = list(GIdf_rCGH, segTables, rCGHResults)
    
    return(res_rCGH)
}


# runs only when script is main
if (sys.nframe() == 0){
    res_rCGH = main()
    GI_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/GI_all_methods"
    GIdf_rCGH = res_rCGH[[1]]
    segTables = res_rCGH[[2]]
    ############## save GI data to file
    source(file.path(working_dir, "crossPackagesFunctions.R"))
    saveGI_ResToFile(GIdf_rCGH, "rCGH")
    ## plot seg tables
    ## -- plotting and saving estimated CN
    # png(paste0(outputFolder, "/segsUsedForGI.png"), width=1500, height=600)
    currSampleName = "14-CJ"
    segTable_rCGH = segTables[[13]]
    generateGrid(paste0(currSampleName," estimated copy number"), mode="CN")
    apply(segTable_rCGH, 1, plotSeg_rCGH, "probes.Sd")
    # dev.off()

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
saveSegTables(segTables_rCGH, rCGHdir)

}
















############ pour slides soutenance
if (FALSE) {
    soutenance_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/res_soutenance"


    ylim = c(-3,3)
    PercentRowsToRemove = 30
    xlab = "position genomique"
    pkgName="rCGH"
    ###### plot dimensions:
    pngWidth = 350
    pngHeight=300

    ### plot raw value per probe
    state="raw"
    rawVals = cgh@cnSet
    rawVals = removePointsForQuickPlotting(rawVals, PercentRowsToRemove)
    savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
    png(savePath, width = pngWidth, height = pngHeight)
    plot(rawVals$Log2Ratio, pch = 20, cex = 0.01, ylab = 'log Ratio', xlab = xlab, ylim = ylim)
    dev.off()

    ### plot adj value per probe
    state = "adj"
    adjVals = cghAdj@cnSet
    adjVals = removePointsForQuickPlotting(adjVals, PercentRowsToRemove)
    savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
    png(savePath, width = pngWidth, height = pngHeight)
    plot(adjVals$Log2Ratio, pch = 20, cex = 0.1, ylab = 'log Ratio', xlab = xlab, ylim = ylim)
    dev.off()

    ### plot seg value per probe
    state = "seg"
    segVals = as.data.frame(cghSeg@cnSet)
    segVals = removePointsForQuickPlotting(segVals, PercentRowsToRemove)
    savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
    png(savePath, width = pngWidth, height = pngHeight)
    plot(y=segVals$Segm, x=(1:length(segVals[[1]])), pch = 20, cex = 0.1, ylab = 'log Ratio', xlab = xlab, ylim = ylim)
    dev.off()

    ### plot Call value per probe
    state = "norm"
    savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
    vals = cghNorm@cnSet
    vals = removePointsForQuickPlotting(vals, PercentRowsToRemove)
    CNvals = vals[["estimCopy"]]
    genomicPos = as.numeric(1:length(CNvals))
    png(savePath, width = pngWidth, height = pngHeight)
    plot(y = CNvals, x=genomicPos, pch=20, cex = 0.1, ylab = 'Nombre de copies', xlab = xlab)
    dev.off()

}





