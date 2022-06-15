

## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
ascatFolder = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/ASCAT"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)
options("max.print" = 100)
#loading libraries
library(ASCAT)
library(dplyr)

source(file.path(working_dir, "EaCoN_functions.R"))
source(file.path(working_dir, "ASCAT_functions.R"))
source(file.path(working_dir, "oncoscanR_functions.R"))
source(file.path(working_dir, "crossPackagesFunctions.R"))
source(file.path(working_dir, "rCGH_functions.R"))


#### functions


# sampleName = "5-LD"
# sampleName = "6-VJ"
# sampleName = "7-DG"
# sampleName = "17-VV"
sampleNames = c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )



#### Save ASCAT objects from .RDS files
if(F) {
    sampleNames = c("1-RV", "3-ES", "9-LA")
    i=1
    for (sampleName in sampleNames) {
        OSCHPFolder = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/premiers_E_recus/all_OSCHP"
        OSCHPPattern = paste0(".*",sampleName,".*","\\.OSCHP")
        nameOSCHP = list.files(OSCHPFolder, pattern=OSCHPPattern)
        pathToOSCHP = file.path(OSCHPFolder,nameOSCHP)
        pathToCel = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/premiers_E_recus/all_CEL"
        pathToATCelFile = file.path(pathToCel, paste0(sampleName,"_AT_(OncoScan_CNV).CEL"))
        pathToGCCelFile = file.path(pathToCel, paste0(sampleName,"_GC_(OncoScan_CNV).CEL"))
        # loading data
        rawData = custom_OS.Process(ATChannelCel = pathToATCelFile, GCChannelCel = pathToGCCelFile, samplename = sampleName, oschp_file=pathToOSCHP, force=T, oschp.keep=T, return.data=TRUE, plot=F, write.data=F) # To do as less data processing as possible, set: wave.renorm = F, gc.renorm=F,
        rawData = removeSexChrData(rawData)
        # cgh@cnSet = removePointsForQuickPlotting(cgh@cnSet, 100
        saveRDS(rawData, paste0(working_dir, "/ASCAT_", i, ".RDS"))
        i=i+1
    }
}


pipelineASCAT = function(sampleName,outputFolder, ascatFolder, penalty=50) {
    ## set paths
    OSCHPFolder = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/premiers_E_recus/all_OSCHP"
    # pathToOSCHP = paste0(OSCHPFolder,sampleName,".OSCHP")
    OSCHPPattern = paste0(".*",sampleName,".*","\\.OSCHP")
    nameOSCHP = list.files(OSCHPFolder, pattern=OSCHPPattern)
    pathToOSCHP = file.path(OSCHPFolder,nameOSCHP)
    pathToCel = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/premiers_E_recus/all_CEL"
    pathToATCelFile = file.path(pathToCel, paste0(sampleName,"_AT_(OncoScan_CNV).CEL"))
    pathToGCCelFile = file.path(pathToCel, paste0(sampleName,"_GC_(OncoScan_CNV).CEL"))
    # outputFolder_plots = file.path(outputFolder, "plots")
    
    # loading data
    rawData = custom_OS.Process(ATChannelCel = pathToATCelFile, GCChannelCel = pathToGCCelFile, samplename = sampleName, oschp_file=pathToOSCHP, force=T, oschp.keep=T, return.data=TRUE, plot=F, write.data=F) # To do as less data processing as possible, set: wave.renorm = F, gc.renorm=F,     
    rawData = removeSexChrData(rawData)
    # segmenting data using ascat.aspcf(ASCAT_obj). This step generates .txt segmentation files
    if(!dir.exists(outputFolder)) dir.create(outputFolder)
    source(file.path(working_dir,"ASCAT_functions.R"))
    segData =  ASCAT::ascat.aspcf(ASCATobj = rawData$data, ascat.gg = rawData$germline, penalty = penalty, out.dir=outputFolder) # 50 is EaCoN default value
    # segData =  custom_ascat.aspcf(ASCATobj = rawData$data, ascat.gg = rawData$germline, penalty = 50, out.dir=outputFolder) # 50 is EaCoN default value
    # estimating copy number & ploidy & cellularity using ASCAT::ascat.runAscat. Also generates rawprofile, ascatprofile and sunrise plots

    ## if sample has not a best gamma value known in the reference file, process it with all gammas and update the file.
    gammaRange = seq(0.35, 0.95, 0.05)
    bestGammas = read.table(file.path(ascatFolder,"bestGammas.txt"), h=T, sep='\t')
    # Retaining most trustworthy results across different gamma values
    if(any(bestGammas$sample==sampleName)) {
        print("using previously defined gamma value")
        bestGamma = dplyr::filter(bestGammas, sample==sampleName)$bestGamma
        callData = ASCAT::ascat.runAscat(ASCATobj=segData,gamma=bestGamma,img.dir=outputFolder,img.prefix=paste0("run_bestGamma=", bestGamma, "_"))
        GOF_all_gammas=NULL
    } else {
        outputFolderGammaRange = file.path(outputFolder, "gammaRange")
        if(!dir.exists(outputFolderGammaRange)) dir.create(outputFolderGammaRange)
        callData = NULL
        GOF_all_gammas = c()
        print("searching for best gamma...")
        for (gamma in gammaRange){
            currCallData = ASCAT::ascat.runAscat(ASCATobj=segData,gamma=gamma,img.dir=outputFolderGammaRange,img.prefix=paste0("normal_run", "_gamma=", gamma, "_"))
            if(is.null(currCallData$goodnessOfFit)) {
                print("currCallData$goodnessOfFit is NULL, adding 0.0 as goodness value")
                GOF_all_gammas = c(GOF_all_gammas, 0.0)
            } else {
                GOF_all_gammas = c(GOF_all_gammas, currCallData$goodnessOfFit)
            }
            # print(c("currCallData's goodness of fit: ': ", currCallData$goodnessOfFit))
            if(!is.null(callData)) {
                # print("callData is NOT null")
                if(!is.null(currCallData$goodnessOfFit)) {
                    # print("currCallData has a goodnessOfFit")
                    if(callData$goodnessOfFit<currCallData$goodnessOfFit) {
                        # print("currCallData has a better goodnessOfFit than callData, we replace the latter by the former")
                        callData = currCallData
                        callData = append(callData, gamma)
                    }
                }
            } else {
                # print("calldata is NULL")
                if(!is.null(currCallData$goodnessOfFit)) {
                    # print("currCallData has a goodnessOfFit, it becomes the first callData")
                        callData = currCallData
                        callData = append(callData, gamma)
                } 
            }
        }
        ## to view best solution: 
        GOF_all_gammas = as.data.frame(GOF_all_gammas)
        GOF_all_gammas$gamma = gammaRange; colnames(GOF_all_gammas)=c("goodness_of_fit", "gamma")
        bestSolution = dplyr::filter(GOF_all_gammas, goodness_of_fit==max(goodness_of_fit))
        print(c("best solution: ", bestSolution))
    }
    return(list(callData, GOF_all_gammas, segData))
}




################### calculate GI





#### main
ascatFolder = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/ASCAT"
sampleNames = c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )
sampleName = "12-BC"
# sampleNames = c("5-LD", "6-VJ", "8-MM")
# GI_ASCAT_df = data.frame(rep(NA, length(sampleNames)))
GI_ASCAT_df = data.frame(matrix(ncol = 5, nrow = length(sampleNames)))
colnames(GI_ASCAT_df) = c("GI", "nbAlter", "nbChr", "runTime", "estimPurity")
rownames(GI_ASCAT_df) = sampleNames
## load GI functions
source(file.path(working_dir, "oncoscanR_functions.R"))
################# start of loop
s=1
penalty = 50
segTables = list()
segsType = "raw"
for (s in 1:length(sampleNames)) {
    sampleName = sampleNames[s]
    before = Sys.time()
    outputFolder = file.path(ascatFolder,sampleName)
    resPipeline = pipelineASCAT(sampleName,outputFolder,ascatFolder, penalty)
    callData = resPipeline[[1]]
    GammasTested = resPipeline[[2]]
    if(!is.null(GammasTested)) {
        message("Writing gamma values along with goodness of fit in a table...")
        write.table(GammasTested, file.path(outputFolder, "gammas_tested.txt"))
    }
    segTable_clean = cleanASCATSegData(callData, trimData=F)
    GI_res = calcGI_ASCAT(segTable_clean)
    
    ## keep this GI in a df along with its sampleName
    after=Sys.time()
    # print(c("GI_res: ", GI_res))
    GI_res = append(GI_res,round(as.numeric(difftime(after, before, units = "secs")), 2)) ## add run time
    # print(c("GI_res with runTime: ", GI_res))
    GI_res = append(GI_res,callData$purity) ## add purity estimate
    print(c("GI_res with runTime and purity: ", GI_res))
    GI_ASCAT_df[sampleName, ] = GI_res
    if(segsType=="raw") {
        segTable = resPipeline[[3]][["Tumor_LogR_segmented"]]
    } else if (segsType=="CN") {
        segTable = callData$segments
    }
    segTables = append(segTables, list(segTable))
}
GI_ASCAT_df_copy = GI_ASCAT_df
################# end of loop

# to use: 
# segTables
# GI_ASCAT_df


############## save GI data to file
source(file.path(working_dir, "crossPackagesFunctions.R"))
saveGI_ResToFile(GI_ASCAT_df_copy, "ASCAT", "estimPurity")

if(F) {
    ### Organizing rows in the same order as oncoscan GI results
    rowsLayout = GI_OncoscanR$sample
    GI_ASCAT_df_copy$sample = rownames(GI_ASCAT_df_copy)
    GI_ASCAT_df_copy = GI_ASCAT_df_copy[rowsLayout,]
    
    ## fusing tables 
    GI_oncosAndAscat = as.data.frame(cbind(GI_ASCAT_df_copy$sample, GI_ASCAT_df_copy$GI_ASCAT, GI_OncoscanR$GI_oncoscanR))
    colnames(GI_oncosAndAscat) = c("sample", "GI_ASCAT", "GI_OncoscanR")
    rownames(GI_oncosAndAscat) = GI_oncosAndAscat$sample

    ## organizing rows in reading order
    rowsLayout = c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )
    GI_oncosAndAscat = GI_oncosAndAscat[rowsLayout,]
    barplot(as.numeric(GI_oncosAndAscat$GI_ASCAT), names.arg = GI_oncosAndAscat$sample)
    ## writing all GIs in a text file
    resultsDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/GI_all_methods"
    write.table(GI_oncosAndAscat,file.path(resultsDir, "gi_results_all_methods.txt"),sep="\t",row.names=FALSE)
}
############## save segments table to file

if(segsType == "raw") {
    source(file.path(working_dir, "CGHcall_functions.R"))
    # segTables_ASCAT = data.frame(segTables[[1]])
    # for (sample in 1:length(segTables)) {
    #     segTables_ASCAT[sample] = segTables[[sample]]
    #     colnames(segTables_ASCAT) = sampleNames
    # }
    segTables_ASCAT = as.data.frame(segTables)
    rowsInfo = segData[["SNPpos"]]
    rowsInfo$ChrEnd = rowsInfo$pos+20
    segTables_ASCAT = cbind(rowsInfo, segTables_ASCAT)
    colnames(segTables_ASCAT) = c("Chromosome", "Start", "End", sampleNames)
    segTables_ASCAT = getSegTables(segTables_ASCAT, sampleNames)
    segTables_ASCAT
}


saveSegTables = function(segTables, outputDir, sampleNames=c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )) {
    segTablesDir = file.path(outputDir, "segTables")
    if(!dir.exists(segTablesDir))dir.create(segTablesDir)
    for (i in 1:length(segTables)) {
        currSegTable = segTables[i]
        currFilePath = file.path(segTablesDir, paste0(sampleNames[i], ".tsv"))
        write.table(currSegTable, currFilePath, sep="\t", row.names=FALSE, quote=F)
    }
}
saveSegTables(segTables_ASCAT, ascatFolder)

##
print("end")


























############ pour slides soutenance

soutenance_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/res_soutenance"

if (TRUE){
    callDataToPlot = callData
    # remove "chr" so segments$chr matches hg38 notation (and is easier to see on plots)
    callDataToPlot$segments$chr = lapply(callDataToPlot$segments$chr,stringr::str_replace, "chr", "")
    callDataToPlot$segments_raw$chr = lapply(callDataToPlot$segments_raw$chr,stringr::str_replace, "chr", "")
    segDataToPlot$chrs = as.vector(lapply(segDataToPlot$chrs,stringr::str_replace, "chr", ""))
}

ylim = c(-1,1)
PercentRowsToRemove = 30
xlab = "position genomique"
pkgName="ASCAT"
pngWidth = 350
pngHeight=300

###### plot dimensions: w=556, h=513
### plot raw value per probe
state="raw"
rawVals = rawData[["data"]][["Tumor_LogR"]]
rawVals = removePointsForQuickPlotting(rawVals, PercentRowsToRemove)
savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
png(savePath, width = pngWidth, height = pngHeight)
plot(rawVals[,1], pch = 20, cex = 0.01, ylab = 'log Ratio', xlab = xlab, ylim = ylim)
dev.off()

### plot seg value per probe
state = "seg"
segVals = as.data.frame(segData[["Tumor_LogR_segmented"]])
segVals = removePointsForQuickPlotting(segVals, PercentRowsToRemove)
savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
png(savePath, width = pngWidth, height = pngHeight)
plot(segVals[,1], pch = 20, cex = 0.1, ylab = 'log Ratio', xlab = xlab, ylim = ylim)
dev.off()

### plot Call value per probe
state = "call"
callVals = as.data.frame(callDataToPlot$nA)
callVals$nB = callDataToPlot$nB
callVals = removePointsForQuickPlotting(callVals, PercentRowsToRemove)
callVals$totCN = callVals[1] + callVals[2]
indexes = c(1:length(callVals$totCN[[1]]))
savePath = paste0(soutenance_dir, "/", pkgName, "_", sampleName, "_", state, ".png")
png(savePath, width = pngWidth, height = pngHeight)
plot(y = callVals$totCN[[1]], pch=20, cex = 0.1, x = indexes, ylab = 'Nombre de copies', xlab = xlab)
dev.off()














### generating ASPCF segmentation visuals
if (F){
    # 2 segments in 1 plot
    vals = c(rnorm(100,2,0.1),rnorm(100,0,0.1))
    # 1 segment in 1 plot
    vals = c(rnorm(200,1,0.1))
    mu = mean(vals)
    plot(vals,pch=20, ylim=c(0.5,1.5), ylab="log ratio", xlab="genomic position")
    abline(h=mu, col="red")
    sumvals=0
    for (i in 1:length(vals)) {
        segments(i,vals[i], i,mu)
        sumvals = sumvals + (vals[i] - mu)**2
    }
    print(sumvals)
}
    




################## run ASCAT with custom psi and rho
if (F) {
    gamma = 0.45 
    psi=1.75 # 1.75
    rho=0.83 # 0.83
    callData = ASCAT::ascat.runAscat(ASCATobj = segData, gamma = gamma, img.dir=outputFolder,
        img.prefix=paste0("custom_run","_gamma=", gamma, "_psi=",psi, "_rho=", rho), psi_manual=psi, rho_manual = rho)
    print(c("goodness of fit ", ", gamma ", "ploidy ", "cellularity"))
    print(c(callData$goodnessOfFit, gamma, callData$ploidy, callData$purity))
}

#################### extract plots and results

if(F){
    ### getting metrics
    # run this line before running ascat.metrics() to avoid error "object 'ascat.bc' not found"
    ascat.bc=segData
    metrics = ascat.metrics(segData,callData)
}


### getting plots
outputFolder_plots = outputFolder
if(!dir.exists(outputFolder_plots)) dir.create(outputFolder_plots)
callDataToPlot = callData
segDataToPlot = segData




### plotting raw values of both alleles of a profile
if (F) {
    # remove small segments that pollute visualization
    segDf_trimmed = callDataToPlot$segments_raw
    segDf_trimmed = dplyr::filter(segDf_trimmed, endpos-startpos>10**7)
    allChrs = unique(as.vector(segDf_trimmed[2]))
    # draw seg data
    graph_title = paste0(sampleName,"      ploidy:",round(callDataToPlot$ploidy,2),", cellularity:",callDataToPlot$purity)
    # dir_customRun = paste0(outputFolder_plots,"/","customRun_gamma=",gamma,"_ploidy=",psi,"_cellularity=",rho,"/")
    dir_customRun = paste0(outputFolder_plots,"/","ploidy=",round(callDataToPlot$ploidy,2),"_cellularity=",callDataToPlot$purity,"/")
    if(!dir.exists(dir_customRun)) dir.create(dir_customRun)
    ### saving plots to png
    path_plot_no_polygon = paste0(dir_customRun, "nopol.png")
    png(path_plot_no_polygon, width=700, height=484)
    # png(path_plot_no_polygon)
    generateGrid(graph_title)
    apply(segDf_trimmed, 1, plotSeg_ASCAT, allChrs)
    dev.off()
    
    path_plot_polygon = paste0(dir_customRun, "polygon.png")
    # png(path_plot_polygon, width=900, height=623)
    png(path_plot_polygon, width=700, height=484)
    generateGrid(graph_title)
    apply(segDf_trimmed, 1, plotSeg_ASCAT, allChrs, drawPolygons=T)
    dev.off()
}



if (F){
    # remove Y chromosome data
    callDataToPlot$segments = dplyr::filter(callDataToPlot$segments,chr!="Y")
    ascat.plotAdjustedAscatProfile(callDataToPlot,REF="hg19",png_prefix=outputFolder_plots)
    ## plot logR and BAF raw Data
    ascat.plotRawData(segDataToPlot,img.dir=outputFolder_plots) #5-LD.tumour.png
    ## plot logR and BAF data before vs after segmentation
    ascat.plotSegmentedData(segDataToPlot,img.dir=outputFolder_plots) #5-LD.ASPCF.png
    ## later
    ascat.plotAscatProfile(segData)          # purity for 5-LD is 90%
    ascat.plotAscatProfile(callDataToPlot[["nA"]], callDataToPlot[["nB"]], ploidy=callDataToPlot[["ploidy"]], rho=0.9,          # purity for 5-LD is 90%
                           goodnessOfFit=callDataToPlot[["goodnessOfFit"]], nonaberrant=F) 
    ascat.plotGenotypes(callData)
    ascat.plotNonRounded(callData)
    ## automatically generated already
    ascat.plotSunrise(callData$distance_matrix)
}


### before and after calling
segDf = callDataToPlot$segments_raw
callDf = callDataToPlot$segments
# remove small segments that pollute visualization
# segDf = callDataToPlot$segments_raw
# segDf_trimmed = dplyr::filter(segDf, endpos-startpos>10**7)
allChrs = unique(as.vector(segDf_trimmed[2]))
lengthOfChrs = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
# draw seg data
graph_title = paste0(sampleName,"      ploidy:",psi,", cellularity:",rho)
dir_customRun = paste0(outputFolder_plots,"/","customRun_gamma=",gamma,"_ploidy=",psi,"_cellularity=",rho,"/")
if(!dir.exists(dir_customRun)) dir.create(dir_customRun)

### saving plots to png
path_plot_no_polygon = paste0(dir_customRun, "nopol.png")
## seg data

generateGrid(graph_title)
apply(segDf, 1, plotSeg_ASCAT, allChrs, lengthOfChrs)
## call data
generateGrid(graph_title)
apply(callDf, 1, plotSeg_ASCAT, allChrs, lengthOfChrs)





### picking BAF and LRR values for one SNP
bafLrr = rawData$data
bafLrr = cbind(bafLrr$Tumor_LogR, bafLrr$Tumor_BAF)
colnames(bafLrr) = c("logR", "BAF")
## S-tag021556 was picked

gamma = 0.45
ri = -0.750
bi = 0.817
rho = 0.8
psi = 4.5
if (T){
    print(paste0("--- ", rho, ", ", psi, "---"))
    term1_A = rho - 1 + 2^(ri/gamma) * (1-bi) * (2*(1-rho) + rho*psi)
    term1_B = rho - 1 + 2^(ri/gamma) * bi * (2*(1-rho) + rho*psi)
    term2 = rho 
    Ncopy_A = term1_A/term2
    print(Ncopy_A)
    Ncopy_B = term1_B/term2
    print(Ncopy_B)
}




distMX = callData[["distance_matrix"]][["5-LD"]]
ascat.plotSunrise(distMX)









