

## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)
#loading libraries
library(ASCAT)
library(dplyr)


#### functions
## functions to plot values for segments of allele A and B
plotSeg = function(seg_df, allChrs, drawPolygons=F, drawDistances=F){
    lengthOfChrs = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
    print(lengthOfChrs)
    ## get nb cols
    nbcols = length(seg_df)
    ## get current chromosome
    currChr = which(allChrs==seg_df["chr"])
    # print(c("currChr: ", currChr))
    # print(c("class(currChr): ", class(currChr)))
    ## get endpos of last segment of previous chromosome
    if((seg_df["chr"]=="X" ) || (seg_df["chr"]=="Y")){
        pos0CurrChr = sum(lengthOfChrs[1:22])
    }
    else {
        print(c("currChr in ASCAT.R: ", currChr))
        if (currChr!=1){
            pos0CurrChr = sum(lengthOfChrs[1:currChr-1])
        } else {
            pos0CurrChr=0
        }
    }
    ## drawing a segment on plot for each segment of the genome
    segStartPos = as.numeric(seg_df[3])
    segStartPos = segStartPos + pos0CurrChr
    segEndPos = as.numeric(seg_df[4])
    # CN values are on last 2 columns
    A_allele_vals = nbcols-1
    deviation = 0.01
    logRAlleleA = as.numeric(seg_df[[nbcols-1]])
    logRAlleleB = as.numeric(seg_df[[nbcols]])
    logRAlleleA_dev = logRAlleleA + deviation
    logRAlleleB_dev = logRAlleleB - deviation
    # print(c("pos0CurrChr, segStartPos: ", pos0CurrChr+segStartPos))
    segments(pos0CurrChr+segStartPos, logRAlleleA_dev, pos0CurrChr+segEndPos, logRAlleleA_dev, col="dark blue", lwd=2)
    segments(pos0CurrChr+segStartPos, logRAlleleB_dev, pos0CurrChr+segEndPos, logRAlleleB_dev, col="dark red", lwd=2)
    ## drawing polygone joining this segment to closest non-null integer
    if(drawPolygons){
        polygon(x=c(pos0CurrChr+segStartPos, pos0CurrChr+segEndPos, pos0CurrChr+segEndPos, pos0CurrChr+segStartPos), 
                y=c(logRAlleleB_dev,logRAlleleB_dev,round(logRAlleleB),round(logRAlleleB)), col = "dark red", density = 10,angle=135)
        polygon(x=c(pos0CurrChr+segStartPos, pos0CurrChr+segEndPos, pos0CurrChr+segEndPos, pos0CurrChr+segStartPos),
                y=c(logRAlleleA_dev,logRAlleleA_dev,round(logRAlleleA),round(logRAlleleA)), col = "dark blue", density = 10, angle=45)
    }
    if(drawDistances) {
        xPosDist = 
        segments()
    }
}

generateGrid = function(graph_title) {
    #create empty plot to add things in
    plot(1, ylim=c(0,2), xlim=c(0,3*10**9),col="white", xaxt="n", yaxt="n", xlab="nombre de bases", ylab="nombre de copies", main=graph_title)
    #  add horizontal grid
    for (CN in c(-3:8)) {
        abline(h=CN, col="dark grey")
    }
    # X-axis
    axis(1, at = c(0, 5*10**8, 1*10**9, 1.5*10**9, 2*10**9, 2.5*10**9, 3*10**9))
    # Y-axis
    axis(2, at = c(-3:8))
}

## functions for GI calculation
getNbChrs = function(segmentsTable) {
    # print(c("segmentsTable: ", segmentsTable))
    chrs = as.vector(segmentsTable$chr)
    # print(c("chrs: ", chrs))
    nbChr = length(unique(chrs))
    # print(c("nbChr: ", nbChr))
    return(nbChr)
}

calcGI_ASCAT = function(segmentsTable) {
    # print(c("segmentsTable: ", segmentsTable))
    nbChr = getNbChrs(segmentsTable)
    # print(c("nbChr: ", nbChr))
    nbAlter = dim(segmentsTable)[1]
    # print(c("nbAlter: ", nbAlter))
    GI = calcGI(nbAlter, nbChr)
    return(list(GI,nbAlter,nbChr))
}





sampleName = "5-LD"
sampleName = "6-VJ"
# sampleName = "7-DG"
# sampleName = "17-VV"
sampleNames = c("2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )



# to load functions from EaCoN_functions.R
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/EaCoN_functions.R")
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/ASCAT_functions.R")

pipelineASCAT = function(sampleName,outputFolder, ascatFolder) {
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
    # segmenting data using ascat.aspcf(ASCAT_obj). This step generates .txt segmentation files
    if(!dir.exists(outputFolder)) dir.create(outputFolder)
    segData =  ASCAT::ascat.aspcf(ASCATobj = rawData$data, ascat.gg = rawData$germline, penalty = 50, out.dir=outputFolder) # 50 is EaCoN default value
    # estimating copy number & ploidy & cellularity using ASCAT::ascat.runAscat. Also generates rawprofile, ascatprofile and sunrise plots

    ## if sample has not a best gamma value known in the reference file, process it with all gammas and update the file.
    gammaRange = seq(0.35, 0.95, 0.05)
    bestGammas = read.table(file.path(ascatFolder,"bestGammas.txt"), h=T, sep='\t')
    # Retaining most trustworthy results across different gamma values
    if(any(bestGammas$sample==sampleName)) {
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
    return(list(callData, GOF_all_gammas))
}




################### calculate GI

cleanASCATSegData = function(callData, trimData = F) {
    ## get seg data from call result
    segTable_raw = callData$segments_raw
    segTable = callData$segments
    
    ## check seg data
    # graph_title = paste0(sampleName, " called data")
    # generateGrid(graph_title)
    # allChrs = unique(as.vector(segTable[2]))
    # apply(segTable, 1, plotSeg, allChrs, lengthOfChrs)
    
    if(trimData) {
        ## removing segments shorter than 300 Kbp
        segTable = dplyr::filter(segTable, endpos-startpos>300000)
    }
    # graph_title = paste0(sampleName, " after removing segments shorter than 300 Kbp")
    # generateGrid(graph_title)
    # apply(segTable, 1, plotSeg, allChrs, lengthOfChrs)
    
    ## removing segments with CN=1 for each allele
    segTableClean = dplyr::filter(segTable, nMajor!=1 | nMinor!=1)
    nbOrigin = dim(segTable)[1]
    nbClean = dim(segTableClean)[1]
    nbRemoved = nbOrigin-nbClean
    print(paste0(nbRemoved, " segments removed out of ", nbOrigin))
    # graph_title = paste0(sampleName, " after removing segments of copy number=2")
    # generateGrid(graph_title)
    # apply(segTable, 1, plotSeg, allChrs, lengthOfChrs)
    return(segTableClean)   
}
    



#### main
# sampleNames = c("17-VV", "10-CB", "15-GG", "20-CJ", "5-LD",  "8-MM",  "11-BG", "16-DD", "21-DC", "9-LA",  "12-BC", "2-AD",  "6-VJ",  "13-VT", "18-JA", "3-ES",  "14-CJ", "19-BF", "4-GM",  "7-DG")
sampleNames = c("2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )
sampleName = "12-BC"
# sampleNames = c("5-LD", "6-VJ", "8-MM")
# GI_ASCAT_df = data.frame(rep(NA, length(sampleNames)))
GI_ASCAT_df = data.frame(matrix(ncol = 4, nrow = length(sampleNames)))
colnames(GI_ASCAT_df) = c("GI_ASCAT", "nbAlter", "nbChr", "runTime")
rownames(GI_ASCAT_df) = sampleNames
## load GI functions
source(file.path(working_dir, "oncoscanR_functions.R"))
segTables = list()
################# start of loop
s=1
for (s in 1:length(sampleNames)) {
    sampleName = sampleNames[s]
    before=Sys.time()
    ascatFolder = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/ASCAT"
    outputFolder = file.path(ascatFolder,sampleName)
    resPipeline = pipelineASCAT(sampleName,outputFolder,ascatFolder)
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
    GI_ASCAT_df[sampleName, ] = append(GI_res,after-before)
    segTable = callData$segments
    segTables = append(segTables, list(segTable))
}
GI_ASCAT_df_copy = GI_ASCAT_df
################# end of loop

# to use: 
# segTables
# GI_ASCAT_df
GI_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/GI_all_methods"
write.table(GI_ASCAT_df_copy,file.path(GI_dir, "gi_results_ASCAT.txt"),sep="\t",row.names=FALSE, quote=F)


library("GGally")
data(iris)
ggpairs(iris[, 1:4], lower=list(continuous="smooth", params=c(colour="blue")),
        diag=list(continuous="bar", params=c(colour="blue")), 
        upper=list(params=list(corSize=6)), axisLabels='show')

## load oncoscanR GI results file
# GI_filePath =  "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/OncoscanR/gi_results_all_methods.txt"
allGi = "C:\Users\e.bordron\Desktop\CGH-scoring\M2_internship_Bergonie\results\GI_all_methods"
GI_OncoscanR = read.table(GI_filePath, h=T)

### Organizing rows in the same order as oncoscan GI results
rowsLayout = GI_OncoscanR$sample
GI_ASCAT_df_copy$sample = rownames(GI_ASCAT_df_copy)
GI_ASCAT_df_copy = GI_ASCAT_df_copy[rowsLayout,]

## fusing tables 
GI_oncosAndAscat = as.data.frame(cbind(GI_ASCAT_df_copy$sample, GI_ASCAT_df_copy$GI_ASCAT, GI_OncoscanR$GI_oncoscanR))
colnames(GI_oncosAndAscat) = c("sample", "GI_ASCAT", "GI_OncoscanR")
rownames(GI_oncosAndAscat) = GI_oncosAndAscat$sample


## organizing rows in reading order
rowsLayout = c("2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )
GI_oncosAndAscat = GI_oncosAndAscat[rowsLayout,]
barplot(as.numeric(GI_oncosAndAscat$GI_ASCAT), names.arg = GI_oncosAndAscat$sample)
## writing all GIs in a text file
resultsDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/GI_all_methods"
write.table(GI_oncosAndAscat,file.path(resultsDir, "gi_results_all_methods.txt"),sep="\t",row.names=FALSE)



##
print("end")








































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

if (TRUE){
# remove "chr" so segments$chr matches hg38 notation (and is easier to see on plots)
callDataToPlot$segments$chr = lapply(callDataToPlot$segments$chr,stringr::str_replace, "chr", "")
callDataToPlot$segments_raw$chr = lapply(callDataToPlot$segments_raw$chr,stringr::str_replace, "chr", "")
segDataToPlot$chrs = as.vector(lapply(segDataToPlot$chrs,stringr::str_replace, "chr", ""))
}



### plotting raw values of both alleles of a profile
if (F) {
    # remove small segments that pollute visualization
    segDf = callDataToPlot$segments_raw
    segDf_trimmed = dplyr::filter(segDf, endpos-startpos>10**7)
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
    apply(segDf_trimmed, 1, plotSeg, allChrs)
    dev.off()
    
    path_plot_polygon = paste0(dir_customRun, "polygon.png")
    # png(path_plot_polygon, width=900, height=623)
    png(path_plot_polygon, width=700, height=484)
    generateGrid(graph_title)
    apply(segDf_trimmed, 1, plotSeg, allChrs, drawPolygons=T)
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
apply(segDf, 1, plotSeg, allChrs, lengthOfChrs)
## call data
generateGrid(graph_title)
apply(callDf, 1, plotSeg, allChrs, lengthOfChrs)





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









