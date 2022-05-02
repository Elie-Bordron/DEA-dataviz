

## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)
#loading libraries
library(ASCAT)
library(dplyr)


sampleName = "5-LD"
# sampleName = "6-VJ"
# sampleName = "7-DG"
pathToOSCHP = paste0("C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_OSCHP/",sampleName,".OSCHP")
pathToCel = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_CEL"
pathToATCelFile = file.path(pathToCel, paste0(sampleName,"_AT_(OncoScan_CNV).CEL"))
pathToGCCelFile = file.path(pathToCel, paste0(sampleName,"_GC_(OncoScan_CNV).CEL"))
outputFolder = paste0("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/ASCAT/",sampleName)
outputFolder_plots = file.path(outputFolder, "plots")

## set variables
lengthOfChrs = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
                 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)


# to load functions from EaCoN_functions.R
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/EaCoN_functions.R")
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/ASCAT_functions.R")
# loading data
rawData = custom_OS.Process(ATChannelCel = pathToATCelFile, GCChannelCel = pathToGCCelFile, samplename = sampleName, oschp_file=pathToOSCHP, force=T, oschp.keep=T, return.data=TRUE, plot=F, write.data=F)
# segmenting data using ascat.aspcf(ASCAT_obj)
segData =  ASCAT::ascat.aspcf(ASCATobj = rawData$data, ascat.gg = rawData$germline, penalty = 50) # 50 is EaCoN default value
# estimating copy number & ploidy & cellularity using ASCAT::ascat.runAscat. Also generates rawprofile, ascatprofile and sunrise plots
if(!dir.exists(outputFolder)) dir.create(outputFolder)
# Retaining most trustworthy results across different gamma values
gammaRange = seq(0.35, 0.95, 0.05)
outputFolderGammaRange = file.path(outputFolder, "gammaRange")
if(!dir.exists(outputFolderGammaRange)) dir.create(outputFolderGammaRange)
callData = NULL
GOF_all_gammas = c()
for (gamma in gammaRange){
    currCallData = ASCAT::ascat.runAscat(ASCATobj = segData, gamma = gamma, img.dir=outputFolderGammaRange,
                                     img.prefix=paste0("normal_run", "_gamma=", gamma, "_"))
    GOF_all_gammas = c(GOF_all_gammas, currCallData$goodnessOfFit)
    if(is.null(callData)) {
        print("first loop")
        callData = currCallData
        callData = append(callData, gamma)
    } else if(currCallData$goodnessOfFit > callData$goodnessOfFit) {
        print("better solution found")
        callData = currCallData
        callData = append(callData, gamma)
    }
}
## to view best solution: 
GOF_all_gammas = as.data.frame(GOF_all_gammas)
GOF_all_gammas$gamma = gammaRange; colnames(GOF_all_gammas)=c("GOF", "gamma")
print(dplyr::filter(GOF_all_gammas, GOF==max(GOF)))



    
## functions to plot values for segments of allele A and B
plotSeg = function(seg_df, allChrs, lengthOfChrs, drawPolygons=F){
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
        if (currChr!=1){
            pos0CurrChr = sum(lengthOfChrs[1:currChr-1])
        } else {
            pos0CurrChr=0
        }
    }
    ## drawing a segment on plot for each segment of the genome
    segStartPos = as.numeric(seg_df[3])
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
### plotting raww values of both alleles of a profile
if (F) {
    # remove small segments that pollute visualization
    segDf = callDataToPlot$segments_raw
    segDf_trimmed = dplyr::filter(segDf, endpos-startpos>10**7)
    allChrs = unique(as.vector(segDf_trimmed[2]))
    # draw seg data
    graph_title = paste0(sampleName,"      ploidy:",psi,", cellularity:",rho)
    dir_customRun = paste0(outputFolder_plots,"/","customRun_gamma=",gamma,"_ploidy=",psi,"_cellularity=",rho,"/")
    if(!dir.exists(dir_customRun)) dir.create(dir_customRun)
    ### saving plots to png
    path_plot_no_polygon = paste0(dir_customRun, "nopol.png")
    png(path_plot_no_polygon, width=700, height=484)
    # png(path_plot_no_polygon)
    generateGrid(graph_title)
    apply(segDf_trimmed, 1, plotSeg, allChrs, lengthOfChrs)
    dev.off()
    
    path_plot_polygon = paste0(dir_customRun, "polygon.png")
    # png(path_plot_polygon, width=900, height=623)
    png(path_plot_polygon, width=700, height=484)
    generateGrid(graph_title)
    apply(segDf_trimmed, 1, plotSeg, allChrs, lengthOfChrs, drawPolygons=T)
    dev.off()
}


################### calculate GI
## get seg data from call result
segTable_raw = callData$segments_raw
segTable = callData$segments

## constants 
allChrs = unique(as.vector(segTable[2]))

## check seg data
graph_title = "default"
generateGrid(graph_title)
apply(segTable_raw, 1, plotSeg, allChrs, lengthOfChrs)
generateGrid(graph_title)
apply(segTable, 1, plotSeg, allChrs, lengthOfChrs)
























































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

if(!dir.exists(outputFolder_plots)) dir.create(outputFolder_plots)
callDataToPlot = callData
segDataToPlot = segData

if (TRUE){
# remove "chr" so segments$chr matches hg38 notation (and is easier to see on plots)
callDataToPlot$segments$chr = lapply(callDataToPlot$segments$chr,stringr::str_replace, "chr", "")
callDataToPlot$segments_raw$chr = lapply(callDataToPlot$segments_raw$chr,stringr::str_replace, "chr", "")
segDataToPlot$chrs = as.vector(lapply(segDataToPlot$chrs,stringr::str_replace, "chr", ""))
}

if (F){
# remove Y chromosome data
callDataToPlot$segments = dplyr::filter(callDataToPlot$segments,chr!="Y")
ascat.plotAdjustedAscatProfile(callDataToPlot,REF="hg38",png_prefix=outputFolder_plots)
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
