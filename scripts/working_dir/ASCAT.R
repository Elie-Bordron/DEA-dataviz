## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)
#loading libraries
library(ASCAT)
library(dplyr)


sampleName = "7-DG"
pathToOSCHP = paste0("C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_OSCHP/",sampleName,".OSCHP")
pathToCel = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_CEL"
pathToATCelFile = file.path(pathToCel, paste0(sampleName,"_AT_(OncoScan_CNV).CEL"))
pathToGCCelFile = file.path(pathToCel, paste0(sampleName,"_GC_(OncoScan_CNV).CEL"))
outputFolder = paste0("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/ASCAT/",sampleName)
outputFolder_plots = file.path(outputFolder, "plots")

# to load functions from EaCoN_functions.R
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/EaCoN_functions.R")
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/ASCAT_functions.R")
# loading data
rawData = custom_OS.Process(ATChannelCel = pathToATCelFile, GCChannelCel = pathToGCCelFile, samplename = sampleName, oschp_file=pathToOSCHP, force=T, oschp.keep=T, return.data=TRUE, plot=F, write.data=F)
# segmenting data using ascat.aspcf(ASCAT_obj)
segData =  ASCAT::ascat.aspcf(ASCATobj = rawData$data, ascat.gg = rawData$germline, penalty = 50) # 50 is EaCoN default value
# estimating copy number & ploidy & cellularity using ASCAT::ascat.runAscat. Also generates rawprofile, ascatprofile and sunrise plots
if(!dir.exists(outputFolder)) dir.create(outputFolder)
gamma = 0.35 # gammaRange = seq(0.35, 0.95, 0.05)

# callData = ASCAT::ascat.runAscat(ASCATobj = segData, gamma = gamma, img.dir=outputFolder) 
# custom psi and rho
psi=1.48
rho=1
callData = ASCAT::ascat.runAscat(ASCATobj = segData, gamma = gamma, img.dir=outputFolder, img.prefix=paste0("custom_psi",psi, "_rho", rho), psi_manual=psi, rho_manual = rho)


### getting metrics
# run this line before running ascat.metrics() to avoid error "object 'ascat.bc' not found"
if(F){
ascat.bc=segData
metrics = ascat.metrics(segData,callData)
}

### getting plots
if(!dir.exists(outputFolder_plots)) dir.create(outputFolder_plots)
callDataToPlot = callData
segDataToPlot = segData

# remove "chr" so segments$chr matches hg38 notation (and is easier to see on plots)
callDataToPlot$segments$chr = lapply(callDataToPlot$segments$chr,stringr::str_replace, "chr", "")
callDataToPlot$segments_raw$chr = lapply(callDataToPlot$segments_raw$chr,stringr::str_replace, "chr", "")
segDataToPlot$chrs = as.vector(lapply(segDataToPlot$chrs,stringr::str_replace, "chr", ""))
if (F){
# remove Y chromosome data
callDataToPlot$segments = dplyr::filter(callDataToPlot$segments,chr!="Y")
ascat.plotAdjustedAscatProfile(callDataToPlot,REF="hg38",png_prefix=outputFolder_plots)
## plot logR and BAF raw Data
ascat.plotRawData(segDataToPlot,img.dir=outputFolder_plots) #5-LD.tumour.png
## plot logR and BAF data before vs after segmentation
ascat.plotSegmentedData(segDataToPlot,img.dir=outputFolder_plots) #5-LD.ASPCF.png
## later
ascat.plotAscatProfile(callData)
ascat.plotGenotypes(callData)
ascat.plotNonRounded(callData)
## automatically generated already
ascat.plotSunrise(callData$distance_matrix)
}

## plot raw values for segments of allele A and B
plotSeg = function(seg_df, allChrs, lengthOfChrs, drawPolygons=F){
    # print(c("seg_df: ", seg_df))
    print(c("seg_df[2]: ", seg_df[2]))
    #get current chromosome
    currChr = which(allChrs==seg_df[2])
    #get endpos of last segment of previous chromosome
    if (currChr!=1){
        pos0CurrChr = sum(lengthOfChrs[1:currChr-1])
    } else {
        pos0CurrChr=0
    }
    # drawing a segment on plot for each segment of the genome
    segStartPos = as.numeric(seg_df[3])
    segEndPos = as.numeric(seg_df[4])
    logRAlleleA = as.numeric(seg_df[7])
    logRAlleleB = as.numeric(seg_df[8])
    segments(pos0CurrChr+segStartPos, logRAlleleA, pos0CurrChr+segEndPos, logRAlleleA, col="dark blue", lwd=2)
    # segments(pos0CurrChr+segStartPos, logRAlleleB, pos0CurrChr+segEndPos, logRAlleleB, col="green", lwd=2)
    # drawing polygone joining this segment to closest non-null integer
    if(drawPolygons){
        polygon(x=c(pos0CurrChr+segStartPos, pos0CurrChr+segEndPos, pos0CurrChr+segEndPos, pos0CurrChr+segStartPos), y=c(logRAlleleA,logRAlleleA,round(logRAlleleA),round(logRAlleleA)), col = "dark blue", density = 10)
        # polygon(x=c(pos0CurrChr+segStartPos, pos0CurrChr+segEndPos, pos0CurrChr+segEndPos, pos0CurrChr+segStartPos), y=c(logRAlleleB,logRAlleleB,round(logRAlleleB),round(logRAlleleB)), col = "green", density = 10)
    }
}

generateGrid = function(graph_title) {
    #create empty plot to add things in
    plot(1, ylim=c(-1,4), xlim=c(0,3*10**9),col="white", xaxt="n", yaxt="n", xlab="nombre de bases", ylab="log ratio", main=graph_title)
    #  add horizontal grid
    for (CN in c(-3:8)) {
    abline(h=CN, col="dark grey")
    }
    # X-axis
    axis(1, at = c(0, 5*10**8, 1*10**9, 1.5*10**9, 2*10**9, 2.5*10**9, 3*10**9))
    # Y-axis
    axis(2, at = c(-3:8))
}


# remove small segments that pollute visualization
segDf = callDataToPlot$segments_raw
segDf_trimmed = dplyr::filter(segDf, endpos-startpos>10**7)
allChrs = unique(as.vector(segDf_trimmed[2]))
lengthOfChrs = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
# draw seg data
graph_title = paste0(sampleName,"      ploidie:",psi,", cellularit�:",rho)
# toSave = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/plot.png"
path_plot_no_polygon = paste0(outputFolder_plots, graph_title, ".png")
# png(path_plot_no_polygon)
generateGrid(graph_title)
apply(segDf_trimmed, 1, plotSeg, allChrs, lengthOfChrs)
# dev.off()

# png(paste0(sampleName,"ploidie:",psi,",cellularit�:",rho, "_polygon.png"))
generateGrid(graph_title)
apply(segDf_trimmed, 1, plotSeg, allChrs, lengthOfChrs, drawPolygons=T)
# dev.off()



### generating ASPCF segmentation visuals
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

