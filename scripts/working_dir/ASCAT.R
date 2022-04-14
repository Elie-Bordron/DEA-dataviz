## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)
#loading libraries
library(ASCAT)
library(dplyr)



pathToOSCHP = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_OSCHP/5-LD.OSCHP"
pathToCel = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_CEL"
pathToATCelFile = file.path(pathToCel, "5-LD_AT_(OncoScan_CNV).CEL")
pathToGCCelFile = file.path(pathToCel, "5-LD_GC_(OncoScan_CNV).CEL")
outputFolder = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/ASCAT/5-LD"
outputFolder_plots = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/ASCAT/5-LD/plots/"
sampleName = "5-LD"

# to load functions from EaCoN_functions.R
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/EaCoN_functions.R")
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/ASCAT_functions.R")
# loading data
rawData = custom_OS.Process(ATChannelCel = pathToATCelFile, GCChannelCel = pathToGCCelFile, samplename = sampleName, oschp_file=pathToOSCHP, force=T, oschp.keep=T, return.data=TRUE, plot=F, write.data=F)
# segmenting data using ascat.aspcf(ASCAT_obj)
segData =  ASCAT::ascat.aspcf(ASCATobj = rawData$data, ascat.gg = rawData$germline, penalty = 50) # 50 is EaCoN default value
# estimating copy number & ploidy & cellularity using ASCAT::ascat.runAscat 
if(!dir.exists(outputFolder)) dir.create(outputFolder)
gamma = 0.35 # gammaRange = seq(0.35, 0.95, 0.05)
callData = ASCAT::ascat.runAscat(ASCATobj = segData, gamma = gamma, img.dir=outputFolder)

### getting metrics
# run this line before running ascat.metrics() to avoid error "object 'ascat.bc' not found"
ascat.bc=segData
metrics = ascat.metrics(segData,callData)

### getting plots
if(!dir.exists(outputFolder_plots)) dir.create(outputFolder_plots)
callDataToPlot = callData
segDataToPlot = segData

# remove "chr" so segments$chr matches hg38 notation (and is easier to see on plots)
callDataToPlot$segments$chr = lapply(callDataToPlot$segments$chr,stringr::str_replace, "chr", "")
segDataToPlot$chrs = as.vector(lapply(segDataToPlot$chrs,stringr::str_replace, "chr", ""))
# remove Y chromosome data
callDataToPlot$segments = dplyr::filter(callDataToPlot$segments,chr!="Y")
ascat.plotAdjustedAscatProfile(callDataToPlot,REF="hg38",png_prefix=outputFolder_plots)
## plot logR and BAF raw Data
ascat.plotRawData(segDataToPlot,img.dir=outputFolder_plots) #5-LD.tumour.png
## plot logR and BAF data before vs after segmentation
ascat.plotSegmentedData(segDataToPlot,img.dir=outputFolder_plots) #5-LD.ASPCF.png
ascat.plotSunrise(callData)
## later
ascat.plotAscatProfile(callData)
ascat.plotGenotypes(callData)
ascat.plotNonRounded(callData)


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

