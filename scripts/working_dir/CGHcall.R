## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
setwd(working_dir)
outputFolder = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/CGHcall/plots"
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)
## to use plotSeg()
source(file.path(working_dir, "rCGH.R"))
source(file.path(working_dir, "rCGH_functions.R"))
source(file.path(working_dir, "CGHcall_functions.R"))


## import libraries
# library(ggplot2)
library(dplyr)
library(CGHcall)

pipelineCGHcall = function(osData, tumor_prop=NULL) {
    # osData = rawProbesData
    before=Sys.time() 
    # osData = s2Probes # to run on one sample
    if(is.null(tumor_prop)) {tumor_prop=1}
    # ACGH_data <- make_cghRaw(Wilting)
    ACGH_data <- make_cghRaw(osData)
    # we want to apply fewest changes possible to data, so we want to do our own preprocess if we have time
    # cghdata = removedNaNProbes = dplyr::filter(ACGH_data, !is.na(ACGH_data[5]))
    cghdata <- preprocess(ACGH_data, maxmiss=95, nchrom=22) # because we don't need sex chromosomes data for GI.
    # plot(cghdata)
    norm.cghdata <- normalize(cghdata, method="median", smoothOutliers=F)
    seg.cghdata <- segmentData(norm.cghdata, method="DNAcopy", undo.splits="sdundo",undo.SD=3, clen=10, relSDlong=5)
    # segTable = segmented(seg.cghdata)
    postseg.cghdata <- postsegnormalize(seg.cghdata)
    # plot(postseg.cghdata, ylimit=c(-2,2))
    rawCghResult <- CGHcall(postseg.cghdata,nclass=5,cellularity=tumor_prop)
    CghResult <- ExpandCGHcall(rawCghResult,postseg.cghdata,CellularityCorrectSeg=F) # use CellularityCorrectSeg=TRUE to correct using cellularity
    after = Sys.time()
    CghResult$processingTime = round(as.numeric(difftime(after, before, units = "secs"))/dim(CghResult)[2], 2) ## If more than 1 sample is processed for this run, the average value is given to each sample. This value is stored in CghResult@phenoData@data[["processingTime"]].
    return(CghResult)
}

# ## to view which probes are removed 
# x = as.data.frame(cghdata@assayData[["copynumber"]])
# x = as.data.frame(ACGH_data@assayData[["copynumber"]])
# removedNaNProbes = dplyr::filter(x, !is.na(x["5-LD"]))
# testDf = data.frame(size=c(5,9,3), weight=c(8,8,NULL))
# dplyr::filter(testDf, testDf["weight"])



#################################### calculate GI
source(file.path(working_dir, "oncoscanR_functions.R"))

## setting paths
dataDirProbesets = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/premiers_E_recus/all_probeset"
resultsDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/CGHcall"
osData = NULL
if(F) {
    sampleNames = c("5-LD", "6-VJ", "8-MM")
    ## for loading data from individual probeset.txt files
    for (sampleName in sampleNames) {
        currSamplePath = paste0(dataDirProbesets, "/", sampleName, ".probeset.txt")
        currSample = read.table(currSamplePath, sep='\t', h=T)
        if(is.null(osData)){
            ## for first sample, initialize osData with columns Probe_id, Chr, startpos and endpos
            osData = currSample[,1:3]
            osData = cbind(osData, currSample[,3]+20)
        }
        ## add log ratio values for current sample, including the first one
        osData = cbind(osData, currSample[,4])
    }
    colnames(osData) = c("probeID", "CHROMOSOME", "START_POS", "END_POS", sampleNames)
}




main = function(runAsCohort=F) {
    sampleNames = c("1-RV", "2-AD", "3-ES", "4-GM", "5-LD",  "6-VJ",  "7-DG",  "8-MM", "9-LA", "10-CB",  "11-BG",  "12-BC",  "13-VT",  "14-CJ", "15-GG", "16-DD", "17-VV", "18-JA", "19-BF", "20-CJ", "21-DC" )
    # sampleNames = c("4-GM", "9-LA", "16-DD")
    
    ## load all-samples probeset.txt file
    allSamplesClean_path = file.path(dataDirProbesets, "allSamplesCleanProbeset_2_3Rec.txt")
    rawProbesData = read.table(allSamplesClean_path, h=T)
    colnames(rawProbesData) = c("probeID", "CHROMOSOME", "START_POS", "END_POS", sampleNames)
    
    # to run pipeline on one sample only    
    # sampleName = "2-AD"
    # s2Probes = dplyr::select(rawProbesData, c("probeID", "CHROMOSOME", "START_POS", "END_POS", sampleName))
    # callAllSamples = pipelineCGHcall(s2Probes)
    
    ### to run pipeline on all samples at the same time
    # rawProbesData = dplyr::select(rawProbesData, c("probeID", "CHROMOSOME", "START_POS", "END_POS", all_of(sampleNames)))
    tumor_prop = c(0.9,0.9,0.9,0.9,0.9,0.8,0.9,0.8,0.9,0.8,0.8,0.9,0.8,0.8,0.8,0.9,0.8,1,0.95,0.8,0.8) # samples 1 to 21
    ## Uncomment layout commands beneath to plot 3 samples within the same plot; also uncomment plot functions in pipelineCGHcall()
    # layout(matrix(c(1,2,3),nrow=3))
    if(runAsCohort) {
        callAllSamples = pipelineCGHcall(rawProbesData, tumor_prop)
    } else {
        callAllSamples = list()
        for (s in 1:length(sampleNames)) {
            currProbesData = dplyr::select(rawProbesData, c("probeID", "CHROMOSOME", "START_POS", "END_POS", sampleNames[s]))
            currCallRes = pipelineCGHcall(currProbesData, tumor_prop[s])
            callAllSamples = append(callAllSamples, currCallRes)
            print(c("callAllSamples: ", callAllSamples))
        }
        if(length(callAllSamples)==1) {
            callAllSamples = callAllSamples[[1]] ## if only one sample was treated, we don't want it to be in a list.
        }
    }
    # layout(matrix(c(1,1)))
    
    if(F) {
        #### to check that nothing went wrong
        rowsInfo = fData(callAllSamples)
        call1 = callAllSamples[,1]
        CN1 = calls(call1)
        complete1 = cbind(rowsInfo, CN1)
        complete1
    }
    ############### to convert probe table to segments table 
    ## retrieve probes information
    rowsInfo = fData(callAllSamples)
    ## retrieve call segments
    getPrbLvSegmentsFromCallObj = function(callRes) {
        ### callRes must be a cghCall object containing results of one or more samples, but can't be a list of cghCall objects
        sampleNames = colnames(callRes@assayData[["calls"]])
        rowsInfo = fData(callRes)
        CGHcall_segments = as.data.frame(calls(callRes))
        CGHcall_segments = cbind(rowsInfo, CGHcall_segments)
        colnames(CGHcall_segments) = c(colnames(rowsInfo), sampleNames)
        return(CGHcall_segments)
    }

# getPrbLvlSegments

    getPrbLvSegments = function(pipelineRes) {
        if(is.list(pipelineRes)) {
            print("result is a list of individual cghCall objects")
            currRes = pipelineRes[[1]]
            probeLevelSegments = getPrbLvSegmentsFromCallObj(currRes)
            for(i in 2:length(pipelineRes)) { ## I must make it so a list containing only 1 cghCall result can't exist
                currRes = pipelineRes[[i]]
                curr_prbLevSegs = getPrbLvSegmentsFromCallObj(currRes)
                probeLevelSegments = cbind(probeLevelSegments, curr_prbLevSegs[,length(curr_prbLevSegs)])
                # print(c("colnames(curr_prbLevSegs)[length(colnames(curr_prbLevSegs))]: ", colnames(curr_prbLevSegs)[length(colnames(curr_prbLevSegs))]))
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
    
    ssRun = getPrbLvSegments(callAllSamples)
    cohortRun = getPrbLvSegments(x) 
    
    
    
    # get segments tables
    allSegTables = getSegTables(CGHcall_segments,sampleNames,rowsInfo)
    # plot called data on all profiles
    source(file.path(working_dir, "CGHcall_functions.R"))
    source(file.path(working_dir, "rCGH_functions.R"))
    plotSegTables(allSegTables,sampleNames,resultsDir)

    # ...or one profile
    # allSegTables = getSegTables(CGHcall_segments,sampleName,rowsInfo)
    # plotSegTables(allSegTables,sampleName,resultsDir)
    
    ## initialize GI df
    
    # GI_ASCAT_df = data.frame(matrix(ncol = 4, nrow = length(sampleNames)))
    # colnames(GI_ASCAT_df) = c("GI", "nbAlter", "nbChr", "runTime")
    # rownames(GI_ASCAT_df) = sampleNames
    
    
    GI_CGHcall_df = data.frame(matrix(ncol = 4, nrow = length(sampleNames)))
    colnames(GI_CGHcall_df) = c("GI", "nbAlter", "nbChr", "runTime")
    rownames(GI_CGHcall_df) = sampleNames
    for (s in 1:length(allSegTables)) {
        print(paste0("======= sample ", sampleNames[s], " ======="))
        GI_res = calcGI_CGHcall(allSegTables[[s]])
        GI_CGHcall_df[sampleNames[s],] = append(GI_res, callAllSamples[,s]$processingTime)
        # print(GI_res)
    }
    ## saving GI table
    # GI_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/GI_all_methods"
    # write.table(GI_CGHcall_df,file.path(GI_dir, "gi_results_CGHcall.txt"),sep="\t",row.names=FALSE, quote=F)
    source(file.path(working_dir, "crossPackagesFunctions.R"))
    saveGI_ResToFile(GI_CGHcall_df, "CGHcall")

    # GI_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/GI_all_methods"
    # GIsAllResults = file.path(GI_dir,"gi_results_all_methods.txt")
    # allGIs = read.table(GIsAllResults, h=T)
    # write.table(allGIs,file.path(GI_dir, "gi_results_all_methods_addedCGHcall.txt"),sep="\t",row.names=FALSE, quote=F)
}

if (sys.nframe() == 0){
    main()
}
# to use:
# allSegTables
# GI_CGHcall_df

if (F) {
    #### Debug
    ## number of probes in total contained in segments
    for (currSample in 1:length(allSegTables)) {
        currData = allSegTables[[currSample]]
        print(paste0(sampleNames[currSample], ": ", sum(currData[5])))
    }
    s2 = allSegTables[[1]]
    s2$segSize = s2$End - s2$Start
    plot(s2$segSize)
}





























############################ To visualize the content of a CGHcall output: a list of *7* elements. see `?CGHcall` for more details
if (FALSE) {
    posteriorfin2 = rawCghResult[1]
    nclone = rawCghResult[2]
    nsamples = rawCghResult[3]
    nclass = rawCghResult[4]
    regionsprof = rawCghResult[5] # 4 cols. prof = profile = Sample the segment belongs to
    df_regions = as.data.frame(regionsprof[[1]])
    nb_segs = dim(df_regions)[1]
    params = rawCghResult$params
    cellularity = rawCghResult$cellularity
}



## divers plots pour explorer les resultats
### in assayData
# log ratio plot
cghRes_1sample = callAllSamples[,1]
# cghRes_1sample = CghResult[,1]
plot(cghRes_1sample)
df_logCN = copynumber(cghRes_1sample)
plot(df_logCN, pch=20, cex=0.001, main="0.001") # not smaller than 0.001
# extracting all data from a cghCall object's AssayData component
probaLoss = probloss(cghRes_1sample)
probaDoubleLoss = probdloss(cghRes_1sample)
probaGain = probgain(cghRes_1sample)
probaNorm = probnorm(cghRes_1sample)
plot(probaNorm, main="probaNorm") # probability of being called at level 0
plot(probaLoss, main="probaLoss") # probability of being called at level -1
plot(probaDoubleLoss, main="probaDoubleLoss") # probability of being called at level -2
plot(probaGain, main="probaGain") # probability of being called at level >0
### in featureData:
chr = chromosomes(cghRes_1sample)
bpstart = bpstart(cghRes_1sample)
bpend = bpend(cghRes_1sample)
# to plot only one chromosome with probability
plot(cghRes_1sample[chromosomes(cghRes_1sample)==1,1])
# extract column names
sample_names = sampleNames(callAllSamples)
# extract row names
probe_ids = featureNames(callAllSamples)
# get/set colnames and rownames
dimnames(callAllSamples)
# get/set nb rows and nb cols
dim(callAllSamples)
# get/set cellularity
pData(callAllSamples)
# not interesting here
varMetadata(callAllSamples)
varLabels(callAllSamples)
featureData(callAllSamples)
fvarMetadata(callAllSamples)
preproc(callAllSamples)
# df with rows id, chr, start and stop
rowsInfo = fData(callAllSamples)
# to show the chromosomes used here (sex chromosomes are included)
tail(chr)





################### to visualize what each step does to the data

if (FALSE) {
    ## plot to compare copynumber before/after maxmiss is applied
    # get log ratio of CN for first sample before and after maxmiss was applied
    s1CNBeforeMaxmiss = copynumber(ACGH_data)[,1]
    s1CNAfterMaxmiss = copynumber(preprocess(ACGH_data, maxmiss=30))[,1]
    lostVals_logicalVec = !(names(s1CNBeforeMaxmiss)%in%names(s1CNAfterMaxmiss))
    s1CNBeforeMaxmiss = as.data.frame(cbind(s1CNBeforeMaxmiss, lostVals_logicalVec))
    # plot(lostVals_logicalVec) ## these values are gathered around 2-3 chromosomes
    ## set colnames
    colnames(s1CNBeforeMaxmiss) = c("CN", "valsLost")
    # use colnames to find values lost after applying maxmiss
    # s1LostVals = dplyr::filter(s1CNBeforeMaxmiss, valsLost==1)
    s1LostVals = dplyr::filter(s1CNBeforeMaxmiss, valsLost==T) ################ try this line if valsLost==1 doesn't work
    s1LostVals[2] = rep("lostval", length(s1LostVals[1]))
    colnames(s1LostVals) = c("value", "valType")
    s1CNBeforeMaxmiss[2] = rep("normal", length(s1CNBeforeMaxmiss[1]))
    colnames(s1CNBeforeMaxmiss) = c("value", "valType")
    # readyToPlot = rbind(s1CNBeforeMaxmiss, s1LostVals)
    s1CNBeforeMaxmiss = cbind(s1CNBeforeMaxmiss, lostVals_logicalVec)
    readyToPlot = s1CNBeforeMaxmiss %>% mutate(valType = replace(valType, lostVals_logicalVec==T, "lostVal"))
    # rownames(readyToPlot)
    # Color selection
    colors <- c("#a52637", "#37a526")
    dotSize = c(5,0.01) # use this the same as color
    # s1CNBeforeMaxmiss = removePointsForQuickPlotting(s1CNBeforeMaxmiss)
    plot(1:length(s1CNBeforeMaxmiss$CN), s1CNBeforeMaxmiss$CN, pch=20, cex=0.01, col=colors[factor(s1CNBeforeMaxmiss$lostVals_logicalVec)], xlab="position sur le genome", ylab="log ratio par sonde")
    ## saving plots
    head(s1CNAfterMaxmiss)
    head(s1CNBeforeMaxmiss[,1])
    jpeg('afterPreprocess.jpg')
    plot(s1CNBeforeMaxmiss[,1], xlab = "genomic position", ylab = "log ratio", main="sample 5-LD before preprocess was applied")
    dev.off()
    jpeg('beforePreprocess.jpg')
    plot(s1CNAfterMaxmiss, xlab = "genomic position", ylab = "log ratio", main="sample 5-LD after preprocess was applied")
    dev.off()
}


if (FALSE) {
    ## Before/After normalization
    ## First, get log ratio of copy number of both states
    lrBeforeNorm = copynumber(cghdata)
    lrAfterNorm = copynumber(norm.cghdata)
    ## then plot it for the first sample
    s1Before = removePointsForQuickPlotting(as.data.frame(lrBeforeNorm[,1]))[[1]]
    s1After = removePointsForQuickPlotting(as.data.frame(lrAfterNorm[,1]))[[1]]
    
    # I used these two plots to make a gif
    plot(1:length(s1Before), s1Before, ylim=c(-1, 1.1), pch=20, cex=0.01)
    plot(1:length(s1After), s1After, ylim=c(-1, 1.1), pch=20, cex=0.01)
}



## before / after segmentation:
if (F) {
    ### base R plots
    segTable = segmented(seg.cghdata)
    logrTable = copynumber(seg.cghdata)
    png("logr.png")
    plot(logrTable, xlab="genomic position", ylab="log ratio", ylim=c(-2.5,1), main="nombre de copies pour un �chantillon")
    dev.off()
    png("seg.png")
    plot(segTable, xlab="genomic position", ylab="log ratio", ylim=c(-2.5,1), main="donn�es segment�es pour un �chantillon")
    dev.off()
}
### CGHcall objects plots
if (F) {
    png("logr.png")
    n = norm.cghdata
    sampleNames(n) = "nombre de copies pour un echantillon"
    plot(n)
    dev.off()
    
    png("seg.png")
    s = seg.cghdata
    sampleNames(s) = "donnees segmentees pour un echantillon"
    plot(s, xlab="genomic position", main="donnees segmentees pour un echantillon")
    dev.off()
}

## undoing splits: testing different values to see their impact
if (FALSE) {
    iter = c(10**-30, 10**-8, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100)
    for (i in 1:length(iter)) {
        undo.SD_ranging = 3
        clen_ranging=10 #10
        relSDlong_ranging=iter[i] #i*2
        print(c("relSDlong_ranging ", iter[i]))
        seg.cghdata <- segmentData(norm.cghdata, method="DNAcopy", undo.splits="sdundo",undo.SD=undo.SD_ranging, clen=clen_ranging, relSDlong=relSDlong_ranging)
        # 1. declare the file under which you want to save a plot
        png(paste('C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/plots/00',toString(i),'.png', sep="")) 
        # 2. run the plot command
        segTable = segmented(seg.cghdata)
        plot(segTable, xlab = "genomic position", ylab = "log ratio", ylim=c(-2.5,1), main=paste("relSDlong_ranging = ", iter[i]))
        # 3. close the process to complete saving 
        dev.off()
    }
}





#plot before and after post-seg normalization
if(F) {
    png("segdata_before_norm.png")
    plot(segmented(seg.cghdata))
    dev.off()
    png("segdata_after_norm.png")
    plot(segmented(postseg.cghdata))
    dev.off()
}




###################################################
### code chunk number 8: CGHcall.Rnw:122-123
###################################################
# plot call probability for sample 1
plot(CghResult[,1])
# plot(1)

###################################################
### code chunk number 9: CGHcall.Rnw:129-130
###################################################
# plot call probability for sample 2
plot(CghResult[,2])
plot(CghResult[,3])
# plot call data
plot(calls(CghResult)[,1], ylab="log ratio", xlab = "genomic position", main="5-LD", ylim=c(-5,5))
###################################################
### code chunk number 10: CGHcall.Rnw:139-140
###################################################
summaryPlot(CghResult)


###################################################
### code chunk number 11: CGHcall.Rnw:149-150
###################################################
frequencyPlotCalls(CghResult)



if (F) {
        

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


    ###### plot histogram with lots of breaks to see gaussians
    df_calls = calls(CghResult)
    df_logCN = copynumber(CghResult)
    segs = segmented(CghResult)
    ## pour voir la distribution de toutes les donnees (la somme de toutes les courbes de gauss)
    logratios_for_hist = c(df_logCN[,1], df_logCN[,2], df_logCN[,3], df_logCN[,4], df_logCN[,5])
    hist(logratios_for_hist, breaks=1000)
    ## pour voir la distribution des segments
    segvals_for_hist = c(segs[,1], segs[,2], segs[,3], segs[,4], segs[,5])
    hist(segvals_for_hist, breaks=1000)
    ## pour voir la distribution des segments called
    callvals_for_hist = c(df_calls[,1], df_calls[,2], df_calls[,3], df_calls[,4], df_calls[,5])
    hist(callvals_for_hist, breaks=1000)

    ## Get log ratio CN this way
    logratioNiveauMoinsUn = getProbesOfCallLevelfromOtherDf(-1, df_calls, df_logCN)
    logratioNiveauPlusUn = getProbesOfCallLevelfromOtherDf(1, df_calls, df_logCN)
    logratioNiveauZero = getProbesOfCallLevelfromOtherDf(0, df_calls, df_logCN)
    hist(logratios_for_hist, breaks=1000, main="gaussians found")
    hist(logratioNiveauZero, add=TRUE, breaks=1000, border="green")
    hist(logratioNiveauPlusUn, add=TRUE, breaks=1000, border="red")
    hist(logratioNiveauMoinsUn, add=TRUE, breaks=1000, border="blue")

    ## Get calls in order to have the plot of called segments: 
    calledNiveauMoinsUn = getProbesOfCallLevel(-1, df_calls)
    calledNiveauPlusUn = getProbesOfCallLevel(1, df_calls)
    calledNiveauZero = getProbesOfCallLevel(0, df_calls)
    # Adding one value to prevent the bars from being too wide
    calledNiveauZero = append(calledNiveauZero, 1)
    calledNiveauZero = append(calledNiveauZero, -1)
    calledNiveauMoinsUn = append(calledNiveauMoinsUn, 0)
    calledNiveauPlusUn = append(calledNiveauPlusUn, 0)
    tail(calledNiveauZero)
    tail(calledNiveauPlusUn)
    tail(calledNiveauMoinsUn)
    # plot hists
    hist(callvals_for_hist, breaks=1000)
    hist(calledNiveauZero, add=F, breaks=1000, border="green", main="segments called ")
    hist(calledNiveauPlusUn, add=T, breaks=1000, border="red")
    hist(calledNiveauMoinsUn, add=T, breaks=1000, border="blue")
    hist(c(-1, 1), add=T, breaks=1000)





    ## to generate random-based data to illustrate what a Gaussian Mixture Model does (it splits data into groups)
    m1 = data.frame(rnorm(200, -3), rnorm(200, -1), rep("minusOne", 200))
    colnames(m1) = c("vals1", "vals2", "grp")
    p1 = data.frame(rnorm(200, 3), rnorm(200, 1),rep("plusOne", 200))
    colnames(p1) = c("vals1", "vals2", "grp")
    zer = data.frame(rnorm(1000, 0), rnorm(1000, 0), rep("zero", 1000))
    colnames(zer) = c("vals1", "vals2", "grp")
    threeGroups = rbind(m1, p1, zer)
    # plot without colors. the data seem homogeneously dispersed
    plot(threeGroups$vals1, threeGroups$vals2, xlab = "values")
    # plotting density of this data. The data appears to be organized in three groups of different means
    # ggplot(threeGroups, aes(x = vals1)) + geom_density(aes(color = grp)) + theme_bw()
    ggplot(threeGroups, aes(x = vals1, fill = grp)) + geom_density(alpha=0.5) + theme_bw()
    # plot with colors. the groups are revealed
    plot(threeGroups$vals1, threeGroups$vals2, 
        pch = 19,
        col = factor(threeGroups$grp))

    ### other test ------------------------------------------
    data(movies)
    a <- density(movies$rating)
    b <- data.frame(a$x, a$y)
    ggplot(b, aes(x=a.x, y=a.y)) + geom_line()

    #### plot of density
    a<-rnorm(100,0,0.3) #component 1
    b<-rnorm(100,-1,0.5) #component 2
    c<-rnorm(100,1,0.8) #component 3
    d<-c(a,b,c) #overall data 

    df<-data.frame(d,id=as.factor(rep(c(1,2,3),each=100))) #add group id
    ggplot(df) +
        # stat_density(aes(x = d),                            position = "stack", geom = "line", show.legend = F, color = "red") +
        stat_density(aes(x = d, alpha=id), position = "stack", geom = "line", show.legend = F, color = "red") +
        stat_density(aes(x = d,  linetype = id), position = "identity", geom = "line", show.legend = F) +
        scale_alpha_manual(values=c(1,0,0)) +
        labs(x = "Segment value (log ratio)", y = "density")
        # ggtitle("Separation des donnees en trois groupes")



    sample1 = df_logCN[,1]

    ### simple density plots
    plot(density(df_logCN[,1]))
    plot(density(df_logCN[,2]))
    plot(density(df_logCN[,3]))

    hist(df_logCN[,1], breaks=200)

    ### use mixtools to visualize histogram + the gaussians
    library(mixtools)
    wait1 <- normalmixEM(df_logCN[,1], lambda = .5, mu = c(-1,0), sigma = 0.5)
    plot(wait1, density=TRUE)


    #- visualize mode in R : use abline()
    ## with density plot
    # Create the mode function.
    getmode <- function(v) {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
    }

    vals = c(round(rnorm(40, 3, 1)), round(rnorm(40, 8, 1)))
    ### using ggplot
    library(ggplot2)
    data_fvec = as.data.frame(vals)
    data_fvec %>% 
        ggplot(aes(x=vals)) +
        geom_density( fill="dodgerblue", alpha=0.5)+
        geom_vline(xintercept=2, size=1.5, color="red") +
        coord_cartesian(xlim = c(1, 6)) 
    ### using base R
    mod = getmode(vals)
    med = median(vals)
    plot(density(vals), xlab="Nombre de chaises par piece",ylab="Densite",main="La mediane (bleu) et le mode (rouge) d'une serie de valeurs")
    abline(v=mod, col="red")
    abline(v=med, col="blue")
    #### to see real values
    plot(factor(vals))
    abline(v=5.5, col="red")
    abline(v=med, col="blue")
}







