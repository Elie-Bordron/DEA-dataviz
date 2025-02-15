#### Debug
if (F) {
    ## number of probes in total contained in segments
    for (currSample in 1:length(allSegTables)) {
        currData = allSegTables[[currSample]]
        print(paste0(sampleNames[currSample], ": ", sum(currData[5])))
    }
    s2 = allSegTables[[1]]
    s2$segSize = s2$End - s2$Start
    plot(s2$segSize)
}

## to view which probes are removed 
x = as.data.frame(cghdata@assayData[["copynumber"]])
x = as.data.frame(ACGH_data@assayData[["copynumber"]])
removedNaNProbes = dplyr::filter(x, !is.na(x["5-LD"]))
testDf = data.frame(size=c(5,9,3), weight=c(8,8,NULL))
dplyr::filter(testDf, testDf["weight"])

## for loading data from individual probeset.txt files
if(F) {
    print("changing value of dataDirProbesets")
    # dataDirProbesets = "C:/Users/User/Desktop/CGH-scoring/M2_internship_Bergonie/testData"
    dataDirProbesets = "C:/Users/User/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
    # sampleNames = c("5-LD", "6-VJ", "8-MM")
    sampleNames = c("1-RV")
    osData = NULL
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

## to create a subset of a probeset.txt file
if (F) {
    dataDirProbesets = "C:/Users/User/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
    sampleName = "1-RV"
    currSamplePath = paste0(dataDirProbesets, "/", sampleName, ".probeset.txt")
    currSample = read.table(currSamplePath, sep='\t', h=T)
    # source(paste0(working_dir, "rCGH_functions.R"))
    source(paste0(GI_scripts_dir, "/rCGH_functions.R"))
    currSample = removePointsForQuickPlotting(currSample, 100)
    write.table(currSample, paste0(dataDirProbesets, "/", sampleName, "_small_complete.probeset.txt"), quote=F, sep='\t', row.names=FALSE)
}

# ProbeSetName	Chromosome	Position	Log2Ratio (1-RV.OSCHP)	WeightedLog2Ratio (1-RV.OSCHP)	AllelicDifference (1-RV.OSCHP)	NormalDiploid (1-RV.OSCHP)	BAF (1-RV.OSCHP)














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
    s1CNBeforeMaxmiss = removePointsForQuickPlotting(s1CNBeforeMaxmiss)
    s1CNAfterMaxmiss = copynumber(preprocess(ACGH_data, maxmiss=30))[,1]
    plot(s1CNBeforeMaxmiss$s1CNBeforeMaxmiss, grp=s1CNBeforeMaxmiss$lostVals_logicalVec)
    
    
    
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
    # library(mixtools)
    # wait1 <- normalmixEM(df_logCN[,1], lambda = .5, mu = c(-1,0), sigma = 0.5)
    # plot(wait1, density=TRUE)


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


