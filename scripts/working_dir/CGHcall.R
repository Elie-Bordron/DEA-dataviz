### R code from vignette source 'CGHcall.Rnw'

## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)

###################################################
### code chunk number 1: CGHcall.Rnw:50-53
###################################################
library(dplyr)
library(CGHcall)
data(Wilting)
# colnames(Wilting)= c("probeID",  "CHROMOSOME", "START_POS", "END_POS", "sample1", "sample2", "sample3", "sample4", "sample5")

Wilting_raw = Wilting
Wilting <- make_cghRaw(Wilting_raw)


###################################################
### code chunk number 2: CGHcall.Rnw:65-66
###################################################
cghdata <- preprocess(Wilting, maxmiss=30, nchrom=22)
## plot to compare copynumber before/after maxmiss is applied
testMaxmissDf <- preprocess(Wilting, maxmiss=30)
CNBeforeMaxmiss = copynumber(Wilting)
CNAfterMaxmiss = copynumber(testMaxmissDf)
s1CNBeforeMaxmiss = CNBeforeMaxmiss[,1]
s1CNAfterMaxmiss = CNAfterMaxmiss[,1]
lostVals_logicalVec = !(names(s1CNBeforeMaxmiss)%in%names(s1CNAfterMaxmiss))
s1CNBeforeMaxmiss = as.data.frame(cbind(s1CNBeforeMaxmiss, lostVals_logicalVec))
colnames(s1CNBeforeMaxmiss) = c("s1CNBeforeMaxmiss", "valsLost")
s1LostVals = dplyr::filter(s1CNBeforeMaxmiss, valsLost==T)
s1LostVals[2] = rep("lostval", length(s1LostVals[1]))
colnames(s1LostVals) = c("value", "valType")
s1CNBeforeMaxmiss[2] = rep("normal", length(s1CNBeforeMaxmiss[1]))
colnames(s1CNBeforeMaxmiss) = c("value", "valType")
# readyToPlot = rbind(s1CNBeforeMaxmiss, s1LostVals)
s1CNBeforeMaxmiss = cbind(s1CNBeforeMaxmiss, lostVals_logicalVec)
readyToPlot = s1CNBeforeMaxmiss %>% mutate(valType = replace(valType, lostVals_logicalVec==T, "lostVal"))
rownames(readyToPlot)
# Color selection
colors <- c("#66BD63", # Orange
            # "#D9EF8B", # Light green
            "#FDAE61") # Darker green
plot(1:length(readyToPlot$value), readyToPlot$value, pch=19, col=colors[factor(readyToPlot$valType)])



plot(s1CNBeforeMaxmiss[,1])
plot(s1CNAfterMaxmiss)

###################################################
### code chunk number 3: CGHcall.Rnw:75-77
###################################################
norm.cghdata <- normalize(cghdata, method="median", smoothOutliers=TRUE)
## Let's do a plot before and after this normalization.
## First, get log ratio of copy number of both states
lrBeforeNorm = copynumber(cghdata)
lrAfterNorm = copynumber(norm.cghdata)
## then plot it for the first sample
s1Before = lrBeforeNorm[,1]
s1After = lrAfterNorm[,1]
plot(1:length(s1Before), s1Before, ylim=c(-1, 1.1))
plot(1:length(s1After), s1After, ylim=c(-1, 1.1))



###################################################
### code chunk number 4: CGHcall.Rnw:89-92
###################################################
norm.cghdata <- norm.cghdata[,1:4] # To save time we will limit our analysis to the first two samples from here on.
seg.cghdata <- segmentData(norm.cghdata, method="DNAcopy", undo.splits="sdundo",undo.SD=3, clen=10, relSDlong=5)



###################################################
### code chunk number 5: CGHcall.Rnw:97-98
###################################################
postseg.cghdata <- postsegnormalize(seg.cghdata)


###################################################
### code chunk number 6: CGHcall.Rnw:106-108
###################################################
# tumor.prop <- c(0.75, 0.9) # one value per sample. proportion of contamination by healthy cells
tumor.prop <- c(0.75, 0.9, 0.6, 0.85) # one value per sample. proportion of contamination by healthy cells
rawResult <- CGHcall(postseg.cghdata,nclass=5,cellularity=tumor.prop)
posteriorfin2 = rawResult[1]
nclone = rawResult[2]
nsamples = rawResult[3]
nclass = rawResult[4]
regionsprof = rawResult[5]
df_regions = as.data.frame(regionsprof[[1]])
nb_segs = dim(df_regions)[1]
plot(1:nb_segs, df_regions$wm)
params = rawResult$params
cellularity = rawResult$cellularity

###################################################
### code chunk number 7: CGHcall.Rnw:113-114
###################################################
result <- ExpandCGHcall(rawResult,postseg.cghdata)
segs = as.data.frame(segmented(result))
head(segs)
segsAsVec = c(t(segs))
summaryAnyVec(segsAsVec)
###################################################
### code chunk number 8: CGHcall.Rnw:122-123
###################################################
plot(result[,1])


###################################################
### code chunk number 9: CGHcall.Rnw:129-130
###################################################
plot(result[,2])


###################################################
### code chunk number 10: CGHcall.Rnw:139-140
###################################################
summaryPlot(result)


###################################################
### code chunk number 11: CGHcall.Rnw:149-150
###################################################
frequencyPlotCalls(result)


## plots en vrac pour explorer les resultats
if(FALSE) {
df_logCN = copynumber(result)
segs = segmented(result)
df_calls = calls(result)
chr = chromosomes(result)
bpstart = bpstart(result)
bpend = bpend(result)
plot(df_logCN)
hist(df_logCN[,1], breaks=1000)
plot(segs[,1])
plot(df_calls[,1])

plot(bpstart)
plot(bpend)

called = result
plot(called[,1])
plot(called[chromosomes(called)==1,1])
log2ratios <- copynumber(called[,1])
sample_names = sampleNames(called)
probe_ids = featureNames(called)

# only the 22 first chromosomes are used here
tail(chr)
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
            library(dplyr)
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

summaryAnyVec = function(vector) {
    vector = as.numeric(vector)
    vecCount = vctrs::vec_count(vector)
    print(c("============ vec count: ============", vecCount))
}



df_calls = calls(result)
df_logCN = copynumber(result)
segs = segmented(result)
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
m1 = data.frame(rnorm(200, -1), rnorm(200, -1), rep("minusOne", 200))
colnames(m1) = c("vals1", "vals2", "class")
p1 = data.frame(rnorm(200, 1), rnorm(200, 1),rep("plusOne", 200))
colnames(p1) = c("vals1", "vals2", "class")
zer = data.frame(rnorm(1000, 0), rnorm(1000, 0), rep("zero", 1000))
colnames(zer) = c("vals1", "vals2", "class")
threeGroups = rbind(m1, p1, zer)
plot(threeGroups$vals1, threeGroups$vals2, 
     pch = 19,
     col = factor(threeGroups$class))





