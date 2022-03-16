### R code from vignette source 'CGHcall.Rnw'

###################################################
### code chunk number 1: CGHcall.Rnw:50-53
###################################################
library(CGHcall)
data(Wilting)
Wilting_raw = Wilting
Wilting <- make_cghRaw(Wilting_raw)


###################################################
### code chunk number 2: CGHcall.Rnw:65-66
###################################################
cghdata <- preprocess(Wilting, maxmiss=30, nchrom=22)


###################################################
### code chunk number 3: CGHcall.Rnw:75-77
###################################################

norm.cghdata <- normalize(cghdata, method="median", smoothOutliers=TRUE)


###################################################
### code chunk number 4: CGHcall.Rnw:89-92
###################################################
# norm.cghdata <- norm.cghdata[,2:3] # To save time we will limit our analysis to the first two samples from here on.
seg.cghdata <- segmentData(norm.cghdata, method="DNAcopy", undo.splits="sdundo",undo.SD=3, clen=10, relSDlong=5)


###################################################
### code chunk number 5: CGHcall.Rnw:97-98
###################################################
postseg.cghdata <- postsegnormalize(seg.cghdata)


###################################################
### code chunk number 6: CGHcall.Rnw:106-108
###################################################
tumor.prop <- c(0.75, 0.9) # one value per sample. proportion of contamination by healthy cells
result <- CGHcall(postseg.cghdata,nclass=5,cellularity=tumor.prop)


###################################################
### code chunk number 7: CGHcall.Rnw:113-114
###################################################
result <- ExpandCGHcall(result,postseg.cghdata)


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
# frequencyPlot(result)

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

## pour voir la distribution de toutes les donnees (la somme de toutes les courbes de gauss)
logratios_for_hist = c(df_logCN[,1], df_logCN[,2], df_logCN[,3], df_logCN[,4], df_logCN[,5])
hist(logratios_for_hist, breaks=1000)
## pour voir la distribution des segments
segvals_for_hist = c(segs[,1], segs[,2], segs[,3], segs[,4], segs[,5])
hist(segvals_for_hist, breaks=1000)
## pour voir la distribution des segments called
callvals_for_hist = c(df_calls[,1], df_calls[,2], df_calls[,3], df_calls[,4], df_calls[,5])
hist(callvals_for_hist, breaks=1000)

## pour visualiser les sondes du groupe -1 du plot precedent (je parle de calling) au sein du plot logratios_for_hist
# allProbes_niveauPlusUn = c()
getProbesOfCallLevelfromOtherDf = function(lvl, df_calls, df_logCN) {
    nb_cols = dim(df_calls)[2]
    allProbes_niveauMoinsUn = c()
    for(i in 1:nb_cols) {
        #1 obtenir la liste des noms de ces sondes (et le sample auquel elles appartiennent)
            #1.1 extraire les 5 cols de df_calls.
            call_currSample = df_calls[,i]
            samplename = colnames(df_calls)[i]
            #1.2 pour chacune, donner la valeur TRUE aux cases qui ont la valeur 1. les autres cases recoivent la valeur FALSE.
            currSample__call_equals_lvl = call_currSample==lvl
            # currSample__call_equals_0 = call_currSample==0
        #2 extraire les valeurs de ces sondes A partir du dataframe logratios CN.
            #2.1 ajouter la colonne call_equals_lvl au df logratios CN.
            df_logCN = cbind(df_logCN, currSample__call_equals_lvl)
            colnames(df_logCN)[i+nb_cols] = paste(samplename, "__call_equals_lvl")
            #2.2 filtrer les lignes pour lesquelles la colonne call_equals_lvl vaut TRUE
            library(dplyr)
            df_logCN = as.data.frame(df_logCN)
            filteredOnCurrSample_callEquals1 = dplyr::filter(df_logCN, df_logCN[i+5]==1)
            #2.3 extraire la colonne de notre sample de ce sous-tableau
            CurrSample__callEquals1 = filteredOnCurrSample_callEquals1[,1]
        #3 ajouter ces sondes a la liste
        allProbes_niveauMoinsUn = append(allProbes_niveauMoinsUn, CurrSample__callEquals1)
        # allProbes_niveauPlusUn = append(allProbes_niveauPlusUn, CurrSample__callEquals1)
    }
    return(allProbes_niveauMoinsUn)
}
## Get log ratio CN this way
allProbes_niveauMoinsUn = getProbesOfCallLevelfromOtherDf(-1, df_calls, df_logCN)
allProbes_niveauPlusUn = getProbesOfCallLevelfromOtherDf(1, df_calls, df_logCN)
allProbes_niveauZero = getProbesOfCallLevelfromOtherDf(0, df_calls, df_logCN)
hist(logratios_for_hist, breaks=1000)
hist(allProbes_niveauZero, add=TRUE, breaks=1000, border="green")
hist(allProbes_niveauPlusUn, add=TRUE, breaks=1000, border="red")
hist(allProbes_niveauMoinsUn, add=TRUE, breaks=1000, border="blue")
## Get calls in order to have this plot: 
allProbes_niveauMoinsUn = getProbesOfCallLevelfromOtherDf(-1, df_calls, df_calls)
allProbes_niveauPlusUn = getProbesOfCallLevelfromOtherDf(1, df_calls, df_calls)
allProbes_niveauZero = getProbesOfCallLevelfromOtherDf(0, df_calls, df_calls)
hist(callvals_for_hist, breaks=1000)
hist(allProbes_niveauZero, add=FALSE, breaks=1000, border="green")
hist(allProbes_niveauPlusUn, add=FALSE, breaks=1000, border="red")
hist(allProbes_niveauMoinsUn, add=FALSE, breaks=1000, border="blue")

## testing and debugging
lvl = 1
call_currSample = df_calls[,2]
samplename = colnames(df_calls)[2]
#1.2 pour chacune, donner la valeur TRUE aux cases qui ont la valeur 1. les autres cases recoivent la valeur FALSE.
currSample__call_equals_lvl = call_currSample==lvl
vctrs::vec_count(call_currSample)
currSample__call_equals_lvl_num = as.numeric(currSample__call_equals_lvl)
vctrs::vec_count(currSample__call_equals_lvl_num)


#####################
#####################
#####################
## FAIRE UNE VERIF ## 
#####################
#####################
#####################






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
