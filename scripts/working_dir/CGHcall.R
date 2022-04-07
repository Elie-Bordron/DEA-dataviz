### R code from vignette source 'CGHcall.Rnw'

## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)

## setting paths
dataDir = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/all_probeset"
resultsDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results"
# paths to probeset.txt files
s5_LD_path = file.path(dataDir, "5-LD.probeset.txt")
s6_VJ_path = file.path(dataDir, "6-VJ.probeset.txt")
s8_MM_path = file.path(dataDir, "8-MM.probeset.txt")
# loading files as tables

s5_LD = read.table(s5_LD_path, sep='\t', h=T)
s6_VJ = read.table(s6_VJ_path, sep='\t', h=T)
s8_MM = read.table(s8_MM_path, sep='\t', h=T)
osData = s5_LD[,1:3]
# osData = cbind(osData, rep(NA, length(osData[,1])))
osData = cbind(osData, s5_LD[,3]+20)
osData = cbind(osData, s5_LD[,4])
osData = cbind(osData, s6_VJ[,4])
osData = cbind(osData, s8_MM[,4])
colnames(osData)= c("probeID",  "CHROMOSOME", "START_POS", "END_POS", "sample1", "sample2", "sample3")



###################################################
### code chunk number 1: CGHcall.Rnw:50-53
###################################################
library(dplyr)
library(CGHcall)
library(ggplot2)

# data(Wilting)
# to set meaningful column names to this df
# colnames(Wilting)= c("probeID",  "CHROMOSOME", "START_POS", "END_POS", "sample1", "sample2", "sample3", "sample4", "sample5")
# Wilting_raw = Wilting
# Wilting <- make_cghRaw(Wilting_raw)
Wilting <- make_cghRaw(osData)

###################################################
### code chunk number 2: CGHcall.Rnw:65-66
###################################################
cghdata <- preprocess(Wilting, maxmiss=5, nchrom=22)

if (FALSE) {
## plot to compare copynumber before/after maxmiss is applied
# get log ratio of CN for first sample before and after maxmiss was applied
s1CNBeforeMaxmiss = copynumber(Wilting)[,1]
s1CNAfterMaxmiss = copynumber(preprocess(Wilting, maxmiss=30))[,1]
lostVals_logicalVec = !(names(s1CNBeforeMaxmiss)%in%names(s1CNAfterMaxmiss))
s1CNBeforeMaxmiss = as.data.frame(cbind(s1CNBeforeMaxmiss, lostVals_logicalVec))
View(as.data.frame(lostVals_logicalVec))
# set colnames
colnames(s1CNBeforeMaxmiss) = c("s1CNBeforeMaxmiss", "valsLost")
# use colnames to find values lost after applying maxmiss
s1LostVals = dplyr::filter(s1CNBeforeMaxmiss, valsLost==T)
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
plot(1:length(readyToPlot$value), readyToPlot$value, pch=20, col=colors[factor(readyToPlot$valType)], xlab="position sur le genome", ylab="log ratio par sonde")
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
###################################################
### code chunk number 3: CGHcall.Rnw:75-77
###################################################
norm.cghdata <- normalize(cghdata, method="median", smoothOutliers=TRUE)
# plotting log ratio of copy number
# plot(norm.cghdata)

if (FALSE) {
## Let's do a plot before and after this normalization.
## First, get log ratio of copy number of both states
lrBeforeNorm = copynumber(cghdata)
lrAfterNorm = copynumber(norm.cghdata)
## then plot it for the first sample
s1Before = lrBeforeNorm[,1]
s1After = lrAfterNorm[,1]
# I used these two plots to make a gif
plot(1:length(s1Before), s1Before, ylim=c(-1, 1.1))
plot(1:length(s1After), s1After, ylim=c(-1, 1.1))
}


###################################################
### code chunk number 4: CGHcall.Rnw:89-92
###################################################
norm.cghdata <- norm.cghdata[,1] # To save time we will limit our analysis to the first two samples from here on.
seg.cghdata <- segmentData(norm.cghdata, method="DNAcopy", undo.splits="sdundo",undo.SD=3, clen=10, relSDlong=5)

## before / after segmentation:
if (F) {
### base R plots
segTable = segmented(seg.cghdata)
logrTable = copynumber(seg.cghdata)
png("logr.png")
plot(logrTable, xlab="genomic position", ylab="log ratio", ylim=c(-2.5,1), main="nombre de copies pour un ï¿½chantillon")
dev.off()
png("seg.png")
plot(segTable, xlab="genomic position", ylab="log ratio", ylim=c(-2.5,1), main="donnï¿½es segmentï¿½es pour un ï¿½chantillon")
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
        png(paste('C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/plots/00',toString(i),'.png', sep="")) # output folder: C:\Users\e.bordron\Desktop\CGH-scoring\M2_internship_Bergonie\scripts\working_dir\plots
        # 2. run the plot command
        segTable = segmented(seg.cghdata)
        plot(segTable, xlab = "genomic position", ylab = "log ratio", ylim=c(-2.5,1), main=paste("relSDlong_ranging = ", iter[i]))
        # 3. close the process to complete saving 
        dev.off()
    }
}


###################################################
### code chunk number 5: CGHcall.Rnw:97-98
###################################################

postseg.cghdata <- postsegnormalize(seg.cghdata)
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
### code chunk number 6: CGHcall.Rnw:106-108
###################################################
# tumor.prop <- c(0.75, 0.9, 0.6, 0.85, 0.65) # one value per sample. proportion of contamination by healthy cells
tumor.prop <- c(0.75, 0.9, 0.6) # one value per sample. proportion of contamination by healthy cells
rawResult <- CGHcall(postseg.cghdata,nclass=5,cellularity=tumor.prop)
## To visualize the content of a CGHcall output: a list of *7* elements. see `?CGHcall` for more details
if (FALSE) {
posteriorfin2 = rawResult[1]
nclone = rawResult[2]
nsamples = rawResult[3]
nclass = rawResult[4]
regionsprof = rawResult[5] # 4 cols. profile = Sample the segment belongs to
df_regions = as.data.frame(regionsprof[[1]])
nb_segs = dim(df_regions)[1]
params = rawResult$params
cellularity = rawResult$cellularity
}
###################################################
### code chunk number 7: CGHcall.Rnw:113-114
###################################################
result <- ExpandCGHcall(rawResult,postseg.cghdata,CellularityCorrectSeg=F) # use CellularityCorrectSeg=TRUE to correct using cellularity

###################################################
### code chunk number 8: CGHcall.Rnw:122-123
###################################################
# plot call probability for sample 1
plot(result[,1])
plot(1)
# plot raw call values
df_calls = calls(result)
plot(df_calls[,1], )

###################################################
### code chunk number 9: CGHcall.Rnw:129-130
###################################################
# plot call probability for sample 2
plot(result[,2])

###################################################
### code chunk number 10: CGHcall.Rnw:139-140
###################################################
summaryPlot(result)


###################################################
### code chunk number 11: CGHcall.Rnw:149-150
###################################################
frequencyPlotCalls(result)


## divers plots pour explorer les resultats
if(FALSE) {
# in assayData
df_logCN = copynumber(result)
plot(df_logCN)
segs = segmented(result)

plot(segs)
png("segs.png")
plot(segs, ylab="log ratio", xlab = "position genomique", main="données de segmentation")
dev.off()

png("calls.png")
plot(calls(result), ylab="log ratio", xlab = "position genomique", main="segments appelés à leurs niveaux estimés de log ratio")
dev.off()
probaLoss = probloss(result)
probaDoubleLoss = probdloss(result)
probaGain = probgain(result)
probaNorm = probnorm(result)
plot(probaNorm)
plot(probaLoss)
plot(probaDoubleLoss)
plot(probaGain)
# in featureData:
chr = chromosomes(result)
bpstart = bpstart(result)
bpend = bpend(result)
# plots
called = result
res_sample1 = called[,1]
plot(called[,1],ampcol="#00ffff",dlcol="#00ff00")
plot(res_sample1)
plot(segs[,1])
plot(df_calls[,1])
plot(bpend)

plot(result)
plot(called[chromosomes(called)==1,1])
log2ratios <- copynumber(res_sample1)

sample_names = sampleNames(called)
probe_ids = featureNames(called)
# get/set colnames and rownames
dimnames(called)
# get/set nb rows and nb cols
dim(called)
# get/set cellularity
pData(called)
# not interesting here
varMetadata(called)
varLabels(called)
featureData(called)
fvarMetadata(called)
preproc(called)
# df with rows id, chr, start and stop
fData(called)
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


# testing cluster circling
library(tidyverse)
library(ggforce)
library(palmerpenguins)
#removing NAs
penguins <- penguins %>%
    drop_na()
#plotting ellipses
threeGroups %>%
    ggplot(aes(x = threeGroups$vals1,
               y = threeGroups$vals2))+
    geom_mark_ellipse(aes(color = threeGroups$grp,
                          label=threeGroups$grp),
                      expand = unit(0.5,"mm"),
                      label.buffer = unit(-5, 'mm'))+
    geom_point(aes(color=threeGroups$grp))+
    theme(legend.position = "none")
ggsave("annotate_groups_clusters_with_ellipse_ggplot2.png")





#### plot of density
a<-rnorm(1000,0,0.3) #component 1
b<-rnorm(1000,-1,0.5) #component 2
c<-rnorm(1000,1,0.8) #component 3
d<-c(a,b,c) #overall data 
df<-data.frame(d,id=rep(c(1,2,3),each=1000)) #add group id
ggplot(df) +
    stat_density(aes(x = d,  linetype = as.factor(id)), position = "stack", geom = "line", show.legend = F, color = "red") +
    stat_density(aes(x = d,  linetype = as.factor(id)), position = "identity", geom = "line", show.legend = F) +
    labs(x = "value", y = "density")
    # ggtitle("")

sample1 = df_logCN[,1]
# ggplot(df_logCN) + stat_density(aes(x = sample1, geom = "line", show.legend = F, color = "red")

### simple density plots
plot(density(df_logCN[,1]))
plot(density(df_logCN[,2]))
plot(density(df_logCN[,3]))

hist(df_logCN[,1], breaks=200)

### use mixtools to visualize histogram + the gaussians
library(mixtools)
wait1 <- normalmixEM(df_logCN[,1], lambda = .5, mu = c(-1,0), sigma = 0.5)
plot(wait1, density=TRUE)











####
library(mclust)
c(rnorm(200, -1),rnorm(200, 1),rnorm(1000, 0))

data(diabetes)
class <- diabetes$class
table(class)
X <- diabetes[,-1]
head(X)
clPairs(X, class)
####
# library(mixtools)
# wait = faithful$waiting
# mixmdl = normalmixEM(wait)
# plot(mixmdl,which=2)
# lines(density(wait), lty=2, lwd=2)




#- visualize mode in R : use abline()
## with density plot
# Create the mode function.
getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

### using ggplot
library(ggplot2)
data_fvec = as.data.frame(v)
data_fvec %>% 
    ggplot(aes(x=v)) +
    geom_density( fill="dodgerblue", alpha=0.5)+
    geom_vline(xintercept=2, size=1.5, color="red") +
    coord_cartesian(xlim = c(1, 6)) 
### using base R
vals = c(round(rnorm(40, 3, 1)), round(rnorm(40, 8, 1)))
mod = getmode(vals)
med = median(vals)
plot(density(vals), xlab="Nombre de chaises par pièce",ylab="Densité",main="La médiane (bleu) et le mode (rouge) d'une série de valeurs")
abline(v=mod, col="red")
abline(v=med, col="blue")
#### to see real values
plot(factor(vals))
abline(v=5.5, col="red")
abline(v=med, col="blue")














## to understand postsegnormalize

customPostsegnormalize = function(segmentData, inter = c(-0.1, 0.1)) 
{
    seg <- segmented(segmentData)
    values <- c()
    for (i in 1:ncol(seg)) {
        values <- c(values, median(seg[, i]))
    }
    matrixValues <- matrix(rep(values, nrow(seg)), ncol = ncol(seg), 
                           byrow = TRUE)
    seg <- seg - matrixValues
    countlevall <- apply(seg, 2, function(x) {
        as.data.frame(table(x))
    })
    print(c("countlevall aka segvec: ", countlevall))
    intcount <- function(int, sv) {
        # print("OO======================intcount=======================OO")
        # print(c("int: ", int))
        # print(c("sv: ", sv)) #table with two columns: the probe log ratio value, and the number of times it appears (== 1 segment)
        sv1 <- as.numeric(as.vector(sv[, 1])) #sv1: log ratio column
        wh <- which(sv1 <= int[2] & sv1 >= int[1]) # position of probes which value is within the interval
        # print(c("sv1: ", sv1))
        # print(c("wh: ", wh)) 
        # print(c("to_sum: ", (sv[wh, 2])))
        # print(c("returned_sum: ", sum(sv[wh, 2])))
        return(sum(sv[wh, 2])) # the number of probes(not their value) within the interval is returned
    }
    postsegnorm <- function(segvec, int = inter, intnr = 3) { #int = c(-0.5, 0.3)
        # print("OO=============================================postsegnorm================================================OO")
        # print(c("interval received: ", int))
        intlength <- (int[2] - int[1])/2
        gri <- intlength/intnr
        intst <- int[1] + (0:intnr) * gri # e.g. (-0.50000000 -0.23333333  0.03333333  0.30000000)
        intend <- intst + intlength # e.g. (0.3000000 0.5666667 0.8333333 1.1000000)
        ints <- cbind(intst, intend) # intervals
        # print(c("intst: ", intst))
        # print(c("intend: ", intend))
        # print(c("intervals list 1/2: ", ints))
        intct <- apply(ints, 1, intcount, sv = segvec) # for each interval, finds the probes contained in it then counts them. we often end up with around 3 logR values repeated each dozens/hundreds of times. The total number of probes is then counted for this interval.
        # intct contains then one value per interval, which represents how much segmented the data is.
        whmax <- which.max(intct) # finds highest value, representing the best interval.
        # print(c("segvec: ", segvec))
        # print(c("interval count: ", intct))# of all 
        # print(c("whmax: ", whmax))
        # print(c("ints: ", ints))
        print(c("best interval found (ints[whmax, ]): ", ints[whmax, ]))
        return(ints[whmax, ]) #returns said best interval.
    }
    postsegnorm_rec <- function(segvec, int, intnr = 3) {
        print("============== postsegnorm_rec =============")
        newint <- postsegnorm(segvec, int, intnr)
        newint <- postsegnorm(segvec, newint, intnr)
        newint <- postsegnorm(segvec, newint, intnr)
        newint <- postsegnorm(segvec, newint, intnr)
        newint <- postsegnorm(segvec, newint, intnr)
        return(newint[1] + (newint[2] - newint[1])/2) # we return the middle point of the interval.
    }
    listres <- lapply(countlevall, postsegnorm_rec, int = inter) # this runs postsegnorm_rec once for every sample in our seg dataset.
    vecres <- c()
    for (i in 1:length(listres)) {
        vecres <- c(vecres, listres[[i]])
    }
    print(c("listres: ", listres))
    print(c("vecres: ", vecres))
    segmented(segmentData) <- t(t(seg) - vecres) # this substracts vecres[i](the middle point of the interval found) to all segmentation data, hence normalizing.
    copynumber(segmentData) <- t(t(copynumber(segmentData) - 
                                       matrixValues) - vecres)
    return(segmentData)
}


#function call
postseg.cghdata <- customPostsegnormalize(seg.cghdata) # argument: objet cghSeg
# plot before/after postsegnormalize
png("./plots/before_postsegnorm.png")
sampleNames(seg.cghdata) = "segmented data"
plot(seg.cghdata, ylim = c(-2,2), ylab="log ratio", xlab="position genomique", main="segmented data")
dev.off()
png("./plots/after_postsegnorm.png")
sampleNames(postseg.cghdata) = "segmented data after normalization"
plot(postseg.cghdata, ylim = c(-2,2), ylab="log ratio", xlab="position genomique", main="segmented data after normalization")
dev.off()

segTable = segmented(seg.cghdata)
# function for plotting 
# interval=65654654
# plot(segTable, ylim=plot_ylim, main="recherche du meilleur intervalle", ylab="log ratio", xlab="position genomique")
plotInterval = function(data, plot_ylim, interval) {
    png(paste("./plots/postsegnorm_recherche_intervalle_", toString(interval), ".png", sep = ""))
    try(plot(data, ylim=plot_ylim, main="recherche du meilleur intervalle", ylab="log ratio", xlab="position genomique"))
    abline(h=interval[1], col="#3b2dbc"); abline(h=interval[2], col="#3b2dbc")
    abline(h=plot_ylim[1], col="#3b2dbc"); abline(h=plot_ylim[2], col="#3b2dbc") #59, 45, 188
    nb_points = length(data)
    x = -10000:nb_points+10000
    polygon(x=c(1, nb_points, nb_points, 1), y=c(interval[1],interval[1],interval[2],interval[2]), col = "#65BFFF", density = 10)
    dev.off()
}
prevInter = c(-2.5,1) # c(-2.5,1) is default ylim for this data
currInter = c(-0.1,0.1)
plotInterval(segTable,prevInter,currInter) 
prevInter = currInter
currInter = c(-0.0333,0.0667) #check
plotInterval(segTable,prevInter,currInter) 
prevInter = currInter
currInter = c(-0.0167,0.0333) #check
plotInterval(segTable,prevInter,currInter) 
prevInter = currInter
currInter = c(-0.0083,0.0167)
plotInterval(segTable,prevInter,currInter) 
prevInter = currInter
currInter = c(-0.0083,0.0042)
plotInterval(segTable,prevInter,currInter) 
prevInter = currInter
currInter = c(-0.0042,0.0021)
plotInterval(segTable,prevInter,currInter) 









intnr=3
int = c(-0.5, 0.3)
intlength = 0.8
gri= intlength/intnr
intst = -0.5 + (0:intnr) * gri
intend = intst + intlength
print(intst)
print(intend)


x1 = -0.033333
x2 = 0.666667

x1 + (x2-x1)/2


