### R code from vignette source 'CGHcall.Rnw'

###################################################
### code chunk number 1: CGHcall.Rnw:50-53
###################################################
library(CGHcall)
data(Wilting)
Wilting <- make_cghRaw(Wilting)


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
norm.cghdata <- norm.cghdata[,1:2]
seg.cghdata <- segmentData(norm.cghdata, method="DNAcopy",undo.splits="sdundo",undo.SD=3,
clen=10, relSDlong=5)


###################################################
### code chunk number 5: CGHcall.Rnw:97-98
###################################################
postseg.cghdata <- postsegnormalize(seg.cghdata)


###################################################
### code chunk number 6: CGHcall.Rnw:106-108
###################################################
tumor.prop <- c(0.75, 0.9)
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


