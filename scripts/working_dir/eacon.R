require(EaCoN)
library(EaCoN)
pathToData = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/data/working_data/2-AD"
pathToATCelFile = file.path(pathToData, "2-AD_AT_(OncoScan_CNV).CEL")
pathToGCCelFile = file.path(pathToData, "2-AD_GC_(OncoScan_CNV).CEL")
# pathToGCCelFile = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/data/working_data/2-AD/2-AD_GC_(OncoScan_CNV).CEL"
print(c("pathToATCelFile: ", pathToATCelFile))
print(c("pathToGCCelFile: ", pathToGCCelFile))
outputFolder = "2-AD"

OS.Process(ATChannelCel = pathToATCelFile, GCChannelCel = pathToGCCelFile, samplename = outputFolder)

