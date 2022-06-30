# setwd( "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results") #ctrl maj H
# ## on teste EaCoN
# 
# library(affyio)
# pathToData = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/data/working_data/2-AD"
# pathToATCelFile = file.path(pathToData, "2-AD_AT_(OncoScan_CNV).CEL")
# pathToGCCelFile = file.path(pathToData, "2-AD_GC_(OncoScan_CNV).CEL")
# celfiles = c(pathToATCelFile, pathToGCCelFile)
# 
# ## various information about header
# celHeader = read.celfile.header(pathToATCelFile, info="full")
# 
# ## Probe intensity matrix. PM = perfect match, MM = mismatch
# PM = read.celfile.probeintensity.matrices
# 
# pmMatrix = read.celfile.probeintensity.matrices(pathToATCelFile, which="pm")
# 
# 
# # 
# library(apt.oncoscan.2.4.0)
# apt.oncoscan.process(ATChannelCel = pathToATCelFile, GCChannelCel = pathToGCCelFile, samplename = "sample5-LD", out.dir = getwd(), force.OS = "windows", apt.build = "na33.r2")
# 
# 
# 
# ############################# view segments files in R 
# dataDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/data/working_data/from_CEL"
# 
# sample11_BG_path = file.path(dataDir, "11-BG.OSCHP.segments.txt")
# sample11_BG_fl_path = file.path(dataDir, "11-BG.OSCHP.segments_clean_fl.txt")
# sample11_BG = read.table(sample11_BG_path, sep='\t', h=TRUE)
# sample11_BG_fl = read.table(sample11_BG_fl_path, sep='\t', h=TRUE)
# 
# 
# sample2_AD_path = file.path(dataDir, "2-AD.OSCHP.segments.txt")
# sample2_AD = read.table(sample2_AD_path, sep='\t', h=TRUE)
# 
# 
# ##### tips
# 
# ## to check if two objects point to the same memory slot
# tracemem(object1)==tracemem(object1)
# 
# 

GI_scripts_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
source(file.path(GI_scripts_dir, "CGHcall.R"))

sampleName="1-RV"
probeData = read.table("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle/1-RV_small.txt", h=TRUE, sep="\t")
probeData = mutate(probeData, END_POS=Position+20)
colnames(probeData) = c("probeID", "CHROMOSOME", "START_POS", sampleName, "END_POS")
probeData = dplyr::select(probeData, c("probeID", "CHROMOSOME", "START_POS", "END_POS", sampleName))
colnames(probeData)[c(2:3)] = c("ChrNum", "ChrStart")
probeData = getAbspos_probeset(probeData)
colnames(probeData)[c(2:3)] = c("CHROMOSOME", "START_POS")
colnames(probeData)[which(colnames(probeData)==sampleName)] <- "Log2Ratio"
plotSegTableForWGV_GG(NULL, probeData)



mtcars
g2 <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point(aes(color = factor(cyl)))
colors <- c("4" = "#D9717D", "6" = "#4DB6D0", "8" = "#BECA55")
g2 = g2 + scale_color_manual(values = colors)
g2
