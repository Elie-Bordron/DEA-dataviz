## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
resDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)
#loading libraries
library(dplyr)
library("GGally")





if (!interactive()) { ## main
    allGIRes = file.path(resDir, "GI_all_methods/GI_all_methods.txt")
    GI_table = read.table(allGIRes, h=T)
    GI_table[,2:5] = as.numeric(GI_table[,2:5])
    GI_table[,5] = as.numeric(GI_table[,5])
    pairs(~ GI_oncoscanR + GI_CGHcall + GI_rCGH + GI_ASCAT, data = GI_table)
    GI_sorted = sort(GI_table$GI_CGHcall)
    plot(GI_sorted)
    
    ggpairs(GI_table[, 2:5], lower=list(continuous=wrap("smooth", colour="blue")),
        diag=list(continuous=wrap("barDiag", fill="dark blue")),
        upper=list(corSize=6), axisLabels='show', cardinality_threshold=20)
    plot(as.data.frame(c(5,5,6), c(8,9,6), c(5,6,6)))
    
    # 
    # 
    # ggpairs(iris[, 1:4], lower=list(continuous="smooth", params=c(colour="blue")),
    #         diag=list(continuous="bar", params=c(colour="blue")), 
    #         upper=list(params=list(corSize=6)), axisLabels='show')
    # 
}













